tic

i1_1=0;
i2_1=0;
i3_1=0;

i1_2=1;
i2_2=0;
i3_2=0;

w2=10;

ro2=1.1;
k2=1/ro2;

thickness =  0.1; %Толщина эллипса (xi_3 = sqr_ro + thickness) 
ellipsoid_condition = 0; %0 - жёсткий |1 - мягкий 

%================указать самому из своих соображений================================
% n, m максимальные
N_min = 0;
M_min = 0;

N_Max = 10;
M_Max = 10;
n_A = 200; %кол-во шагов сетки  

% начальные шаги поиска пары собственных значений
delta_o= 5000;
delta_t = 5000;
%Точность поиска собственных значений
e = 1e-3;
%точность метода Рунге-Кутты 
options_h_l = odeset('RelTol', 1e-8, 'AbsTol', 1e-10,'Refine', '8');
options_A = odeset('RelTol', 1e-8, 'AbsTol', 1e-10, 'Refine', '8');
options_W = odeset('RelTol', 1e-5, 'AbsTol', 1e-7,'Refine', '4');


sqr_ro = 1/k2; %ro^2
xi_3 = sqr_ro + thickness; %Толщина эллипса.
h_A = 1/n_A;
h_A_1 = 1/(n_A);
h_A_2 = (sqr_ro-1)/(2*n_A);
w=sqrt(w2);
k=sqrt(k2);
i_ = 1i;
n_theta = n_A + 1;
n_fi = 2*n_A + 1;
term_1 = zeros(n_theta, n_fi);
term_2 = zeros(n_theta, n_fi);

for n = N_min:N_Max
    e_1 = exp(2*i_*(-pi*(i1_1 + i2_1 + (1 - i3_1))/2 - pi*n));
    e_2 = exp(2*i_*(-pi*(i1_2 + i2_2 + (1 - i3_2))/2 - pi*n));
    
    M_N = min(n, M_Max);
    
    for m = M_min:M_N
        Xi_1 = cell(2, 1);
        Xi_2 = cell(2, 1);
        Psi = zeros(2, 1);
        w_xi3 = zeros(2, 1);

       parfor q = 1:2
            if q == 1
                [h_1, l_1] = find_h_l(i1_1, i2_1,i3_1,w2, k2, w, k, m, n, delta_o, delta_t, e, options_h_l);
                [Xi_1{q}, Xi_2{q}] = find_A(h_1, l_1, i1_1, i2_1, i3_1, w2, k2, w, k, m, n, options_A, h_A_1, h_A_2);
                w_xi3(q) = RWEF_13_v2(sqr_ro,w2,l_1,h_1,i3_1,ellipsoid_condition,xi_3,options_W);
                Psi(q) = Xi_1{q}(n_theta) * Xi_2{q}(n_fi);
            end
            if q == 2
                [h_2, l_2] = find_h_l(i1_2, i2_2,i3_2,w2, k2, w, k, m, n, delta_o, delta_t, e, options_h_l);
                [Xi_1{q}, Xi_2{q}] = find_A(h_2, l_2, i1_2, i2_2, i3_2, w2, k2, w, k, m, n, options_A, h_A_1, h_A_2);
                w_xi3(q) = RWEF_13_v2(sqr_ro,w2,l_2,h_2,i3_2,ellipsoid_condition,xi_3,options_W);
                Psi(q) = Xi_1{q}(n_theta) * Xi_2{q}(n_fi);
            end
       end
        
        parfor q = 1:n_theta
            for j = 1:n_fi
                term_1(q, j) = term_1(q, j) + (e_1 * Xi_1{1}(q) * Xi_2{1}(j) * Psi(1) * w_xi3(1));
                term_2(q, j) = term_2(q, j) + (e_2 * Xi_1{2}(q) * Xi_2{2}(j) * Psi(2) * w_xi3(2));
            end
        end
    end
end

F = zeros(n_theta, n_fi);
parfor q = 1 : n_theta
    for j = 1 : n_fi
        F(q, j) = abs(16 * pi * i_ * (term_1(q, j)+term_2(q, j)));
    end
end
 

X_1 = 0:1:n_theta-1;
theta_ = (pi/2/(n_theta-1)*X_1);
X_2 = 0:1:double(n_fi)-1;
fi_ = (pi/double(n_fi-1)*X_2);

theta = repmat(theta_', 1, n_fi);
fi = repmat(fi_, n_theta, 1);

X = F .* sin(theta) .* cos(fi);
Y = F .* sin(theta) .* sin(fi);
Z = F .* cos(theta);

figure
hold on;

mesh(X, Y, -Z);
mesh(X, Y, Z);
mesh(X, -Y, Z);
mesh(X, -Y, -Z);

xlabel('x') 
ylabel('y')
zlabel('z')
title('Падение по X')

grid;
axis square;
hidden on;

hold off;

elapsed_time = toc;
disp(['Время выполнения программы: ' num2str(elapsed_time) ' секунд']);

save('onX_')
%size = sum(min((0:N_Max)', M_Max) + 1)