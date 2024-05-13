%занчени€ даны
i1_1=0
i2_1=0
i3_1=0

i1_2=1
i2_2=0
i3_2=0

w2=100
k2=0.5

% указать самому из своих соображений

% n, m максимальные
N_Max = 0
M_Max = 0

%количество вызовов сумм
%size = sum(min((0:N_Max)', M_Max) + 1)

%шаг сетки
h_A  = 1e-3
% начальные шаги поиска пары собственных значений
delta_o= 100000
delta_t = 100000
%“очность поиска собственных значений
e = 1e-5
%точность метода –унге- утты 
options_h_l = odeset('RelTol', 1e-5, 'AbsTol', 1e-7, 'Refine', 10)
options_A = odeset('RelTol', 1e-5, 'AbsTol', 1e-7, 'Refine', 10)

%начало вычислений амплитуды рассе€ни€
w=sqrt(w2)
k=sqrt(k2)

n_theta = 1/h_A + 1
n_fi = int32((k^(-2)-1)/h_A + 1)
term_1 = zeros(n_theta, n_fi)
term_2 = zeros(n_theta, n_fi)   

for n = 0:N_Max
    e_1 = exp(2i*(-pi*(i1_1 + i2_1 + (1 - i3_1))/2 - pi*n));
    e_2 = exp(2i*(-pi*(i1_2 + i2_2 + (1 - i3_2))/2 - pi*n));
    
    M_N = min(n, M_Max);
    
    for m = 0:M_N
        [h_1, l_1] = find_h_l(i1_1, i2_1, i3_1, w2, k2, w, k, m, n, delta_o, delta_t, e, options_h_l);
        [Xi_1_1, Xi_2_1] = find_A(h_1, l_1, i1_1, i2_1, i3_1, w2, k2, w, k, m, n, options_A, h_A);
        Psi_1 = Xi_1_1(n_theta) * Xi_2_1(n_fi);
        
        [h_2, l_2] = find_h_l(i1_1, i2_1, i3_1, w2, k2, w, k, m, n, delta_o, delta_t, e, options_h_l);
        [Xi_1_2, Xi_2_2] = find_A(h_2, l_2, i1_2, i2_2, i3_2, w2, k2, w, k, m, n, options_A, h_A);
        Psi_2 = Xi_1_2(n_theta) * Xi_2_2(n_fi);
        
        parfor q = 1:n_theta
            for j = 1:n_fi
                term_1(q, j) = term_1(q, j) + (e_1 * Xi_1_1(q) * Xi_2_1(j) * Psi_1);
                term_2(q, j) = term_2(q, j) + (e_2 * Xi_1_2(q) * Xi_2_2(j) * Psi_2);
            end
        end
    end
end

F = zeros(n_theta, n_fi)
parfor q = 1 : n_theta
    for j = 1 : n_fi
        F(q, j) = abs(16 * pi * i * (term_1(q, j)+term_2(q, j)))
    end
end
 
X_1 = 0:1:n_theta-1
theta_ = 360/(n_theta-1)*X_1
X_2 = 0:1:double(n_fi)-1
fi_ = 360/double(n_fi-1)*X_2

theta = repmat(theta_', 1, n_fi);
fi = repmat(fi_, n_theta, 1);

%перевести в декартовыв
X = F .* sin(theta) .* cos(fi)
Y = F .* sin(theta) .* sin(fi)
Z = F .* cos(theta)


save Res_Data.mat

surf(X, Y, Z)
%plot3(X, Y, Z)

%встроеный переводчик в дкартовы от матлаба, переводит по другому
%[x,y,z] = sph2cart(deg2rad(theta), deg2rad(fi), F);
%surf(x,y,z)

