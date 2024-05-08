%занчени€ даны
i1=0
i2=0
i3=0
w2=100
k2=0.5

% указать самому из своих соображений

% n, m максимальные
N_Max = 0
M_Max = 0
%шаг сетки
h_A  = 1e-2
% начальные шаги поиска пары собственных значений
delta_o=3000
delta_t =2000
%“очность поиска собственных значений
e = 1e-5
%точность метода –унге- утты 
options_h_l = odeset('RelTol', 1e-5, 'AbsTol', 1e-7, 'MaxStep',1e-2)
options_A = odeset('RelTol', 1e-4, 'AbsTol', 1e-5, 'MaxStep',1e-2)

%начало вычислений амплитуды рассе€ни€
w=sqrt(w2)
k=sqrt(k2)
n_theta = 1/h_A + 1
n_fi = int8((k^(-2)-1)/h_A + 1)
%n_fi=2
F = zeros(n_theta, n_fi)
for n = 0 : N_Max
    for m = 0 : M_Max
        [h, l] = find_h_l(i1, i2, i3, w2, k2, w, k, m, n, delta_o, delta_t, e, options_h_l)
        [X_1, Xi_1, X_2, Xi_2] = find_A(h, l, i1, i2, i3, w2, k2, w, k,  m, n, options_A, h_A)
        for i = 1 : n_theta
            for j = 1 : n_fi
                F(i, j) = F(i,j) + Xi_1(i)+Xi_2(j) 
            end
        end
    end
end

save Res_Data.mat