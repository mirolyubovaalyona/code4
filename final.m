%�������� ����
i1_1=0
i2_1=0
i3_1=0

i1_2=1
i2_2=0
i3_2=0

w2=100
k2=0.5

% ������� ������ �� ����� �����������

% n, m ������������
N_Max = 0
M_Max = 0
%��� �����
h_A  = 1e-2
% ��������� ���� ������ ���� ����������� ��������
delta_o=3000
delta_t =2000
%�������� ������ ����������� ��������
e = 1e-5
%�������� ������ �����-����� 
options_h_l = odeset('RelTol', 1e-5, 'AbsTol', 1e-7, 'MaxStep',1e-2)
options_A = odeset('RelTol', 1e-4, 'AbsTol', 1e-5, 'MaxStep',1e-2)

%������ ���������� ��������� ���������
w=sqrt(w2)
k=sqrt(k2)
n_theta = 1/h_A + 1
n_fi = int8((k^(-2)-1)/h_A + 1)
%n_fi=2
term_1 = zeros(n_theta, n_fi)
term_2 = zeros(n_theta, n_fi)
for n = 0 : N_Max
    for m = 0 : M_Max
        e_1 = exp(i*(-pi*(i1_1+i2_1+(1-i3_1))/2-pi*n))
        [h_1, l_1] = find_h_l(i1_1, i2_1, i3_1, w2, k2, w, k, m, n, delta_o, delta_t, e, options_h_l)
        [X_theta_1, Xi_1_1, X_fi_1, Xi_2_1] = find_A(h_1, l_1, i1_1, i2_1, i3_1, w2, k2, w, k,  m, n, options_A, h_A)
        w_xi3_1_placeholder = 7
        A_xi3_1_placeholder = 9
        Psi_1 = Xi_1(n_theta) * Xi_2(n_fi)
        
        e_2 = exp(i*(-pi*(i1_2+i2_2+(1-i3_2))/2-pi*n))
        [h_2, l_2] = find_h_l(i1_2, i2_2, i3_2, w2, k2, w, k, m, n, delta_o, delta_t, e, options_h_l)
        [X_theta_2, Xi_1_2, X_fi_2, Xi_2_2] = find_A(h_2, l_2, i1_2, i2_2, i3_2, w2, k2, w, k,  m, n, options_A, h_A)
        w_xi3_2_placeholder = 7
        A_xi3_2_placeholder = 9
        Psi_2 = Xi_1(n_theta) * Xi_2(n_fi)
        
        for i = 1 : n_theta
            for j = 1 : n_fi
                term_1 = e_1 * Xi_1_1(i)*Xi_2_1(j) * Psi_1 * w_xi3_1_placeholder * A_xi3_1_placeholder
                term_2 = e_2 * Xi_1_2(i)*Xi_2_2(j) * Psi_2 * w_xi3_2_placeholder * A_xi3_2_placeholder
            end
        end
    end
end

F = zeros(n_theta, n_fi)
save Res_Data.mat