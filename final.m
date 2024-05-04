%занчения даны
i1=0
i2=1
i3=0
w2=100
k2=0.5
m=50
n=105

% указать самому из своих соображений
delta_o=10000
delta_t =5000
e = 1e-4
options_h_l = odeset('RelTol' , 1e-7, 'AbsTol' , 1e-10 ,'Refine','8', 'MaxStep',1e-2)
options_A = odeset('RelTol' , 1e-7, 'AbsTol' , 1e-10 ,'Refine','8', 'MaxStep',1e-2)
h_A  = 1e-3
    
[h, l] = find_h_l(i1, i2,i3,w2, k2, m, n, delta_o, delta_t, e, options_h_l)
res = find_A(h, l, i1, i2,i3,w2, k2, m, n, options_A, h_A)
 


