% указать самому из своих соображений
delta_o = 300
delta_t = 300
e = 0.00001
 
%занчения даны
i1=0
i2=0
i3=0
w2=1000
k2=0.5
w=sqrt(w2)
k=sqrt(k2)
m=10
n=20





[h, l]=initial_h_l(n, m, i1, i2, i3, w, k)

t=l
o=h-l



if i1==1
    th0_1=0
else
    th0_1=pi/2
end

if i2==1
    th0_2= pi+m*pi
    th0_3=0
else
    th0_2= pi/2+m*pi
    th0_3= pi/2
end

if i3==1
    th0_4=pi+(n-m)*pi
else
    th0_4=pi/2+(n-m)*pi
end


c1=1/2
c2=(1+k.^(-2))/2
[o, t]=find_o_t(th0_1, th0_2, th0_3, th0_4 ,o, t, c1, c2, w, k,  delta_o, delta_t, e)

l_=t
h_=o+t


%%%%%%

th_1_L=ode45(@(x, th) ODE_th_1_L( x, th, w, k, o, t), [0, 1/2], th0_1)
th_1_R=ode45(@(x, th) ODE_th_1_R( x, th, w, k, o, t), [1, 1/2], th0_2)
th_2_L=ode45(@(x, th) ODE_th_2_L( x, th, w, k, o, t), [1,(1+k.^(-2))/2], th0_3)
th_2_R=ode45(@(x, th) ODE_th_2_R( x, th, w, k, o, t), [k.^(-2), (1+k.^(-2))/2], th0_4)

figure
plot(th_1_L.x,th_1_L.y)
hold on
plot(th_1_R.x,th_1_R.y)
plot(th_2_L.x,th_2_L.y)
plot(th_2_R.x,th_2_R.y)
title('th')
hold off;


l
l_

l-l_

h
h_
h-h_

th_1_L=deval(th_1_L,c1)
th_1_R=deval(th_1_R,c1)
th_2_L=deval(th_2_L,c2)
th_2_R=deval(th_2_R,c2)
