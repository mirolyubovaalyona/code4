%вспомогательная фрзовая функция $3

i1=1
i2=0
i3=1
k2=0.9
w2=1
k=sqrt(k2)
w=sqrt(w2)
m=9
n=10
h=465.062
l=456.820



t=l
o=h-l

%t =  456.777
%o = 8.2459

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

th_1_L=ode45(@(x, th) ODE_th_1_L( x, th, w, k, o, t), [0, 1/2], th0_1)
th_1_R=ode45(@(x, th) ODE_th_1_R( x, th, w, k, o, t), [1, 1/2], th0_2)
th_2_L=ode45(@(x, th) ODE_th_2_L( x, th, w, k, o, t), [1,(1+k.^(-2))/2], th0_3)
th_2_R=ode45(@(x, th) ODE_th_2_R( x, th, w, k, o, t), [k.^(-2), (1+k.^(-2))/2], th0_4)


th_1_L_=deval(th_1_L,c1)
th_1_R_=deval(th_1_R,c1)
th_2_L_=deval(th_2_L,c2)
th_2_R_=deval(th_2_R,c2)

th_1_L_>th_1_R_
th_2_L_>th_2_R_


figure
plot(th_1_L.x,th_1_L.y)
hold on
plot(th_1_R.x,th_1_R.y)
title('[0, 1]')
hold off;

figure
plot(th_2_L.x,th_2_L.y)
hold on
plot(th_2_R.x,th_2_R.y)
title('[1, k^2]')
hold off;

%%%%%%%
figure
plot(th_1_L.x,th_1_L.y)
hold on
plot(th_1_R.x,th_1_R.y)
plot(th_2_L.x,th_2_L.y)
plot(th_2_R.x,th_2_R.y)
title('[0, k^2]')
hold off;