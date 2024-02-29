
delta_o=50
delta_t =50
e =0.00001


i1=1
i2=0
i3=1
k2=0.9
w2=1
k=sqrt(k2)
w=sqrt(w2)
m=9
n=10



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


% параграф 4

o = 8.2459
t = 456.7771 



%================_____th_____=====================
th_1_L=ode45(@(x, th) ODE_th_1_L( x, th, w, k, o, t), [0, 1/2], th0_1)
th_1_R=ode45(@(x, th) ODE_th_1_R( x, th, w, k, o, t), [1, 1/2], th0_2)
th_2_L=ode23tb(@(x, th) ODE_th_2_L( x, th, w, k, o, t), [1,(1+k.^(-2))/2], th0_3)
th_2_R=ode45(@(x, th) ODE_th_2_R( x, th, w, k, o, t), [k.^(-2), (1+k.^(-2))/2], th0_4)


th_1_L=deval(th_1_L,c1)
th_1_R=deval(th_1_R,c1)
th_2_L=deval(th_2_L,c2)
th_2_R=deval(th_2_R,c2)


 %================_____fi_0______=====================
fi0_1=(th_1_L+pi/2+th_1_R+pi/2)/2
fi0_2=(th_2_L+pi/2+th_2_R+pi/2)/2


 %================_____fi_____=====================
fi_1_L=ode45(@(x, fi) ODE_fi_1_L( x, fi, w, k, o, t), [1/2, 0], fi0_1)
fi_1_R=ode45(@(x, fi) ODE_fi_1_R( x, fi, w, k, o, t), [1/2, 1], fi0_1)
fi_2_L=ode45(@(x, fi) ODE_fi_2_L( x, fi, w, k, o, t), [(1+k.^(-2))/2, 1], fi0_2)
fi_2_R=ode45(@(x, fi) ODE_fi_2_R( x, fi, w, k, o, t), [(1+k.^(-2))/2, k.^(-2)], fi0_2)


fi_1_L=deval(fi_1_L,c1)
fi_1_R=deval(fi_1_R,c1)
fi_2_L=deval(fi_2_L,c2)
fi_2_R=deval(fi_2_R,c2)

sin(th_1_L-fi_1_L)
sin(th_1_R-fi_1_R)
sin(th_2_L-fi_2_L)
sin(th_2_R-fi_2_R)
