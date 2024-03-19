i1=1
i2=0
i3=1
w2=1
k2=0.9
w=sqrt(w2)
k=sqrt(k2)
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
h = 465.0230
l =  456.7771
t = l
o = h-l




%================_____th_____=====================
th_1_L=ode45(@(x, th) ODE_th_1_L( x, th, w, k, o, t), [0, 1/2], th0_1)
th_1_R=ode45(@(x, th) ODE_th_1_R( x, th, w, k, o, t), [1, 1/2], th0_2)
th_2_L=ode23tb(@(x, th) ODE_th_2_L( x, th, w, k, o, t), [1,(1+k.^(-2))/2], th0_3)
th_2_R=ode45(@(x, th) ODE_th_2_R( x, th, w, k, o, t), [k.^(-2), (1+k.^(-2))/2], th0_4)

figure
plot(th_1_L.x,th_1_L.y)
hold on
plot(th_1_R.x,th_1_R.y)
plot(th_2_L.x,th_2_L.y)
plot(th_2_R.x,th_2_R.y)
title('th')
hold off;


 %================_____fi_0______=====================
fi0_1=deval(th_1_L,c1)+pi/2
fi0_2=deval(th_1_R,c1)+pi/2
fi0_3=deval(th_2_L,c2)+pi/2
fi0_4=deval(th_2_R,c2)+pi/2

 %================_____fi_____=====================
fi_1_L=ode45(@(x, fi) ODE_fi_1_L( x, fi, w, k, o, t), [1/2, 0], fi0_1)
fi_1_R=ode45(@(x, fi) ODE_fi_1_R( x, fi, w, k, o, t), [1/2, 1], fi0_1)
fi_2_L=ode45(@(x, fi) ODE_fi_2_L( x, fi, w, k, o, t), [(1+k.^(-2))/2, 1], fi0_2)
fi_2_R=ode45(@(x, fi) ODE_fi_2_R( x, fi, w, k, o, t), [(1+k.^(-2))/2, k.^(-2)], fi0_2)


%================_____H_____=====================
H_1_L=ode45(@(x, H) ODE_H_1_L( x, H, w, k, o, t, th_1_L), [0, 1/2], 0)
H_1_R=ode45(@(x, H) ODE_H_1_R( x, H, w, k, o, t, th_1_R), [1, 1/2], 0)
H_2_L=ode45(@(x, H) ODE_H_2_L( x, H, w, k, o, t, th_2_L), [1, (1+k.^(-2))/2], 0)
H_2_R=ode45(@(x, H) ODE_H_2_R( x, H, w, k, o, t, th_2_R), [k.^(-2), (1+k.^(-2))/2], 0)

figure
plot(H_1_L.x,H_1_L.y)
hold on
plot(H_1_R.x,H_1_R.y)
plot(H_2_L.x,H_2_L.y)
plot(H_2_R.x,H_2_R.y)
title('H')
hold off;


 %================_____y_0____=====================
H_1_L=deval(H_1_L,c1)
H_1_R=deval(H_1_R,c1)
H_2_L=deval(H_2_L,c2)
H_2_R=deval(H_2_R,c2)


%================_____y_____=====================
y_0_1=(H_1_L-H_1_R).^(-1/2)
y_0_2=(H_2_L-H_2_R).^(-1/2)

y_1_L=ode45(@(x, y) ODE_y_1_L( x, y, w, k, o, t, fi_1_L), [1/2, 0],  y_0_1)
y_1_R=ode45(@(x, y) ODE_y_1_R( x, y, w, k, o, t, fi_1_R), [1/2, 1], y_0_1)
y_2_L=ode45(@(x, y) ODE_y_2_L( x, y, w, k, o, t, fi_2_L), [(1+k.^(-2))/2, 1], y_0_2)
y_2_R=ode45(@(x, y) ODE_y_2_R( x, y, w, k, o, t, fi_2_R), [(1+k.^(-2))/2, k.^(-2)], y_0_2)

figure
plot(y_1_L.x,y_1_L.y)
hold on
plot(y_1_R.x,y_1_R.y)
plot(y_2_L.x,y_2_L.y)
plot(y_2_R.x,y_2_R.y)
title('Y')
hold off;

 %================_____A_____=====================
h = 0.001


Y_1_L=[]

for i =  0: h: 1/2
   Y_1_L(end+1) =  I_2_1_L( i , w, k, o, t,  th_1_L, y_1_L, fi_1_L)
end
I_2_1=trapz(Y_1_L)


Y_1_R=[]

for i = 1/2: h :1
   Y_1_R(end+1) = I_2_1_R( i , w, k, o, t,  th_1_R, y_1_R, fi_1_R)
end

I_2_2=trapz(Y_1_R)

I_2 = I_2_1 + I_2_2

X_2_L=[]
Y_2_L=[]

for i = 1: (((1+k.^(-2))/2+1)/2)/100:(1+k.^(-2))/2
   X_2_L(end+1) = i
   Y_2_L(end+1) =  J_2_2_L( i , w, k, o, t,  th_2_L, y_2_L, fi_2_L)
 end

J_2_1=trapz(Y_2_L)

X_2_R=[]
Y_2_R=[]

 for i = (1+k.^(-2))/2: (((1+k.^(-2))/2+k.^(-2))/2)/100:(k^(-2))
   X_2_R(end+1) = i
   Y_2_R(end+1) =  J_2_2_R( i , w, k, o, t,  th_2_R, y_2_R, fi_2_R)
 end

J_2_2=trapz(Y_2_R)

J_2 = J_2_1 + J_2_2

C = (J_2-I_2).^(-1/4)

if i2==0
    a1 = (abs(A_2_L( 1 , w, k, o, t,  th_2_L, y_2_L, fi_2_L)/ A_1_R( 1 , w, k, o, t,  th_1_R, y_1_R, fi_1_R)).^(1/2))*C
    a2 = ((a1)^(-1))*C^2
    if A_2_L( 1 , w, k, o, t,  th_2_L, y_2_L, fi_2_L)* A_1_R( 1 , w, k, o, t,  th_1_R, y_1_R, fi_1_R)< 0
        a2 = -a2
    end
else
    a1 = 1
    a2 = 1
end

X_1=[]
A_1=[]

for i =  0: h: 1/2
    X_1(end+1) = i
    A_1(end+1) =  a1 * A_1_L( i , w, k, o, t,  th_1_L, y_1_L, fi_1_L)
end

for i = 1/2: h :1
   X_1(end+1) = i
   A_1(end+1) = a1 * A_1_R( i , w, k, o, t,  th_1_R, y_1_R, fi_1_R)
end

for i =  1: h: (1+k.^(-2))/2
    X_1(end+1) = i
    A_1(end+1) =  a1 * A_2_L( i , w, k, o, t,  th_2_L, y_2_L, fi_2_L)
end

for i = (1+k.^(-2))/2: h : k^(-2)
   X_1(end+1) = i
   A_1(end+1) = a1 * A_2_R( i , w, k, o, t,  th_2_R, y_2_R, fi_2_R)
end

figure
plot(X_1, A_1)
title('A')
