i1=0
i2=1
i3=0
w2=500
k2=0.5
w=sqrt(w2)
k=sqrt(k2)
m=25
n=45

h = 4.6506e+03
l = 4.6937e+03


t = l
o = h-l

options=odeset('AbsTol',1e-06,'RelTol',1e-03,'Refine','8')

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

%================_____H_____=====================
H_1_L=ode45(@(x, H) ODE_H_0_1_L( x, H, w, k, o, t), [0, 1/2], [th0_1, 0, 0], options)
H_1_R=ode45(@(x, H) ODE_H_0_1_R( x, H, w, k, o, t), [1, 1/2], [th0_2, 0, 0], options)
H_2_L=ode45(@(x, H) ODE_H_0_2_L( x, H, w, k, o, t), [1, (1+k.^(-2))/2], [th0_3 ,0, 0], options)
H_2_R=ode45(@(x, H) ODE_H_0_2_R( x, H, w, k, o, t), [k.^(-2), (1+k.^(-2))/2], [th0_4, 0, 0], options)


figure
plot(H_1_L.x,H_1_L.y)
hold on
plot(H_1_R.x,H_1_R.y)
plot(H_2_L.x,H_2_L.y)
plot(H_2_R.x,H_2_R.y)
title('H_0')
hold off;


 %================_____fi_0______=====================
th = deval(H_1_L,c1)
fi0_1=th(1)+pi/2
th = deval(H_1_R,c1)
fi0_2=th(1)+pi/2
th = deval(H_2_L,c2)
fi0_3=th(1)+pi/2
th = deval(H_2_R,c2)
fi0_4=th(1)+pi/2


 %================_____y_0____=====================
H_0_1_L = deval(H_1_L,c1)
H_0_1_R = deval(H_1_R,c1)
H_0_2_L = deval(H_2_L,c2)
H_0_2_R = deval(H_2_R,c2)

y_0_1=(H_0_1_L(2)-H_0_1_R(2))^(-1/2)
y_0_2=(H_0_2_L(2)-H_0_2_R(2))^(-1/2)

 %================_____fi_____=====================
fi_1_L=ode45(@(x, fi) ODE_fi_1_L( x, fi, w, k, o, t), [1/2, 0], [fi0_1, y_0_1], options)
fi_1_R=ode45(@(x, fi) ODE_fi_1_R( x, fi, w, k, o, t), [1/2, 1], [fi0_1, y_0_1], options)
fi_2_L=ode45(@(x, fi) ODE_fi_2_L( x, fi, w, k, o, t), [(1+k.^(-2))/2, 1], [fi0_2, y_0_2], options)
fi_2_R=ode45(@(x, fi) ODE_fi_2_R( x, fi, w, k, o, t), [(1+k.^(-2))/2, k.^(-2)], [fi0_2, y_0_2], options)

figure
plot(fi_1_L.x,fi_1_L.y)
hold on
plot(fi_1_R.x,fi_1_R.y)
plot(fi_2_L.x,fi_2_L.y)
plot(fi_2_R.x,fi_2_R.y)
title('fi')
hold off;


 %================_____A_____=====================
H_1_1_L = deval(H_1_L,c1)
H_1_1_R = deval(H_1_R,c1)
H_1_2_L = deval(H_2_L,c2)
H_1_2_R = deval(H_2_R,c2)

h = 0.001
I2 = (H_1_1_L(3) - H_1_1_R(3)) / (H_0_1_L(2) - H_0_1_R(2))
J2 = (H_1_2_L(3) - H_1_2_R(3)) / (H_0_2_L(2) - H_0_2_R(2))


%================_____a1, a2_____=====================
C = (J2-I2)^(-1/4)
if i2==0
    a1 = ((abs(A_2_L( 1 , w, k, o, t,  H_2_L, fi_2_L)/ A_1_R( 1 , w, k, o, t,  H_1_R, fi_1_R)))^(1/2))*C
    a2 = (a1^(-1))*(C^2)
    if A_2_L( 1 , w, k, o, t,  H_2_L,  fi_2_L)* A_1_R( 1 , w, k, o, t,  H_1_R, fi_1_R) < 0
        a2 = -a2
    end
else
    a1 = ((abs(diff_A_2_L( 1 , w, k, o, t, H_2_L, fi_2_L) * ((diff_A_1_R( 1, w, k, o, t,  H_1_R, fi_1_R))^(-1))))^(1/2))*C 
    a2 = (a1^(-1))*(C^2)
    if A_2_L( 1 , w, k, o, t, H_2_L, fi_2_L)* A_1_R( 1 , w, k, o, t,  H_1_R, fi_1_R)< 0
        a2 = -a2
    end
end

X_1=[]
A_1=[]

for i =  0: h: 1/2
    X_1(end+1) = i
    A_1(end+1) =  a1 * A_1_L( i , w, k, o, t,  H_1_L, fi_1_L)
end

for i = 1/2: h :1
   X_1(end+1) = i
   A_1(end+1) = a1 * A_1_R( i , w, k, o, t,  H_1_R,  fi_1_R)
end

for i =  1: h: (1+k.^(-2))/2
    X_1(end+1) = i
    A_1(end+1) =  a1 * A_2_L( i , w, k, o, t,  H_2_L, fi_2_L)
end

for i = (1+k.^(-2))/2: h : k^(-2)-h
   X_1(end+1) = i
   A_1(end+1) = a1 * A_2_R( i , w, k, o, t,  H_2_R,  fi_2_R)
end

figure
plot(X_1, A_1)
grid on
title('A')


