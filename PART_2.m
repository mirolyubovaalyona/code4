%8

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
 

A=[]
x=[]
 
 for i = 0: (1/2)/100:1/2
     x(end+1) = i 
     y=deval(y_1_L,i)
     A(end+1) = -y*sin(deval(th_1_L,i))/(v_1(i, w, k, o, t)*sin(deval(th_1_L,i)-deval(fi_1_L,i)))
 end
 
  for i = 1/2: (1/2)/100:1
     x(end+1) = i
     y=deval(y_1_R,i)
     A(end+1) = - y*sin(deval(th_1_R,i))/(v_2(i, w, k, o, t)*sin(deval(th_1_R,i)-deval(fi_1_R,i)))
 end
 
  for i = 1: (((1+k.^(-2))/2+1)/2)/100:(1+k.^(-2))/2
     x(end+1) = i 
     y=deval(y_2_L,i)
     A(end+1) = -y*sin(deval(th_2_L,i))/(v_3(i, w, k, o, t)*sin(deval(th_2_L,i)-deval(fi_2_L,i)))
 end
 
  for i = (1+k.^(-2))/2: (((1+k.^(-2))/2+k.^(-2))/2)/100:(k^(-2))
     x(end+1) = i
     y=deval(y_2_R,i)
     A(end+1) = - y*sin(deval(th_2_R,i))/(v_4(i, w, k, o, t)*sin(deval(th_2_R,i)-deval(fi_2_R,i)))
  end
 
figure
plot(x, A)
title('A')