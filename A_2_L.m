function A = A_2_L( x,  w, k, o, t, H, fi)
    v =(exp(-4/abs(w^2*((2*(x - 1)^2)/(1/k^2 - 1) + 1)^2 + o/k^2 - (2*t*(x - 1)^2)/(k^2*(1/k^2 - 1))))*((abs(w^2*((2*(x - 1)^2)/(1/k^2 - 1) + 1)^2 + o/k^2 - (2*t*(x - 1)^2)/(k^2*(1/k^2 - 1)))/4)^(1/2) - 1) + 1)^(1/2)
    th= deval(H, x)
    y= deval(fi, x)
    fi= deval(fi, x)
    A = - y(2)*sin(th(1))/(v*sin(th(1)-fi(1)))
end