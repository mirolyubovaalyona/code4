function A = A_1_R( x,  w, k, o, t, H, fi, y)
    v =(exp(-4/abs(w^2*(2*(x - 1)^2 - 1)^2 + o/k^2 + (2*t*(x - 1)^2)/k^2))*((abs(w^2*(2*(x - 1)^2 - 1)^2 + o/k^2 + (2*t*(x - 1)^2)/k^2)/4)^(1/2) - 1) + 1)^(1/2)
    th= deval(H, x)
    y= deval(y, x)
    fi= deval(fi, x)
    A = - y*sin(th(1))/(v*sin(th(1)-fi))
end