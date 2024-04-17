function A = A_1_L( x,  w, k, o, t, H, fi, y)
    v =(exp(-4/abs(4*w^2*x^4 + o/k^2 - (t*(2*x^2 - 1))/k^2))*((abs(4*w^2*x^4 + o/k^2 - (t*(2*x^2 - 1))/k^2)/4)^(1/2) - 1) + 1)^(1/2)
    th= deval(H, x)
    y= deval(y, x)
    fi= deval(fi, x)
    A = - y*sin(th(1))/(v*sin(th(1)-fi))
end