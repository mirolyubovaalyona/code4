function diff_A = diff_A_1_R( x,  w, k, o, t,  H, fi)
    v =(exp(-4/abs(w^2*(2*(x - 1)^2 - 1)^2 + o/k^2 + (2*t*(x - 1)^2)/k^2))*((abs(w^2*(2*(x - 1)^2 - 1)^2 + o/k^2 + (2*t*(x - 1)^2)/k^2)/4)^(1/2) - 1) + 1)^(1/2)
    th= deval(H, x)
    y= deval(fi, x)
    fi= deval(fi, x)
    diff_A = -v*y(2)*cos(th(1))*(sin(th(1)-fi(1))^(-1))
end