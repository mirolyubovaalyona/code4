function diff_A = diff_A_1_L( x,  w, k, o, t, th,  H, fi)
    v =(exp(-4/abs(4*w^2*x^4 + o/k^2 - (t*(2*x^2 - 1))/k^2))*((abs(4*w^2*x^4 + o/k^2 - (t*(2*x^2 - 1))/k^2)/4)^(1/2) - 1) + 1)^(1/2)
    th= deval(H, x)
    y= deval(fi, x)
    fi= deval(fi, x)
    diff_A = -v*y(2)*cos(th(1))*(sin(th(1)-fi(1))^(-1))
end