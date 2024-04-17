function diff_A = diff_A_2_R( x,  w, k, o, t, H, fi, y)
    v =(exp(-4/abs(w^2*((2*(x - 1/k^2)^2)/(1/k^2 - 1) - 1/k^2)^2 + o/k^2 + (t*((2*(x - 1/k^2)^2)/(1/k^2 - 1) - 1/k^2 + 1))/k^2))*((abs(w^2*((2*(x - 1/k^2)^2)/(1/k^2 - 1) - 1/k^2)^2 + o/k^2 + (t*((2*(x - 1/k^2)^2)/(1/k^2 - 1) - 1/k^2 + 1))/k^2)/4)^(1/2) - 1) + 1)^(1/2)
    th= deval(H, x)
    y= deval(y, x)
    fi= deval(fi, x)
    diff_A = -v*y*cos(th(1))*(sin(th(1)-fi)^(-1))
end