function A = A_2_R( x,  w, k, o, t, th, y, fi)
    q=o*(k^(-2))-t*(k^(-2))*((k^(-2))-1-2*((x-(k^(-2)))^2)*(((k^(-2))-1)^(-1)))+(w^2)*(((k^(-2))-2*(((x-(k^(-2)))^2))/((k^(-2))-1))^2)
    v2=(((abs(q)/4)^(1/2))-1)*exp(-4/abs(q))+1
    v=v2^(1/2)
    th= deval(th, x)
    y= deval(y, x)
    fi= deval(fi, x)
    A = - y*sin(th)/(v*sin(th-fi))

end