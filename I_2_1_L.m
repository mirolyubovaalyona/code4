function [I_2] = I_2_1_L( x,  w, k, o, t, th, y, fi)
    f=(2*(x^2)-1)*(2*(x^2)-(k^(-2)))/8
    q=o*(k^(-2))+(1-2*(x^2))*t*(k^(-2))+4*(w^2)*(x^4)
    v2 =exp(-4/abs(4*w^2*x^4 + o/k^2 - (t*(2*x^2 - 1))/k^2))*((abs(4*w^2*x^4 + o/k^2 - (t*(2*x^2 - 1))/k^2)/4)^(1/2) - 1) + 1
    v =(exp(-4/abs(4*w^2*x^4 + o/k^2 - (t*(2*x^2 - 1))/k^2))*((abs(4*w^2*x^4 + o/k^2 - (t*(2*x^2 - 1))/k^2)/4)^(1/2) - 1) + 1)^(1/2)
    th= deval(th, x)
    y= deval(y, x)
    fi= deval(fi, x)
    A = - y*sin(th)/(v*sin(th-fi))
    I_2 = x*(A.^2)/(f.^(1/2))
end