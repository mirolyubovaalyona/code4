function res=RWEF_13_v2(sqr_ro,sqr_w,l,h,tet_condition,ellips_condition,end_xi,options)
    ro=sqrt(sqr_ro);
    w=sqrt(sqr_w);
    inf_xi=xi_infinity(l,h);
    xi=[ro^2,end_xi,inf_xi]';
    X=set_x(xi,ro);
    n=length(X);

    %% setup 
    temp_R=set_R(w,ro,l,h); %Set R1 R2 ^R2
    P1=temp_R(1);
    P2=temp_R(2);
    P3=temp_R(3);
    %-------------
    arg=[0 ro w l h temp_R];
    %-------------
    teta=zeros(n,1);
    omega=zeros(n,1);
    omega(1)=1;
    if tet_condition==1
        teta(1)=pi/2-v2(arg);
    else
        teta(1)=-v2(arg);
    end
    %-------------
    teta_=zeros(n,1);
    omega_=zeros(n,1);
    omega_(n)=1;

    %%
    %TETA_OMEGA ->
    for i=1:n-1
        y0=[teta(i) 1];
        tspan=[X(i) X(i+1)];

        sol=ode45(@(x,y) set_teta_omega(x,y,ro,w,l,h,P1,P2,P3),tspan,y0,options);
        temp=deval(sol,sol.x);
        omega(i+1)=temp(2,end);
        teta(i+1)=temp(1,end);
    end

    %r <-
    r=R(omega,n,1);
    %L L' ->
    tic
    L=Lambda(X,ro,w,l,h,P1,P2,P3,teta,r,1,n);
    dL=dLambda(X,ro,w,l,h,P1,P2,P3,teta,L,1,n);

    teta_(n)=teta(n)-pi/2;

    %%
    %TETA_OMEGA <-
    for i=n-1:-1:1
        y0=[teta_(i+1) 1];
        tspan=[X(i+1) X(i)];

        sol_=ode45(@(x,y) set_teta_omega(x,y,ro,w,l,h,P1,P2,P3),tspan,y0,options);
        temp=deval(sol_,sol_.x);
        omega_(i)=temp(2,end);
        teta_(i)=temp(1,end);
    end

    %r ->
    r_=R(omega_,1,n);
    %L L' <-
    L_=Lambda(X,ro,w,l,h,P1,P2,P3,teta_,r_,n,1);
    dL_=dLambda(X,ro,w,l,h,P1,P2,P3,teta_,L_,n,1);

    if ellips_condition==1
        Lr1=L(2);
        sqr_Lr1=Lr1^2;
        Lr2=L_(2);
        
        res=(sqr_Lr1-j*Lr2*Lr1)/(sqr_Lr1*Lr2^2);
    else
        dLr1=dL(2);
        sqr_dLr1=dLr1^2;
        dLr2=dL_(2);
        
        res=(sqr_dLr1-j*dLr2*dLr1)/(sqr_dLr1*dLr2^2);
    end


end
%% TETA|OMEGA
function res=set_teta_omega(x,y,ro,w,l,h,P1,P2,P3)   
    arg=[x ro w l h P1 P2 P3];
    %v1|v2----------------------------------
    V1=(v1(arg));
    V2=(v2(arg));
    %dv1|dv2--------------------------------
    DV1=dv1(arg);    
    DV2=dv2(arg);
    %psi------------------------------------
    psi=w*x+V2+y(1);
    %U--------------------------------------
    f_val = f(x,ro);
    df1_val = df1(x,ro);
    df2_val = df2(x,ro);
    K1=f_1(x,ro);
    
    q=(h*ro^2)-(l*(ro^2))*(x^2+ro^2)+(w^2)*(x^2+ro^2)^2;
    
    U=(q/K1)-(df2_val/(2*f_val))+((df1_val^2)/(4*K1));
%TETA-----------------------------------    
    cos_psi = cos(psi);
    sin_psi = sin(psi);
    sin_2psi = sin(2*psi);
    
    F1=(((U/V1)-w-DV2)*(cos_psi)^2)+((V1-w-DV2)*(sin_psi)^2)+(-(DV1/(2*V1))*sin_2psi);     %t11+t12+t13;
%V--------------------------------------        
    %Om11=(DV1/V1)*cos(2*psi);
    %Om12=((U/V1)-V1)*sin(2*psi);
    
    F2=-(((((DV1/V1)*cos(2*psi))+(((U/V1)-V1)*sin(2*psi)))/2)*y(2));
%---------------------------------------
    res=[F1; F2]; 
end
%% R|f*Lambda|f'*Lambda' !! start-1 <-  !!!  start+1 -> !!
function res=R(Omega,start,endl)
    res=zeros(length(Omega),1);
    res(start)=1;
    if start>endl
        k=1; %start=n-1 endl=1
        for i=start-k:-k:endl
            res(i)=Omega(i+k)*res(i+k);
        end
    else
        k=-1; %start=2 endl=n
        for i=start-k:-k:endl
            res(i)=Omega(i+k)*res(i+k);
        end
    end
end

function res=Lambda(x,ro,w,l,h,P1,P2,P3,tet,rr,start,endl)
    res=zeros(length(x),1);
    arg=[0 ro w l h P1 P2 P3];
    if start>endl
        k=-1;
    else
        k=1;
    end
    for i=start:k:endl
        arg(1)=x(i);
        psi=w*x(i)+v2(arg)+tet(i);
        
        temp=(rr(i)/sqrt(v1(arg)))*cos(psi);        %v1 != v^2
        res(i)=abs(temp/sqrt(f(x(i),ro)))*sign(temp/sqrt(f(x(i),ro)));
    end
end

function res=dLambda(x,ro,w,l,h,P1,P2,P3,tet,LL,start,endl)
    res=zeros(length(x),1);
    arg=[0 ro w l h P1 P2 P3];
    if start>endl
        k=-1;
    else
        k=1;
    end
    for i=start:k:endl
       arg(1)=x(i);
       psi=w*x(i)+v2(arg)+tet(i);
       df=df12(x(i),ro);
       
       p1=(LL(i)*sqrt(f(x(i),ro)))/(df);
       p2=-(v1(arg))*tan(psi);
       res(i)=abs(p1*p2)*sign(p1*p2);
    end    
end
%% error check 
function res=check(r,r_,tet,tet_)
    res=zeros(length(r),1);
    for i=1:length(r)
        res(i)=r(i)*r_(i)*sin(tet(i)-tet_(i));
    end
end

%% f|f^2|f'|d(sqrt(f))|f''                          
function res=f(x, ro)
    x2_plus_ro2 = x^2 + ro^2;
    res = sqrt(x2_plus_ro2 * (x2_plus_ro2 - 1));
end

function res=f_1(x, ro)
    x2_plus_ro2 = x^2 + ro^2;
    res = x2_plus_ro2 * (x2_plus_ro2 - 1);
end
function res=df1(x, ro)
    x2_plus_ro2 = x^2 + ro^2;
    res = (2 * x * x2_plus_ro2 + 2 * x * (x2_plus_ro2 - 1)) / (2 * sqrt(x2_plus_ro2 * (x2_plus_ro2 - 1)));
end
function res=df12(x, ro)
    x2_plus_ro2 = x^2 + ro^2;
    res = (2 * x * x2_plus_ro2 + 2 * x * (x2_plus_ro2 - 1)) / (4 * ((x2_plus_ro2 * (x2_plus_ro2 - 1))^(3/4)));
end
function res=df2(x, ro)
    ro2_plus_x2 = ro^2 + x^2;
    ro2_plus_x2_minus_1 = ro2_plus_x2 - 1;
    sqrt_ro2_plus_x2_minus_1 = sqrt(ro2_plus_x2_minus_1);
    
    a=(4*ro^2 + 12*x^2-2)/(2*sqrt_ro2_plus_x2_minus_1);
    b=-(2*x*ro2_plus_x2+2*x*ro2_plus_x2_minus_1)^2/(4*(ro2_plus_x2*ro2_plus_x2_minus_1)^(3/2));
    res=a+b;
end
%v1|v2
function res=v1(arg)
    x=arg(1); ro=arg(2); w=arg(3); P1=arg(6); P3=arg(8);
    x2_plus_ro2 = x^2 + ro^2;
    res = w*(1 + P1/x2_plus_ro2 + P3/x2_plus_ro2^2);
end
function res=v2(arg)
    x=arg(1); ro=arg(2); w=arg(3); P1=arg(6); P3=arg(8);
    x_over_ro = x / ro;
    atan_x_over_ro = atan(x_over_ro);
    x2_plus_ro2 = x^2 + ro^2;
    res = w*((P1/ro)*atan_x_over_ro + (P3/(2*ro^2))*(x/x2_plus_ro2 + (atan_x_over_ro)/ro));
end

%%  d(v1^2)|d(v2)
%dv1=(v1^2)'
function res=dv1(arg)
    x=arg(1); ro=arg(2); w=arg(3); P1=arg(6); P3=arg(8);   
    x2_plus_ro2 = x^2 + ro^2;
    res = -w*((2*P1*x)/(x2_plus_ro2)^2+(4*P3*x)/(x2_plus_ro2)^3);
end
function res=dv2(arg)
    x=arg(1); ro=arg(2); w=arg(3); P1=arg(6); P3=arg(8);
    x_over_ro = x / ro;
    x2_over_ro2_plus_1 = x_over_ro^2 + 1;
    x2_plus_ro2 = x^2 + ro^2;
    res = w*(P1/(ro^2*x2_over_ro2_plus_1) + (P3*(1/(x2_plus_ro2)-(2*x^2)/(x2_plus_ro2)^2 + 1/(ro^2*x2_over_ro2_plus_1)))/(2*ro^2));
end

%% R1|R2|R2_
function res=R_1(w,ro,l)
    res=(w^2-(l*ro^2))/(2*w^2);
end
function res=R_2(w,ro,h,p)
    P1=(h*ro^2)/(2*w^2);
    P2=p*(1-ro^2-p/2);
    P3=-((2*ro^2-1)/(4*w^2));
    res=P1+P2+P3;
end
function res=R_2_(ro,p1,p2)
    res=p2+p1*ro^2;
end   
function res=set_R(w,ro,l,h)
    P1=R_1(w,ro,l);
    P2=R_2(w,ro,h,P1);
    P3=R_2_(ro,P1,P2);
    
    res=[P1 P2 P3];
end
function res=set_x(xi,ro)
    res = sqrt(xi - ro^2);
end
function res = xi_infinity(l,h)
    pop = (h+l)/2;
    if pop>1000
        res=10;
    elseif pop>100
        res=100;
    else
        res=1000;
    end
end
%%                              end of real functions