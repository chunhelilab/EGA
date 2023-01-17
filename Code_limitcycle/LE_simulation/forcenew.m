%%% force of the synthetic oscillatory network

function F=forcenew(x)

    a=10;
    a0=1e-3*a;
    b=0.5;
    n=2;

    F=zeros(6,1);
    F(1) = -x(1)+a/(1+x(6)^n)+a0;
    F(2) = -b*(x(2)-x(1));
    F(3) = -x(3)+a/(1+x(2)^n)+a0;
    F(4) = -b*(x(4)-x(3));
    F(5) = -x(5)+a/(1+x(4)^n)+a0;
    F(6) = -b*(x(6)-x(5));

end