%%% force for tristable genetic circuit
%%% a,b,k: non-negative number
%%% s: 2*2 matrix
%%% x: 2*1 matrix
%%% n: non-negative integer

function F=force(a,b,k,s,x,n)
    F=zeros(2,1);
    F(1)=a*x(1)^n/(s(1,1)^n+x(1)^n)+b*s(2,1)^n/(s(2,1)^n+x(2)^n)-k*x(1);
    F(2)=a*x(2)^n/(s(2,2)^n+x(2)^n)+b*s(1,2)^n/(s(1,2)^n+x(1)^n)-k*x(2);
    
end