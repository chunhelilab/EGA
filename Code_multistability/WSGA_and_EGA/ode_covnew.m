%%% ODE for the third moments

function dx = ode_covnew(~, x)
a1 = 1; a2 = 1; b1 = 1; b2 = 1;
k1 = 1; k2 = 1; n = 4; S = 0.5;

x1 = x(9);
x2 = x(10);

dx  = zeros(10,1);
dx(9) = a1 * x1^n / (S^n + x1^n) + b1 * S^n / (S^n + x2^n) - k1 * x1;
dx(10) = a2 * x2^n / (S^n + x2^n) + b2 * S^n / (S^n + x1^n) - k2 * x2;


%% jacobian matrix
A = [ (a1*n*x1^(n - 1))/(S^n + x1^n) - k1 - (a1*n*x1^n*x1^(n - 1))/(S^n + x1^n)^2,  -(S^n*b1*n*x2^(n - 1))/(S^n + x2^n)^2;
    -(S^n*b2*n*x1^(n - 1))/(S^n + x1^n)^2, (a2*n*x2^(n - 1))/(S^n + x2^n) - k2 - (a2*n*x2^n*x2^(n - 1))/(S^n + x2^n)^2];

%% the third moments
y1=x(1);%111
y2=x(2);%211/121/112
y5=x(5);%122/212/221
y8=x(8);%222
dx(1)=3*(A(1,1)*y1+A(1,2)*y2);
dx(2)=2*(A(1,1)*y2+A(1,2)*y5)+A(2,1)*y1+A(2,2)*y2;
dx(3)=2*(A(1,1)*y2+A(1,2)*y5)+A(2,1)*y1+A(2,2)*y2;
dx(4)=2*(A(1,1)*y2+A(1,2)*y5)+A(2,1)*y1+A(2,2)*y2;
dx(5)=A(1,1)*y5+A(1,2)*y8+2*(A(2,1)*y2+A(2,2)*y5);
dx(6)=A(1,1)*y5+A(1,2)*y8+2*(A(2,1)*y2+A(2,2)*y5);
dx(7)=A(1,1)*y5+A(1,2)*y8+2*(A(2,1)*y2+A(2,2)*y5);
dx(8)=3*(A(2,1)*y5+A(2,2)*y8);

end