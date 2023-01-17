%%% ODE model of GA

function dx = ode_cov_limit3(~, x)
%% Jacobi
x1 = x(1);
x2 = x(2);
x3 = x(3);
x4 = x(4);
x5 = x(5);
x6 = x(6);
n=2;
a=10;
a0=1e-3*a;
b=0.5;

dx(1) = -x1+a/(1+x6^n)+a0;
dx(2) = -b*(x2-x1);
dx(3) = -x3+a/(1+x2^n)+a0;
dx(4) = -b*(x4-x3);
dx(5) = -x5+a/(1+x4^n)+a0;
dx(6) = -b*(x6-x5);

A = [-1,0,0,0,0,-a*2*x6/(1+x6^n)^2;
    b,-b,0,0,0,0;
    0,-a*2*x2/(1+x2^n)^2,-1,0,0,0;
    0,0,b,-b,0,0;
    0,0,0,-a*2*x4/(1+x4^n)^2,-1,0;
    0,0,0,0,b,-b];

%% Covariance
x_sig = reshape(x(7:end), 6, 6);
weak=x_sig * A' + A * x_sig +2*eye(6);
% weak relation from MFA approach
for i=1:2:5
    for j=2:2:6
        weak(i,j)=0;
    end
end
for i=2:2:6
    for j=1:2:5
        weak(i,j)=0;
    end
end

dx = [dx(:); reshape(weak, 36, 1)];
end