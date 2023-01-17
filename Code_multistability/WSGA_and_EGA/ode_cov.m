%%% ODE for the first two moments

function dx = ode_cov(~, x)
a1 = 1; a2 = 1; b1 = 1; b2 = 1;
k1 = 1; k2 = 1; n = 4; S = 0.5;

x1 = x(1);
x2 = x(2);

dx = zeros(2,1);
dx(1) = a1 * x1^n / (S^n + x1^n) + b1 * S^n / (S^n + x2^n) - k1 * x1;
dx(2) = a2 * x2^n / (S^n + x2^n) + b2 * S^n / (S^n + x1^n) - k2 * x2;


%% jacobian matrix
A = [ (a1*n*x1^(n - 1))/(S^n + x1^n) - k1 - (a1*n*x1^n*x1^(n - 1))/(S^n + x1^n)^2,  -(S^n*b1*n*x2^(n - 1))/(S^n + x2^n)^2;
    -(S^n*b2*n*x1^(n - 1))/(S^n + x1^n)^2, (a2*n*x2^(n - 1))/(S^n + x2^n) - k2 - (a2*n*x2^n*x2^(n - 1))/(S^n + x2^n)^2];
 
%% variance
x_sig = reshape(x(3:end), 2, 2);
dx = [dx; reshape(x_sig * A.' + A * x_sig + 2 * eye(2), 4, 1)]; % from Omega-expansion
end