%%% ODE model of EGA

function dx = ode_covnew_limit3(t, x)
dim=6;

%% Jacobi
y1 = x(dim^3+1);
y2 = x(dim^3+2);
y3 = x(dim^3+3);
y4 = x(dim^3+4);
y5 = x(dim^3+5);
y6 = x(dim^3+6);
n=2;
a=10;
a0=1e-3*a;
b=0.5;

dx(dim^3+1,1) = -y1+a/(1+y6^n)+a0;
dx(dim^3+2,1) = -b*(y2-y1);
dx(dim^3+3,1) = -y3+a/(1+y2^n)+a0;
dx(dim^3+4,1) = -b*(y4-y3);
dx(dim^3+5,1) = -y5+a/(1+y4^n)+a0;
dx(dim^3+6,1) = -b*(y6-y5);

A = [-1,0,0,0,0,-a*2*y6/(1+y6^n)^2;
    b,-b,0,0,0,0;
    0,-a*2*y2/(1+y2^n)^2,-1,0,0,0;
    0,0,b,-b,0,0;
    0,0,0,-a*2*y4/(1+y4^n)^2,-1,0;
    0,0,0,0,b,-b];

%% the third moments
x_phi=reshape(x(1:dim^3),dim,dim,dim); % the third moments

dphi=permute(pagemtimes(A,permute(x_phi,[3 2 1])),[1 3 2])...
    +permute(pagemtimes(A,permute(x_phi,[3 2 1])),[2 1 3])...
    +permute(pagemtimes(A,permute(x_phi,[3 2 1])),[3 2 1]);


% weak relation from MFA approach
for i=1:dim
    for j=1:dim
        for k=1:dim
            if mod(i-j,2)==0 && mod(j-k,2)==0
                continue
            else
                dphi(i,j,k)=0;
            end
        end
    end
end

dx(1:dim^3,1)=reshape(dphi,dim^3,1);
end