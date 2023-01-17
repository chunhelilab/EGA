%%% the characteristic function with first three moments

function f = fUn(x,y,u,v,a1,a2,a3,a4,s1,s2,s3)
D=0.05;
f = cos(u.*x+v.*y-D^(1.5)*(a1.*x.^3+3*a2.*x.^2.*y+3*a3.*x.*y.^2+a4.*y.^3)/6).*exp((-D.*(s1.*x.^2+2*s2.*x.*y+s3.*(y.^2)))/2)/(2*pi)^2;
end