%%% function to realize the fast integration method of 27-node interpolation
%%% in dimension 3 (be similar to 'intnew.m')

function  S = int3new(fun,x1,x2,y1,y2,z1,z2,m,n,p)

%%% integral range: [x1,x2][y1,y2][z1,z2]
%%% partition number: m,n,p

dx = (x2-x1)/m/2;
dy = (y2-y1)/n/2;
dz = (z2-z1)/p/2;
x = x1 + dx*(0:(2*m));
y = y1 + dy*(0:(2*n));
z = z1 + dz*(0:(2*p));
[X,Y,Z] = meshgrid(x,y,z);
s=permute(fun(X,Y,Z),[2 1 3]);


s1 = sum(s([1,2*m+1],[1,2*n+1],[1,2*p+1]),'all');
s2 = sum(s([1,2*m+1],[1,2*n+1],3:2:(2*p-1)),'all');
s3 = sum(s([1,2*m+1],3:2:(2*n-1),[1,2*p+1]),'all');
s4 = sum(s(3:2:(2*m-1),[1,2*n+1],[1,2*p+1]),'all');
s5 = sum(s([1,2*m+1],3:2:2*n-1,3:2:2*p-1),'all');
s6 = sum(s(3:2:2*m-1,[1,2*n+1],3:2:2*p-1),'all');
s7 = sum(s(3:2:2*m-1,3:2:2*n-1,[1,2*p+1]),'all');
s8 = sum(s(3:2:2*m-1,3:2:2*n-1,3:2:2*p-1),'all');

s9 = sum(s(2:2:2*m,[1,2*n+1],[1,2*p+1]),'all');
s10 = sum(s([1,2*m+1],2:2:2*n,[1,2*p+1]),'all');
s11 = sum(s([1,2*m+1],[1,2*n+1],2:2:2*p),'all');
s12 = sum(s(2:2:2*m,[1,2*n+1],3:2:(2*p-1)),'all');
s13 = sum(s(3:2:(2*m-1),2:2:2*n,[1,2*p+1]),'all');
s14 = sum(s([1,2*m+1],3:2:(2*n-1),2:2:2*p),'all');
s15 = sum(s(2:2:2*m,3:2:(2*n-1),[1,2*p+1]),'all');
s16 = sum(s([1,2*m+1],2:2:2*n,3:2:(2*p-1)),'all');
s17 = sum(s(3:2:(2*m-1),[1,2*n+1],2:2:2*p),'all');
s18 = sum(s(2:2:2*m,3:2:(2*n-1),3:2:(2*p-1)),'all');
s19 = sum(s(3:2:(2*m-1),2:2:2*n,3:2:(2*p-1)),'all');
s20 = sum(s(3:2:(2*m-1),3:2:(2*n-1),2:2:2*p),'all');

S=dx*dy*dz*(-s1-2*(s2+s3+s4)-4*(s5+s6+s7)-8*s8+4/3*(s9+s10+s11+2*(s12+s13+s14+s15+s16+s17)+4*(s18+s19+s20)));

end