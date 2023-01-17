%%% function to realize the fast integration method of 8-node interpolation
%%% in dimension 2

function  S = intnew(fun,x1,x2,y1,y2,n,m)

%%% integral range: [x1,x2][y1,y2]
%%% partition number: n,m

dx = (x2-x1)/n/2;
dy = (y2-y1)/m/2;
x = x1 + dx*(0:(2*n));
y = y1 + dy*(0:(2*m));
[X,Y] = meshgrid(x,y);
s = fun(X,Y)';% (2n+1) x (2m+1)
s1 = sum(s(3:2:(2*n-1),3:2:(2*m-1)),'all'); % The "inside" vertex with coefficient -4/3 (counted four times)
s2 = sum(s(3:2:(2*n-1),2:2:(2*m)),'all'); % The "inside" midpoint with coefficient 8/3 (counted twice)
s3 = sum(s(2:2:(2*n),3:2:(2*m-1)),'all'); % The "inside" midpoint with coefficient 8/3 (counted twice)
s4 = sum(s([1,2*n+1],2:2:(2*m)),'all'); % The "boundary" midpoint with coefficient 4/3 (counted once)
s5 = sum(s(2:2:(2*n),[1,2*m+1]),'all'); % The "boundary" midpoint with coefficient 4/3 (counted once)
s6 = sum(s([1,2*n+1],[1,2*m+1]),'all'); % The vertex with coefficient -1/3 (counted once)
s7 = sum(s([1,2*n+1],3:2:(2*m-1)),'all'); % The "boundary" vertex with coefficient -2/3 (counted twice)
s8 = sum(s(3:2:(2*n-1),[1,2*m+1]),'all'); % The "boundary" vertex with coefficient -2/3 (counted twice)
S = dx*dy/3*(8*(s2+s3)+4*(s4+s5-s1)-s6-2*(s7+s8)); % weighted summation
end