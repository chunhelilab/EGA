clear;close all;clc

dt = 0.01;
steps = 1e4;
par_num = 3*1e3 ; % Parallel quantity
x = zeros(2,steps,par_num);

%% initial values

ini1 = repmat([2.2;0.3],1,par_num/3)+0.2*normrnd(0,1,2,par_num/3);
ini2 = repmat([0.55;0.8],1,par_num/3)+0.4*normrnd(0,1,2,par_num/3);
ini3 = repmat([0.3;2.2],1,par_num/3)+0.2*normrnd(0,1,2,par_num/3);

ini=[ini1 ini2 ini3]; 

%% parameters in force
a=1;
b=1;
k=1;
s=[0.5,0.5;0.5,0.5];
n=4;
D=0.01;

%%% parallel pool is required
%% Explict Euler
parfor p = 1:par_num
    y = zeros(2,steps)
    y(:,1) = ini(:,p)
    for i = 1:steps-1
        y(:,i+1)=y(:,i)+force(a,b,k,s,y(:,i),n)*dt...
            +normrnd(0,sqrt(2*D*dt),2,1);
    end
    x(:,:,p) = y
end

s = x(:,:,[1]);
for i = 1:par_num-1
    s = [s x(:,:,[i+1])];
end

[landscape,xedge,yedge]=histcounts2(s(1,:),s(2,:),201);

landscape=landscape/sum(sum(landscape));
landscape=landscape.';

%% Figure
figure(1)
imshow(landscape,[])
axis on
colorbar;

y=linspace(min(xedge),max(xedge),size(xedge,2)-1);
z=linspace(min(yedge),max(yedge),size(yedge,2)-1);
[Y,Z]=meshgrid(y,z);
figure(2)
surf(Y,Z,landscape)
shading interp;
xlim([0 2.5]);
ylim([0 2.5]);
