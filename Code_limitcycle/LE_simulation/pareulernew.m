clear;close all;clc
%%
dt=0.1;
steps=3*1e3;
par_num = 1e4;
x=zeros(6,steps,par_num);

ini=zeros(6,par_num); % initial values
ini(1,:)=ones(1,par_num)+0.1*normrnd(0,1,1,par_num);
ini(2,:)=rand*ones(1,par_num)+0.1*normrnd(0,1,1,par_num);
ini(3,:)=ones(1,par_num)+0.1*normrnd(0,1,1,par_num);
ini(4,:)=rand*ones(1,par_num)+0.1*normrnd(0,1,1,par_num);
ini(5,:)=ones(1,par_num)+0.1*normrnd(0,1,1,par_num);
ini(6,:)=rand*ones(1,par_num)+0.1*normrnd(0,1,1,par_num);

%%
D=0.02;

%%% paraller pool is required
parfor p = 1:par_num
    y = zeros(6,steps);
    y(:,1) = ini(:,p)
    for i = 1:steps-1
        y(:,i+1)=y(:,i)+forcenew(y(:,i))*dt...
            +normrnd(0,sqrt(2*D*dt),6,1);
    end
    x(:,:,p) = y
end

% x=x(1:2:5,:,:);

%%
s = x(:,:,1);
for i = 1:par_num-1
    s = [s x(:,:,i+1)];
end

t1=1e2;
t2=dt*steps;
T1=int64(t1/dt)*par_num+1;
T2=int64(t2/dt)*par_num; % starting time and ending time

[a,xedge,yedge]=histcounts2(s(1,T1:T2),s(3,T1:T2),500);
a=a/sum(sum(a));
a=a+eps;

figure(1)
imshow(-log(a),[])

axis on
colorbar;

%%
y=linspace(min(xedge),max(xedge),size(xedge,2)-1);
z=linspace(min(yedge),max(yedge),size(yedge,2)-1);
[Y,Z]=meshgrid(y,z);

figure(2)
surf(Y,Z,a)
shading interp;
xlim([0,max(xedge)]);
ylim([0,max(yedge)]);
xlabel('lacl','FontSize',14)
ylabel('tetR','Fontsize',14)
zlabel('P','FontSize',14)
