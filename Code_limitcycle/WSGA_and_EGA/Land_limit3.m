close all;clear;clc

load momentdata.mat % moment data from LE simulation, see detail in Text S5 of Supplementery Material
%% moments

[t,x_c]=ode45(@ode_cov_limit3,0:0.01:150,[mean;sigma]);

figure(1)
subplot(1, 2, 1)
plot(t, x_c(:, 1),'r', 'linewidth', 1)
hold on
plot(t, x_c(:, 3),'b', 'linewidth', 1)
hold on
plot(t, x_c(:, 5),'y', 'linewidth', 1)

xlabel('t(min)','FontSize',12);
ylabel('x','FontSize',12);
legend('m_1', 'm_2','m_3')
% ylim([0,150])

subplot(1, 2, 2)
sigmat=zeros(6,6,length(t));
for i=1:6
    for j=1:6
        for k=1:length(t)
            sigmat(i,j,k)=x_c(k,6*i+j);
        end
    end
end
sigmat=sigmat(1:2:5,1:2:5,:);

plot(t, reshape(sigmat(1,1,:),length(t),1),'r', 'linewidth', 1)
hold on
plot(t, reshape(sigmat(1,2,:),length(t),1),'b--', 'linewidth', 1)
hold on
plot(t, reshape(sigmat(1,3,:),length(t),1),'m--', 'linewidth', 1)
hold on
plot(t, reshape(sigmat(2,2,:),length(t),1),'b', 'linewidth', 1)
hold on
plot(t, reshape(sigmat(2,3,:),length(t),1),'c--', 'linewidth', 1)
hold on
plot(t, reshape(sigmat(3,3,:),length(t),1),'y', 'linewidth', 1)


xlabel('t(min)','FontSize',12);
ylabel('\sigma','FontSize',12);
legend('\sigma_{11}', '\sigma_{12}','\sigma_{13}',...
 '\sigma_{22}','\sigma_{23}','\sigma_{33}');

%% Landscape的计算
D=0.02;

g = 101; % length of step
pps=zeros(g,g,g);
step=8/(g-1);
for T=100:0.2:116.8 % a period lasts for about 17s
    t=int32(100*T+1);
    p_c = zeros(g,g,g);
    cov=reshape(sigmat(:,:,t),3,3);

    R1=round(x_c(t,1)/step); % approximate range of Gaussian peak to reduce computation
    R2=round(x_c(t,3)/step);
    R3=round(x_c(t,5)/step);
    r=ceil(0.64/step);
    
    for i=max(R1-r,1):min(R1+r,g)
        for j=max(R2-r,1):min(R2+r,g)
            for k=max(R3-r,1):min(R3+r,g)
                p1=(i-1)*step;
                p2=(j-1)*step;
                p3=(k-1)*step;
                p_c(i,j,k) = 1 / sqrt((2*pi)^2 * det( D*cov ))...
                    * exp((-1/2) * ([p1;p2;p3]-[x_c(t,1);x_c(t,3);x_c(t,5)])'...
                    * ( D*cov )^(-1) * ([p1;p2;p3]-[x_c(t,1);x_c(t,3);x_c(t,5)]));
            end
        end
    end
    p_c = p_c / max(p_c,[],'all');
    pps=pps+p_c;
end

pps=pps / sum(sum(sum(pps)));
pps12=sum(pps,3); % projection to 2-D place
pps12=pps12+eps;

figure(2) % plot the landscape
surf(0:step:8, 0:step:8, pps12)
shading interp
xlabel('m_2','FontSize',14);
ylabel('m_1','FontSize',14);
zlabel('P','FontSize',14);
xlim([0, 8])
ylim([0, 8])
view([-29, 59]);