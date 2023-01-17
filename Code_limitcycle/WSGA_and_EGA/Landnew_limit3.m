close all;clear;clc

dim=6;

load momentdata.mat % moment data from LE simulation

%% The third moments

[t,y_c]=ode45(@ode_covnew_limit3,0:0.01:150,[phi;mean]);

%% Covariance

[t,x_c]=ode45(@ode_cov_limit3,0:0.01:150,[mean;sigma]);

sigmat=zeros(6,6,length(t));
for i=1:6
    for j=1:6
        for k=1:length(t)
            sigmat(i,j,k)=x_c(k,6*i+j);
        end
    end
end
sigmat=sigmat(1:2:5,1:2:5,:); % projection

figure(1) % figure of the first three moments
font=18;

subplot(1, 3, 1) % mean
plot(t, y_c(:, dim^3+1),'color',[0.95,0.05,0], 'linewidth', 2.5)
hold on
plot(t, y_c(:, dim^3+3),'color',[0.25,0,0.95], 'linewidth', 2.5)
hold on
plot(t, y_c(:, dim^3+5),'color',[0.95,0.7,0], 'linewidth', 2.5)

xlim([0,110])
% xlabel('T(min)','FontSize',font);
% ylim([0,10])
ylabel('\mu','FontSize',font);
set(gca,'FontSize',font);
legend('\mu_1', '\mu_2','\mu_3','FontSize',font,'Location','northwest','Box','off')

subplot(1,3,2) % covariance
plot(t, reshape(sigmat(1,1,:),length(t),1),'color',[0.95,0.05,0], 'linewidth', 2.5)
hold on
plot(t, reshape(sigmat(2,2,:),length(t),1),'color',[0.25,0,0.95], 'linewidth', 2.5)
hold on
plot(t, reshape(sigmat(3,3,:),length(t),1),'color',[0.95,0.7,0], 'linewidth', 2.5)

xlim([0,110])
% xlabel('T(min)','FontSize',font);
% ylim([0,20])
ylabel('\sigma','FontSize',font);
set(gca,'FontSize',font);
legend('\sigma_{11}', '\sigma_{22}','\sigma_{33}','FontSize',font,'Location','northwest','Box','off');

subplot(1, 3, 3)
plot(t, y_c(:, 1),'color',[0.95,0.05,0], 'linewidth', 2.5)
hold on
plot(t, y_c(:, 87),'color',[0.25,0,0.95], 'linewidth', 2.5)
hold on
plot(t, y_c(:, 173),'color',[0.95,0.7,0], 'linewidth', 2.5)

xlim([0,110])
% xlabel('T(min)','FontSize',font);
% ylim([0,15])
ylabel('\phi','FontSize',font);
set(gca,'FontSize',font);
legend('\phi_{111}', '\phi_{222}', '\phi_{333}','FontSize',font,'Location','northwest','Box','off')

%% Landscape for three mRNA

D=0.02;
g = 101;
step=8/(g-1); % the length of step
pps=zeros(g);
for T=100:0.2:116.8 % a period lasts for about 17s

    t=int32(100.*T+1);
    p_c = zeros(g,g,g); % a Gaussian part of the limit cycle

    X=[x_c(t,1),x_c(t,3),x_c(t,5)]; % the first moments (mean)
    
    sigma=reshape(x_c(t,7:42),6,6).'; % the second moments (variance)
    sigma=sigma(1:2:5,1:2:5);
    
    phi=zeros(6,6,6); % the third moments (skewness)
    for i=1:6
        for j=1:6
            for k=1:6
                phi(i,j,k)=y_c(t,36*i+6*j+k-42);
            end
        end
    end
    phi=phi(1:2:5,1:2:5,1:2:5);
    
    R1=round(x_c(t,1)/step); % approximate range of one peak to reduce computation
    R2=round(x_c(t,3)/step);
    R3=round(x_c(t,5)/step);
    r=ceil(0.64/step);
    
    for i=max(R1-r,1):min(R1+r,g) % local calculation to reduce computation
        for j=max(R2-r,1):min(R2+r,g)
            for k=max(R3-r,1):min(R3+r,g)
                p1=(i-1)*step;
                p2=(j-1)*step;
                p3=(k-1)*step;
                p=[p1,p2,p3];
                % characteristic function
                fun = @(x,y,z)fUn3(x,y,z,x_c(t,1)-p1,x_c(t,3)-p2,x_c(t,5)-p3,y_c(t,1),y_c(t,87),y_c(t,173),y_c(t,3),y_c(t,89),y_c(t,169),...
                    y_c(t,15),y_c(t,101),y_c(t,145),y_c(t,17),sigma(1,1),sigma(2,2),sigma(3,3),sigma(2,1),sigma(3,1),sigma(3,2));
                % Fourier inverse transformation
                p_c(i,j,k)=int3new(fun,-20,20,-20,20,-20,20,100,100,100);
                
            end
        end
    end
    pps=pps+p_c; % summation
end

pps=pps / sum(sum(sum(pps))); % normalization
pps12=sum(pps,3); % project to 2-D place

figure(3) % plot the landscape
surf(0:step:8, 0:step:8, pps12)
shading interp
xlabel('m_2','FontSize',14);
ylabel('m_1','FontSize',14);
zlabel('P','FontSize',14);
xlim([0, 8])
ylim([0, 8])
view([-29, 59]);
