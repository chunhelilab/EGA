close all;clear;clc
%% different initial conditions
[t,x_lh]=ode45(@ode_cov,0:0.01:15,[0.3;2.2; 2;-0.5*ones(2,1);1]);
figure(1)
subplot(1, 2, 1)
plot(t, x_lh(:, 1),'r', 'linewidth', 1)
hold on
plot(t, x_lh(:, 2), 'g--', 'linewidth', 1)
xlabel('t(s)','FontSize',12);
ylabel('x','FontSize',12);
legend('x_1', 'x_2')
subplot(1, 2, 2)
plot(t, x_lh(:, 3),'r', 'linewidth', 1)
hold on
plot(t, x_lh(:, 4),'b', 'linewidth', 1)
hold on
plot(t, x_lh(:, 6),'g--', 'linewidth', 1)
xlabel('t(s)','FontSize',12);
ylabel('\sigma','FontSize',12);
legend('\sigma_{11}', '\sigma_{12}','\sigma_{22}')

[t,x_m]=ode45(@ode_cov,0:0.01:15,[0.55;0.8; 1;0*ones(2,1);1]);
figure(2)
subplot(1, 2, 1)
plot(t, x_m(:, 1),'r', 'linewidth', 1)
hold on
plot(t, x_m(:, 2), 'g--', 'linewidth', 1)
xlabel('t(s)','FontSize',12);
ylabel('x','FontSize',12);
legend('x_1', 'x_2')
subplot(1, 2, 2)
plot(t, x_m(:, 3),'r', 'linewidth', 1)
hold on
plot(t, x_m(:, 4),'b', 'linewidth', 1)
hold on
plot(t, x_m(:, 6),'g--', 'linewidth', 1)
xlabel('t(s)','FontSize',12);
ylabel('\sigma','FontSize',12);
legend('\sigma_{11}', '\sigma_{12}','\sigma_{22}')

[t,x_hl]=ode45(@ode_cov,0:0.01:15,[2.2; 0.3; 1;-0.5*ones(2,1);2]);
figure(3)
subplot(1, 2, 1)
plot(t, x_hl(:, 1),'r', 'linewidth', 1)
hold on
plot(t, x_hl(:, 2), 'g--', 'linewidth', 1)
xlabel('t(s)','FontSize',12);
ylabel('x','FontSize',12);
legend('x_1', 'x_2')
subplot(1, 2, 2)
plot(t, x_hl(:, 3),'r', 'linewidth', 1)
hold on
plot(t, x_hl(:, 4),'b', 'linewidth', 1)
hold on
plot(t, x_hl(:, 6),'g--', 'linewidth', 1)
xlabel('t(s)','FontSize',12);
ylabel('\sigma','FontSize',12);
legend('\sigma_{11}', '\sigma_{12}', '\sigma_{22}')

%% Landscape

weight = [0.2;0.6;0.2];
D=0.05;
T=10; %observed time

cov_lh = [x_lh(100*T+1, 3), x_lh(100*T+1, 4); x_lh(100*T+1, 5), x_lh(100*T+1, 6)]; % covariance
cov_m = [x_m(100*T+1, 3), x_m(100*T+1, 4); x_m(100*T+1, 5), x_m(100*T+1, 6)];
cov_hl = [x_hl(100*T+1, 3), x_hl(100*T+1, 4); x_hl(100*T+1, 5), x_hl(100*T+1, 6)];

g = 201; % length of step
step=2.5/(g-1);
        
p_lh = zeros(g);
p_m =zeros(g);
p_hl = zeros(g);
for i=1:g
    for j=1:g
        p1=(i-1)*step;
        p2=(j-1)*step;
        p_lh(i, j) = 1 / sqrt((2*pi)^2 * D*det(cov_lh)) * exp((-1/2) * ([p1;p2]-[x_lh(100*T+1,1);x_lh(100*T+1,2)])' * (D*cov_lh)^(-1) * ([p1;p2]-[x_lh(100*T+1,1);x_lh(100*T+1,2)]));
        p_m(i, j) = 1 / sqrt((2*pi)^2 * D*det(cov_m)) * exp((-1/2) * ([p1;p2]-[x_m(100*T+1,1);x_m(100*T+1,2)])' * (D*cov_m)^(-1) * ([p1;p2]-[x_m(100*T+1,1);x_m(100*T+1,2)]));
        p_hl(i, j) = 1 / sqrt((2*pi)^2 * D*det(cov_hl)) * exp((-1/2) * ([p1;p2]-[x_hl(100*T+1,1);x_hl(100*T+1,2)])' * (D*cov_hl)^(-1) * ([p1;p2]-[x_hl(100*T+1,1);x_hl(100*T+1,2)]));
    end
end
p_lh = p_lh / sum(sum(p_lh));
p_m = p_m / sum(sum(p_m));
p_hl = p_hl / sum(sum(p_hl));

pps = weight(1, 1) * p_lh + weight(2, 1) * p_m + weight(3,1) * p_hl; % weighted summation
pps=pps';

figure(4);
surf(0:step:2.5, 0:step:2.5, pps)
shading interp
xlabel('m_1','FontSize',14);
ylabel('m_2','FontSize',14);
zlabel('P','FontSize',14);
xlim([0, 2.5])
ylim([0, 2.5])
view([-29, 59]);