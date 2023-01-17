close all;clear;clc
%% The third moments
[t,y_lh]=ode45(@ode_covnew,0:0.01:8,[1.31;-0.75*ones(3,1);0.40*ones(3,1);-0.21;0.3;2.2]); % the first representative, from the data of LE simulation when t = 0.6s

[~,y_m]=ode45(@ode_covnew,0:0.01:8,[-0.45;-0.69*ones(3,1);2.18*ones(3,1);-3.92;0.55;0.8]); % the second representative, from the data of LE simulation when t = 0.6s

[~,y_hl]=ode45(@ode_covnew,0:0.01:8,[-0.25;0.44*ones(3,1);-0.78*ones(3,1);1.29;2.2;0.3]); % the third representative, from the data of LE simulation when t = 0.6s

figure(1)
font=18;
subplot(1,3,1)
plot(t, y_lh(:, 1),'r--', 'linewidth', 1)
hold on
plot(t, y_lh(:, 2),'b', 'linewidth', 1)
hold on
plot(t, y_lh(:, 5),'g', 'linewidth', 1)
hold on
plot(t, y_lh(:, 8),'black--', 'linewidth', 1)
xlabel('T(s)','FontSize',font);
ylabel('\phi','FontSize',font);
set(gca,'FontSize',font);
legend('\phi_{111}', '\phi_{112}', '\phi_{122}','\phi_{222}','FontSize',font)

subplot(1,3,2)
plot(t, y_m(:, 1),'r--', 'linewidth', 1)
hold on
plot(t, y_m(:, 2),'b', 'linewidth', 1)
hold on
plot(t, y_m(:, 5),'g', 'linewidth', 1)
hold on
plot(t, y_m(:, 8),'black--', 'linewidth', 1)
xlabel('T(s)','FontSize',font);
ylabel('\phi','FontSize',font);
set(gca,'FontSize',font);
legend('\phi_{111}', '\phi_{112}', '\phi_{122}','\phi_{222}','FontSize',font)

subplot(1,3,3)
plot(t, y_hl(:, 1),'r--', 'linewidth', 1)
hold on
plot(t, y_hl(:, 2),'b', 'linewidth', 1)
hold on
plot(t, y_hl(:, 5),'g', 'linewidth', 1)
hold on
plot(t, y_hl(:, 8),'black--', 'linewidth', 1)
xlabel('T(s)','FontSize',font);
ylabel('\phi','FontSize',font);
set(gca,'FontSize',font);
legend('\phi_{111}', '\phi_{112}', '\phi_{122}','\phi_{222}','FontSize',font)

%% Covariance
[~,x_lh] = ode45(@ode_cov,0:0.01:15,[0.3;2.2; 0.2;0;0;0.2]);

[~,x_m] = ode45(@ode_cov,0:0.01:15,[0.55;0.8; 0.4;0;0;0.4]);

[~,x_hl] = ode45(@ode_cov,0:0.01:15,[2.2;0.3; 0.2;0;0;0.2]);

%% Landscape, the same as 'Land.m'
weight = [0.2;0.6;0.2];
D = 0.05;
times = [1.25 1.5 2 2.25 2.5 2.75 3 3.25 3.75 5 6.25 7.5 8.75 10];

DATA = cell(size(times,2),1); % data for different time

for u = 1:size(times,2)
    T = times(u);

    % covariance at stable state
    cov_lh = [x_lh(100*T+1, 3), x_lh(100*T+1, 4); x_lh(100*T+1, 5), x_lh(100*T+1, 6)];
    cov_m = [x_m(100*T+1, 3), x_m(100*T+1, 4); x_m(100*T+1, 5), x_m(100*T+1, 6)];
    cov_hl = [x_hl(100*T+1, 3), x_hl(100*T+1, 4); x_hl(100*T+1, 5), x_hl(100*T+1, 6)];

    g = 201;
    step = 2.5/(g-1);

    x_lh_ = x_lh(100*T+1,1:2);
    y_lh_ = y_lh(100*T+1,[1 2 5 8]);
    x_m_ = x_m(100*T+1,1:2);
    y_m_ = y_m(100*T+1,[1 2 5 8]);
    x_hl_ = x_hl(100*T+1,1:2);
    y_hl_ = y_hl(100*T+1,[1 2 5 8]);

    p_lh = zeros(g^2,1);
    p_m = zeros(g^2,1);
    p_hl = zeros(g^2,1);

    parfor k = 1:g^2
        i = mod(k-1,g)+1;
        j = floor((k-1)/g)+1;
        p1 = (i-1)*step;
        p2 = (j-1)*step;
        fun1 = @(x,y)fUn(x,y,x_lh_(1,1)-p1,x_lh_(1,2)-p2,y_lh_(1,1),y_lh_(1,2),y_lh_(1,3),y_lh_(1,4),cov_lh(1,1),cov_lh(2,1),cov_lh(2,2));
        p_lh(k) = intnew(fun1,-20,20,-20,20,100,100); % Fourier inverse transformation
    end
    p_lh = p_lh / sum(p_lh,"all");

    parfor k = 1:g^2
        i = mod(k-1,g)+1;
        j = floor((k-1)/g)+1;
        p1 = (i-1)*step;
        p2 = (j-1)*step;
        fun2 = @(x,y)fUn(x,y,x_m_(1,1)-p1,x_m_(1,2)-p2,y_m_(1,1),y_m_(1,2),y_m_(1,3),y_m_(1,4),cov_m(1,1),cov_m(2,1),cov_m(2,2));
        p_m(k) = intnew(fun2,-20,20,-20,20,100,100); % Fourier inverse transformation
    end
    p_m = p_m / sum(p_m,"all");

    parfor k = 1:g^2
        i = mod(k-1,g)+1;
        j = floor((k-1)/g)+1;
        p1 = (i-1)*step;
        p2 = (j-1)*step;
        fun3 = @(x,y)fUn(x,y,x_hl_(1,1)-p1,x_hl_(1,2)-p2,y_hl_(1,1),y_hl_(1,2),y_hl_(1,3),y_hl_(1,4),cov_hl(1,1),cov_hl(2,1),cov_hl(2,2));
        p_hl(k) = intnew(fun3,-20,20,-20,20,100,100); % Fourier inverse transformation
    end
    p_hl = p_hl / sum(p_hl,"all");

    p_lh = reshape(p_lh,[g g]);
    p_m = reshape(p_m,[g g]);
    p_hl = reshape(p_hl,[g g]);

    p_lh = p_lh / sum(p_lh,"all");
    p_m = p_m / sum(p_m,"all");
    p_hl = p_hl / sum(p_hl,"all");

    pps = weight(1, 1) * p_lh +  weight(2, 1) * p_m + weight(3,1) * p_hl;
    pps = pps.';

    figure();
    surf(0:step:2.5, 0:step:2.5, pps)
    shading interp
    xlabel('x_1','FontSize',14);
    ylabel('x_2','FontSize',14);
    zlabel('P','FontSize',14);
    xlim([0, 2.5])
    ylim([0, 2.5])
    view([-29, 59]);

    DATA{u} = pps;
end