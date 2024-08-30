%function P_k_store = DualControl(HL)
%该代码是在原始DCEE和RHC结合的基础上输出各种目的图
% 其他代码
close all
clear all
clc
% 
% function P_k_store=DualControl(HL);

addpath('Functions')
UAVVel=2; %UAV velocity
sampleTime=10;
bLim=400;
% Simulated source parameters
% true source
s.Q = 20; % Release rate per t
% source coodinates
s.x = 5;
s.y = 8;
s.z= 1;
s.u = 10;
s.phi = 1;
s.ci = 0.1;%0.14;  % Also s.D for Pasquil model
s.cii = 0.6;%0.53; % Also s.t or tau for Pasquil Model
% s.ci = 0.14;
% s.cii = 0.53;
% s.duration = 0;
m.thresh = 5e-4;

% Create rectangular domain area
xmin = 0;
xmax = 10;
ymin = 0;
ymax = 10;
zmin = 0;
zmax = 5;
domain = [xmin xmax ymin ymax]; % Size of search area

% example data
stepsize = 0.1; %horisontal (x and y) spacing

x_coord = xmin : stepsize : xmax;
y_coord = ymin : stepsize : ymax;
z_coord = zmin : stepsize : zmax;

% Create 3D grid
[X,Y,Z] = meshgrid(x_coord,y_coord,z_coord);

% hold grid information and is passed as an argument to the RadioactiveDispersionModel 
% function to calculate the concentration of the dispersion at different grid points
ex.x_matrix = X;
ex.y_matrix = Y;
ex.z_matrix = Z;

% Initialisation and parameters
StartingPosition = [2 2,4]; % Starting position [x,y,z]
moveDist = 0.5; % How far to move


P_k = StartingPosition; % Current position
P_k_store = [];
P_k_store = [P_k_store; P_k];

pos.x_matrix = P_k(1);
pos.y_matrix = P_k(2);
pos.z_matrix = P_k(3);

% Plot example dispersion from true source
% conc = simpleGaussianPlume(s,m,ex);
conc = RadioactiveDispersionModel(s,ex);
figure(1)
hold off
height = 1;
concSurf=conc(:,:,height)/s.Q;
concSurf(concSurf<=m.thresh)=NaN;
pcolor(ex.x_matrix(:,:,height),ex.y_matrix(:,:,height),concSurf);
axis([xmin xmax ymin ymax]);
% slice(conc,28,17,1);
shading interp
% axis equal
% view(0,90)
xlab = xlabel('x');
ylab = ylabel('y');
colorbar
hold on
plot(s.x,s.y,'k.','MarkerSize',30)
plot(P_k(1),P_k(2),'rx')
plot(P_k_store(:,1),P_k_store(:,2),'r')
S = [];
colorbar

% initialise PF
N = 20000; %20000
PF_Memory=10;
resample = 0;
% Uniform prior for location
theta.x = xmin + (xmax-xmin) * rand(N,1);
theta.y = ymin + (ymax-ymin) * rand(N,1);
theta.z = ones(N,1)*s.z;
a = ones(N,1)*2;
b = ones(N,1)*5;
theta.Q = gamrnd(a,b);%200*rand(N,1);%
figure;histogram(theta.Q)
theta.u =s.u*ones(N,1);%2+6*rand(N,1);%0.75+0.5*rand(N,1);0 + randn(N,1)*0.5;%
theta.phi = s.phi*ones(N,1);%(10 + 30*rand(N,1)).*pi/180;
theta.ci = s.ci*ones(N,1);%0.12+0.1*rand(N,1);
theta.cii =s.cii*ones(N,1);%0.5+ 0.1*rand(N,1);
%Wp refers to particle weights
Wp = ones(N,1);
Wpnorm = Wp./sum(Wp);
Wp = Wpnorm;
  
figure(4)
hold on
% preprocess(s,theta,Wpnorm);

timestamp(1)=0;
D=[];
sampleHistory=P_k(1:2);
% step_no=0;%calculat the No. of steps;
% distance=0;%calculat the distance from start position;
% 初始化绘图所需的变量
total_steps = [];
total_distances = [];
distance_errors = [];
stepSizes = [];
covariances = [];
timestamp = zeros(1, 101); % 假设最多100次迭代

for i = 1:100
    % simulated data
    Dsim=PredictMeasurement(s,pos,m.thresh);
    D(i)=Dsim;

    thetaPrev=theta;
    [theta,Wpnorm]=UpdatePFPlume(D,theta,Wpnorm,pos,P_k_store,m,N,PF_Memory,domain);
    
    figure(1)
    hold off
    concSurf=conc(:,:,height)/s.Q;
    concSurf(concSurf<=0.045)=NaN;
    pcolor(ex.x_matrix(:,:,1),ex.y_matrix(:,:,1),log10(concSurf))
    c1=min(min(concSurf));
    c2=max(max(concSurf));
    clim([-1.3 log10(c2)]);
    c=colorbar;
    c.Label.String='Concentration log(kg/m^3)';
    c.Limits=[-1.3,0];
    shading interp
    grid on
    hold on
    S(i) = 5+ceil(D(i)*5e4);
    scatter3(theta.x,theta.y,theta.z,3,'g','filled')
    plot3(pos.x_matrix,pos.y_matrix,pos.z_matrix,'ro','MarkerFaceColor','r','MarkerSize',5)
    plot3(P_k_store(:,1),P_k_store(:,2),P_k_store(:,3),'r--','LineWidth',2)
    plot3(s.x,s.y,s.z,'k.','markersize',20)
    xlab = xlabel('x (m)');
    ylab = ylabel('y (m)');
    set(xlab,'FontSize',16);
    set(ylab,'FontSize',16);
    set(gca,'fontsize',16)

    view(0,90)
    
    axis([xmin xmax ymin ymax])
    drawnow
    
%     if bLim<=810
%         if bLim<=630
%             if bLim<=400
%                 if bLim<=0
%                     true
%                 end
%             end
%         end
%     end
    tic
    
    M= 40;%16;%5
    MM = 1;
    
    % horizon length
    HL = 3;

    % calculate the next position using RHC
    [~, pos] = RHC(HL,1,moveDist,pos,pos,theta,Wpnorm,D,m,s,M,MM,i,P_k_store,N,PF_Memory,domain,xmin,xmax,ymin,ymax,zmin,zmax);
    
    P_k = [pos.x_matrix pos.y_matrix pos.z_matrix];
    sampleHistory=[sampleHistory;P_k(1:2)];
    
    % 记录数据
distance_error = norm(P_k(1:2) - [s.x, s.y]);
Covar = cov(theta.x, theta.y);
Spread = sqrt(Covar(1,1) + Covar(2,2));

total_steps(i) = i;
total_distances(i) = sum(sqrt(sum(diff(P_k_store).^2, 2)));
distance_errors(i) = distance_error;
stepSizes(i) = moveDist;
covariances(i) = Spread;

    move_time=floor(norm(P_k-P_k_store(i,:))/UAVVel)+sampleTime;
    bLim=bLim-move_time
    if bLim<=0
        break
    end
    timestamp(i+1)=timestamp(i)+move_time;
    P_k_store = [P_k_store; P_k];
    
    
    
    Covar = cov(theta.x,theta.y);
    Spread = sqrt(Covar(1,1)+Covar(2,2))

    %  if norm(P_k(1:2)-[s.x s.y])==0
    %     break
    % end
    
    if Spread<0 % 3.5? 4? 5 default
        break
    end
    
end



indx = resampleStratified(Wpnorm);

theta.x = theta.x(indx);
theta.y = theta.y(indx);
theta.z = theta.z(indx);
theta.Q = theta.Q(indx);
theta.u = theta.u(indx);
theta.phi = theta.phi(indx);
theta.ci = theta.ci(indx);
theta.cii = theta.cii(indx);
figure(4)
hold on
preprocess(s,theta,Wpnorm);
% 绘图部分
% 绘制距离误差随迭代次数变化的图
figure;
plot(1:length(total_steps), distance_errors, '-o', 'Color', [0.6350 0.0780 0.1840], 'LineWidth', 2);
xlabel('Iteration', 'FontSize', 14, 'FontWeight', 'bold');
ylabel('Distance Error (meters)', 'FontSize', 14, 'FontWeight', 'bold');
title('Distance Error to Source over Iterations', 'FontSize', 16, 'FontWeight', 'bold');
grid on;
set(gca, 'FontSize', 12, 'Box', 'on', 'Color', 'w');
legend('Distance Error', 'Location', 'northeast');
%text(length(total_steps)*0.9, max(distance_errors)*0.9, 'Parameter: Distance Error', 'FontSize', 12, 'Color', 'k', 'HorizontalAlignment', 'right');

% 绘制 "Step Size over Time" 图表
figure;
plot(timestamp(1:length(stepSizes)), stepSizes, '-o', 'Color', [0 0.4470 0.7410], 'LineWidth', 1.5);
xlabel('Time', 'FontSize', 14, 'FontWeight', 'bold');
ylabel('Step Size', 'FontSize', 14, 'FontWeight', 'bold');
title('Step Size over Time', 'FontSize', 16, 'FontWeight', 'bold');
grid on;
set(gca, 'FontSize', 12, 'Box', 'on', 'Color', 'w'); % 纯白色背景
legend('Step Size', 'Location', 'northeast');
text(max(timestamp)*0.9, max(stepSizes)*0.9, 'Parameter: Step Size', 'FontSize', 12, 'Color', 'k', 'HorizontalAlignment', 'right');

% 绘制 "Covariance over Time" 图表
figure;
plot(timestamp(1:length(covariances)), covariances, '-^', 'Color', [0.8500 0.3250 0.0980], 'LineWidth', 1.5);
xlabel('Iteration', 'FontSize', 14, 'FontWeight', 'bold');
ylabel('Covariance', 'FontSize', 14, 'FontWeight', 'bold');
title('Covariance over Time', 'FontSize', 16, 'FontWeight', 'bold');
grid on;
set(gca, 'FontSize', 12, 'Box', 'on', 'Color', 'w');
legend('Covariance', 'Location', 'northeast');
%text(max(timestamp)*0.9, max(covariances)*0.9, 'Parameter: Covariance', 'FontSize', 12, 'Color', 'k', 'HorizontalAlignment', 'right');

% 绘制 "Time to Arrive at Source" 图表
figure;
plot(timestamp, '-s', 'Color', [0.9290 0.6940 0.1250], 'LineWidth', 2);
xlabel('Iteration', 'FontSize', 14, 'FontWeight', 'bold');
ylabel('Time (seconds)', 'FontSize', 14, 'FontWeight', 'bold');
title('Time to Arrive at Source', 'FontSize', 16, 'FontWeight', 'bold');
grid on;
set(gca, 'FontSize', 12, 'Box', 'on', 'Color', 'w');
legend('Time', 'Location', 'northeast');
%text(max(timestamp)*0.9, max(timestamp)*0.9, 'Parameter: Time', 'FontSize', 12, 'Color', 'k', 'HorizontalAlignment', 'right');

% 绘制总步数（到达源位置）的图
figure;
plot(1:i, total_steps, '-o', 'Color', [0.3010 0.7450 0.9330], 'LineWidth', 2);
xlabel('Iteration', 'FontSize', 14, 'FontWeight', 'bold');
ylabel('Total Steps', 'FontSize', 14, 'FontWeight', 'bold');
title('Total Steps to Arrive at Source', 'FontSize', 16, 'FontWeight', 'bold');
grid on;
set(gca, 'FontSize', 12, 'Box', 'on', 'Color', 'w');
legend('Total Steps', 'Location', 'northeast');
%text(i*0.9, total_steps(i)*0.9, 'Parameter: Total Steps', 'FontSize', 12, 'Color', 'k', 'HorizontalAlignment', 'right');




% After your last plot code, add this section to plot the total distance traveled

% Calculate the total distance traveled at each iteration
total_distance_traveled = cumsum(sqrt(sum(diff(P_k_store).^2, 2)));

% Plot the total distance traveled over iterations
figure;
plot(1:length(total_distance_traveled), total_distance_traveled, '-d', 'Color', [0.4660 0.6740 0.1880], 'LineWidth', 2);
xlabel('Iteration', 'FontSize', 14, 'FontWeight', 'bold');
ylabel('Total Distance Traveled (meters)', 'FontSize', 14, 'FontWeight', 'bold');
title('Total Distance Traveled to Arrive at Source', 'FontSize', 16, 'FontWeight', 'bold');
grid on;
set(gca, 'FontSize', 12, 'Box', 'on', 'Color', 'w');
legend('Total Distance Traveled', 'Location', 'northeast');
%text(length(total_distance_traveled)*0.9, max(total_distance_traveled)*0.9, 'Parameter: Total Distance Traveled', 'FontSize', 12, 'Color', 'k', 'HorizontalAlignment', 'right');
