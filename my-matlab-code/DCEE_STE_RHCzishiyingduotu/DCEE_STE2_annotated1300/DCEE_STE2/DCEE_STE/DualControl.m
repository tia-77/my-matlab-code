%该代码是在自适应步长+DCEE+RHC结合的基础上输出各种目的图
close all
clear all
clc
% 
% function P_k_store=DualControl(HL);

addpath('Functions')

UAVVel=2; %UAV velocity
sampleTime=10;
bLim=700;
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
%stepSize = 0.5; % How far to move


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


% 初始化变量
stepSizes = [];
covariances = [];

% 初始化HL记录数组
HL_values = zeros(1, 100);  % 假设最大迭代次数为100次
total_steps = 0;  % 初始化总步数
total_steps_array = zeros(1, 100);  % 初始化存储步数的数组
first_arrival = false;  % 标志是否已经到达源头

timestamp(1)=0;
D=[];
sampleHistory=P_k(1:2);
% step_no=0;%calculat the No. of steps;
% distance=0;%calculat the distance from start position;



% 定义一个接近源头的阈值和微调步长
close_threshold = 0.5;
fine_tune_step = 0.1;
total_distance = 0; % 用于累计总距离



  total_distances = zeros(1, 100);  % 假设最大迭代次数为100次
  distance_errors = zeros(1, 100);  % 假设最大迭代次数为100次
first_arrival = false;  % 标志是否已经到达源头
% 在循环外部初始化 total_steps
   total_steps = 0;
    close_threshold = 0.1; %距离污染源的距离
for i = 1:100

  
 
distance_to_source = norm([s.x, s.y, s.z] - P_k);

if distance_to_source < close_threshold && ~first_arrival
        first_arrival = true;
    end

% 更新位置后的时间计算
time_to_move = floor(norm(P_k - P_k_store(end, :)) / UAVVel) + sampleTime;
bLim = bLim - time_to_move;
if bLim <= 0
    break;
end
% 仅当尚未到达源头时计算步数
    if ~first_arrival
        if i > 1
            step_distance = norm(P_k - P_k_store(end,:));
            total_steps = total_steps + step_distance;
        total_distance = total_distance + step_distance;
        else
            total_steps = total_steps + stepsize;
        end
        total_steps_array(i) = total_steps;
     total_distances(i) = total_distance;
    else
        % 一旦到达源头，步数不再增加
        total_steps_array(i) = total_steps_array(i-1);
% 一旦到达源头，停止累积移动距离
        total_distances(i) = total_distances(i-1);
    end

% 更新时间戳
    timestamp(i+1) = timestamp(i) + time_to_move;


% 计算步长：使用当前位置与前一位置的差值来动态调整步长
    if i > 1
        step_distance = norm(P_k - P_k_store(end,:));
        total_steps = total_steps + step_distance;
    else
        total_steps = total_steps + stepsize; % 首次迭代时，直接使用初始步长
    end
    

 % 存储当前迭代的总步数
    total_steps_array(i) = total_steps;
      % 在此处更新 P_k_store 中的 P_k
    P_k_store = [P_k_store; P_k];
        
       
    % 模拟数据
    Dsim = PredictMeasurement(s, pos, m.thresh);
    D(i) = Dsim;

    % 更新粒子滤波器
    thetaPrev = theta;
    [theta, Wpnorm] = UpdatePFPlume(D, theta, Wpnorm, pos, P_k_store, m, N, PF_Memory, domain);


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
    
    % % horizon length
     %HL = 3;

    % Calculation step size, inversely proportional to the variance
    Covar = cov(theta.x, theta.y);
    variance = sqrt(Covar(1,1) + Covar(2,2));
    k_step = 1; 
    minStep = 0.01;
    maxStep = 0.5;
    threshold = 0.1; % 根据实际情况调整阈值
    % stepSize = adjustStepSize(variance, minStep, maxStep, m.thresh, k_step);
   stepSize = 0.5;


 % 更新步长和协方差
    stepSizes = [stepSizes; stepSize];
    covariances = [covariances; variance];



     % Calculate the prediction length
    k_hl = 1.5; 
    minLength = 1; % Minimum prediction length
    maxLength = 20; % Maximum prediction length
    HL = adjustPredictionHorizon(variance, minLength, maxLength, k_hl);
% 记录每次迭代的HL值
HL_values(i) = HL;


    % calculate the next position using RHC
    [~, pos] = RHC(HL,1,stepSize,pos,pos,theta,Wpnorm,D,m,s,M,MM,i,P_k_store,N,PF_Memory,domain,xmin,xmax,ymin,ymax,zmin,zmax);
    
    P_k = [pos.x_matrix pos.y_matrix pos.z_matrix];
    sampleHistory=[sampleHistory;P_k(1:2)];

    % 计算当前位置与源位置的距离误差
    distance_error = norm([s.x, s.y, s.z] - P_k);
    
    % 存储当前的误差
    distance_errors(i) = distance_error;

    % 插入的代码部分：计算步距并累计总距离
% 计算步距并累计总距离
if i > 1
    step_distance = norm(P_k - P_k_store(i-1,:));
    total_distance = total_distance + step_distance;
end

% 存储每一步的累积总距离
total_distances(i) = total_distance;
% % 添加一个变量来存储总步数（到达源位置）
% total_steps = 0;






    % 存储当前位置
    P_k_store = [P_k_store; P_k];

    
    % 更新步长和时间戳
    move_time=floor(norm(P_k-P_k_store(i,:))/UAVVel)+sampleTime;
    bLim=bLim-move_time;
    if bLim <= 0
        break;
    end
    timestamp(i+1)=timestamp(i)+move_time;

    % 存储每一步的累积总距离
    total_distances(i) = total_distance;

    % 继续更新位置
    P_k_store = [P_k_store; P_k];
 

% 在更新完 bLim 后，更新 timestamp
timestamp(i+1) = timestamp(i) + move_time;
    
    Covar = cov(theta.x,theta.y);
    Spread = sqrt(Covar(1,1)+Covar(2,2))

    %  if norm(P_k(1:2)-[s.x s.y])==0
    %     break
    % end
    
    if Spread<0 % 3.5? 4? 5 default
        break
    end
    
end


% 在主循环结束后，计算最终位置与源头的距离误差
final_position = P_k_store(end, :);  % 机器人最终位置
distance_error = norm([s.x, s.y, s.z] - final_position);  % 计算误差

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
% Function definitiono
% function stepSize = adjustStepSize(variance, minStep, maxStep, threshold, k)
%     epsilon = 1e-6;
%     variance = max(variance, epsilon);
%     stepSize = k / variance;
%     stepSize = max(minStep, min(maxStep, stepSize));
% end
% 
% function HL = adjustPredictionHorizon(covariance, minHorizon, maxHorizon, k)
%     epsilon = 1e-6;
%     covariance = max(covariance, epsilon);
%     HL = k * covariance;
%     HL = round(max(minHorizon, min(maxHorizon, HL)));
%     HL = max(HL, 1); % Make sure HL is at least 1
% end


% %Function definition 原始变方差平方根
% 
% function stepSize = adjustStepSize(variance, minStep, maxStep, threshold, k)
%     epsilon = 1e-6;
%     variance = max(variance, epsilon);
%     stepSize = k / sqrt(variance);
%     stepSize = max(minStep, min(maxStep, stepSize));
% end
% 
% function HL = adjustPredictionHorizon(covariance, minHorizon, maxHorizon, k)
%     epsilon = 1e-6;
%     covariance = max(covariance, epsilon);
%     HL = k * sqrt(covariance);
%     HL = round(max(minHorizon, min(maxHorizon, HL)));
%     HL = max(HL, 1); % Make sure HL is at least 1
% end



% function stepSize = adjustStepSize(variance, minStep, maxStep, threshold, k)
%     epsilon = 1e-6;
%     variance = max(variance, epsilon);
% 
%     if variance < threshold
%         stepSize = minStep; % 如果方差大于阈值，使用最小步长
%     else
%         stepSize = k * variance; % 否则根据方差动态调整步长
%     end
%     stepSize = max(minStep, min(maxStep, stepSize));
% end
% 


function HL = adjustPredictionHorizon(covariance, minHorizon, maxHorizon, k)
    epsilon = 1e-6;
    covariance = max(covariance, epsilon);

    if covariance < 1 % 设定一个阈值
        HL = minHorizon; % 如果协方差大于阈值，使用最小预测长度
    else
        HL = k * covariance; % 否则根据协方差动态调整预测长度
    end
    HL = round(max(minHorizon, min(maxHorizon, HL)));
    HL = max(HL, 1); % 确保 HL 至少为 1
end



% 绘制误差随迭代次数的变化曲线
figure;
plot(1:i, distance_errors(1:i), '-o', 'Color', [0.6350 0.0780 0.1840], 'LineWidth', 2);
xlabel('Iteration', 'FontSize', 14, 'FontWeight', 'bold');
ylabel('Distance Error (meters)', 'FontSize', 14, 'FontWeight', 'bold');
title('Distance Error to Source over Iterations', 'FontSize', 16, 'FontWeight', 'bold');
grid on;
set(gca, 'FontSize', 12, 'Box', 'on', 'Color', 'w');
legend('Distance Error', 'Location', 'northeast');
%text(i*0.9, distance_errors(i)*0.9, 'Parameter: Distance Error', 'FontSize', 12, 'Color', 'k', 'HorizontalAlignment', 'right');

% % 绘制 "Step Size over Time" 图表
% figure;
% plot(timestamp(1:length(stepSizes)), stepSizes, '-o', 'Color', [0 0.4470 0.7410], 'LineWidth', 1.5);
% xlabel('Time', 'FontSize', 14, 'FontWeight', 'bold');
% ylabel('Step Size', 'FontSize', 14, 'FontWeight', 'bold');
% title('Step Size over Time', 'FontSize', 16, 'FontWeight', 'bold');
% grid on;
% set(gca, 'FontSize', 12, 'Box', 'on', 'Color', 'w'); % 纯白色背景
% legend('Step Size', 'Location', 'northeast');
% % 右上角标注参数量名称
% text(max(timestamp)*0.9, max(stepSizes)*0.9, 'Parameter: Step Size', 'FontSize', 12, 'Color', 'k', 'HorizontalAlignment', 'right');

% 绘制 "Covariance over Time" 图表
figure;
plot(timestamp(1:length(covariances)), covariances, '-^', 'Color', [0.8500 0.3250 0.0980], 'LineWidth', 1.5);
xlabel('Iteration', 'FontSize', 14, 'FontWeight', 'bold');
ylabel('Covariance', 'FontSize', 14, 'FontWeight', 'bold');
title('Covariance over Iteration', 'FontSize', 16, 'FontWeight', 'bold');
grid on;
set(gca, 'FontSize', 12, 'Box', 'on', 'Color', 'w');
legend('Covariance', 'Location', 'northeast');
%text(max(timestamp)*0.9, max(covariances)*0.9, 'Parameter: Covariance', 'FontSize', 12, 'Color', 'k', 'HorizontalAlignment', 'right');

% 绘制 "Time to Arrive at Source" 图表
figure;
plot(1:length(timestamp), timestamp, '-s', 'Color', [0.9290 0.6940 0.1250], 'LineWidth', 2);
xlabel('Iteration', 'FontSize', 14, 'FontWeight', 'bold');
ylabel('Time (seconds)', 'FontSize', 14, 'FontWeight', 'bold');
title('Time to Arrive at Source', 'FontSize', 16, 'FontWeight', 'bold');
grid on;
set(gca, 'FontSize', 12, 'Box', 'on', 'Color', 'w');
legend('Time', 'Location', 'northeast');



figure;
plot(1:i, total_distances(1:i), '-d', 'Color', [0.4940 0.1840 0.5560], 'LineWidth', 2);
xlabel('Iteration', 'FontSize', 14, 'FontWeight', 'bold');
ylabel('Total Distance (meters)', 'FontSize', 14, 'FontWeight', 'bold');
title('Total Distance Traveled', 'FontSize', 16, 'FontWeight', 'bold');
grid on;
set(gca, 'FontSize', 12, 'Box', 'on', 'Color', 'w');
legend('Distance', 'Location', 'northeast');
%text(i*0.9, total_distances(i)*0.9, 'Parameter: Distance', 'FontSize', 12, 'Color', 'k', 'HorizontalAlignment', 'right');


figure;
plot(1:i, total_steps_array(1:i), '-o', 'Color', [0.3010 0.7450 0.9330], 'LineWidth', 2);
xlabel('Iteration', 'FontSize', 14, 'FontWeight', 'bold');
ylabel('Total Steps', 'FontSize', 14, 'FontWeight', 'bold');
title('Total Steps to Arrive at Source', 'FontSize', 16, 'FontWeight', 'bold');
grid on;
set(gca, 'FontSize', 12, 'Box', 'on', 'Color', 'w');
legend('Total Steps', 'Location', 'northeast');


% 在主循环结束后，绘制HL随迭代次数的变化图
figure;
plot(1:i, HL_values(1:i), '-o', 'Color', [0.4940 0.1840 0.5560], 'LineWidth', 2);
xlabel('Iteration', 'FontSize', 14, 'FontWeight', 'bold');
ylabel('Prediction Horizon Length (HL)', 'FontSize', 14, 'FontWeight', 'bold');
title('Prediction Horizon Length (HL) over Iterations', 'FontSize', 16, 'FontWeight', 'bold');
grid on;
set(gca, 'FontSize', 12, 'Box', 'on', 'Color', 'w');
legend('HL', 'Location', 'northeast');
%text(i*0.9, HL_values(i)*0.9, 'Parameter: HL', 'FontSize', 12, 'Color', 'k', 'HorizontalAlignment', 'right');
