    t = [0.235,0.435,1.556,1.756,2.24];
    thetaA1 = 28.2/180*pi;   % 抓取点  0.4922  
    thetaD1 = 40.7/180*pi;   % 放置点  0.7103
    g = [28.04, 29.49, 44, 44.03]/180*pi; % 一自由度轨迹点 0.48939	0.5147	0.7679	0.76847
    
    thetaA2 = 74.7/180*pi;   % 抓取点  0.4922  
    thetaD2 = 81.9/180*pi;   % 放置点  0.7103
    g2 = [73.8689, 29.49, 44, 79.3939]/180*pi;
    
    %% 第一段五次
    t0 = 0; t02 = 0;
    theta0 = thetaA1;theta02 = thetaA2;
    theta1 = g(1);theta12 = g2(1);
    v0 = 0;v01 = 0;
    v1 = (g(2)-g(1))/0.2;v12 = (g2(2)-g2(1))/0.2;
    a0 = 0;a02 = 0;
    a1 = 0;a12 = 0;
    figure(1);
    [time1, position1, velocity1, acceleration1] = quintic_trajectory(t0,t(1),theta0,theta1,v0,v1,a0,a1);
    [time12, position12, velocity12, acceleration12] = quintic_trajectory(t02,t(1),theta02,theta12,v01,v12,a02,a12);
    
    subplot(3, 1, 1);
    plot(time1, position1,'LineWidth', 1.5);
    hold on;
    plot(time12, position12,'LineWidth', 1.5);
    hold on;
    plot(0,0);
    hold on;
    
    subplot(3, 1, 2);
    plot(time1, velocity1,'LineWidth', 1.5);
    hold on;
    plot(time12, velocity12,'LineWidth', 1.5);
    hold on;
    plot(0,0);
    hold on;
    
    subplot(3, 1, 3);
    plot(time1, acceleration1,'LineWidth', 1.5);
    hold on;
    plot(time12, acceleration12,'LineWidth', 1.5);
    hold on;
    plot(0,0);
    hold on;
    
    %% 第二段 一次
    vc1 = v1;
    [time2, position2, velocity2, acceleration2] = curveLine(theta1,t(1),t(2),vc1);
    
    subplot(3, 1, 1);
    plot(time2, position2,'LineWidth', 1.5);
    plot(0,0);
    hold on;
    plot(0,0);
    hold on;
    
    hold on;
    subplot(3, 1, 2);
    plot(time2, velocity2,'LineWidth', 1.5);
    hold on;
    plot(0,0);
    hold on;
    plot(0,0);
    hold on;
    
    subplot(3, 1, 3);
    plot(time2, acceleration2,'LineWidth', 1.5);
    hold on;
    plot(0,0);
    hold on;
    plot(0,0);
    hold on;
    
    %% 第三段：五次
    theta2 = g(2);
    theta3 = g(3);
    vc2 = (g(4)-g(3))/0.2;
    
    theta13 = 52/180*pi;   % 抓取点  0.4922  
    theta23 = 0/180*pi;   % 放置点  0.7103
    
    [time3, position3, velocity3, acceleration3] = quintic_trajectory(t(2),t(3),theta2,theta3,vc1,vc2,a0,a1);
    [time33, position33, velocity33, acceleration33] = quintic_trajectory(t(2),t(3),theta13,theta23,0,0,0,0);
    subplot(3, 1, 1);
    plot(time3, position3,'LineWidth', 1.5);
    hold on;
    plot(0,0);
    hold on;
    plot(time33, position33,'LineWidth', 1.5);
    hold on;
    
    subplot(3, 1, 2);
    plot(time3, velocity3,'LineWidth', 1.5);
    hold on;
    plot(0,0);
    hold on;
    plot(time33, velocity33,'LineWidth', 1.5);
    hold on;
    
    subplot(3, 1, 3);
    plot(time3, acceleration3,'LineWidth', 1.5);
    hold on; 
    plot(0,0);
    hold on;
    plot(time33, acceleration33,'LineWidth', 1.5);
    hold on; 
    
    %% 第四段：一次
    theta3 = g(3);
    [time4, position4, velocity4, acceleration4] = curveLine(theta3,t(3),t(4),vc2);
    
    subplot(3, 1, 1);
    plot(time4, position4,'LineWidth', 1.5);
    hold on;
    plot(0,0);
    hold on;
    plot(0,0);
    hold on;
    
    subplot(3, 1, 2);
    plot(time4, velocity4,'LineWidth', 1.5);
    hold on;
    plot(0,0);
    hold on;
    plot(0,0);
    hold on;
    
    subplot(3, 1, 3);
    plot(time4, acceleration4,'LineWidth', 1.5);
    hold on; 
    plot(0,0);
    hold on;
    plot(0,0);
    hold on;
    
    %% 第五段：五次
    
    theta4 = g(4);theta42 = g2(4);
    theta5 = thetaD1;theta52 = thetaD2;
    [time5, position5, velocity5, acceleration5] = quintic_trajectory(t(4),t(5),theta4,theta5,vc2,0,a0,a1);
    [time52, position52, velocity52, acceleration52] = quintic_trajectory(t(4),t(5),theta42,theta52,0,0,0,0);
    subplot(3, 1, 1);
    xlabel('Time');
    ylabel('Position');
    
    plot(time5, position5,'LineWidth', 1.5);
    hold on;
    plot(time52, position52,'LineWidth', 1.5);
    hold on;
    plot(0,0);
    hold on;
    
    subplot(3, 1, 2);
    xlabel('Time');
    ylabel('Velocity');
    
    plot(time5, velocity5,'LineWidth', 1.5);
    hold on;
    plot(time52, velocity52,'LineWidth', 1.5);
    hold on;
    plot(0,0);
    hold on;
    
    subplot(3, 1, 3);
    xlabel('Time');
    ylabel('Acceleration');
    
    plot(time5, acceleration5,'LineWidth', 1.5);
    hold on; 
    plot(time52, acceleration52,'LineWidth', 1.5);
    hold on;
    plot(0,0);
    hold on;
    legend('axes 1', 'axes 2','axes 3'); % Add legend to differentiate the trajectories

%% 运动学正解图像
% 第一段
tt1 = linspace(0,t(1),50);
coefficients11 = quintic_coeff(t0,t(1),theta0,theta1,v0,v1,a0,a1);
theta1t1= coefficients11(1) + coefficients11(2)*tt1 + coefficients11(3)*tt1.^2 + coefficients11(4)*tt1.^3 + coefficients11(5)*tt1.^4 + coefficients11(6)*tt1.^5;
coefficients12 = quintic_coeff(t02,t(1),theta02,theta12,v01,v12,a02,a12);
theta1t12= coefficients12(1) + coefficients12(2)*tt1 + coefficients12(3)*tt1.^2 + coefficients12(4)*tt1.^3 + coefficients12(5)*tt1.^4 + coefficients12(6)*tt1.^5;

% 第二段
tv1=linspace(0,0.2,50);
theta1c1=g(1)+vc1*tv1;

% 第三段
tt2 = linspace(0,t(3)-t(2),50);
coefficients13 = quintic_coeff(t(2),t(3),theta2,theta3,vc1,vc2,a0,a1);
theta1t2 = coefficients13(1) + coefficients13(2)*tt2 + coefficients13(3)*tt2.^2 + coefficients13(4)*tt2.^3 + coefficients13(5)*tt2.^4 + coefficients13(6)*tt2.^5;

xt1 = zeros(1,50);yt1 = zeros(1,50);
xc1 = zeros(1,50);yc1 = zeros(1,50);
xt2 = zeros(1,50);yt2 = zeros(1,50);
xc2 = zeros(1,50);yc2 = zeros(1,50);
xt3 = zeros(1,50);yt3 = zeros(1,50);

for i = 1:50
    [xt1(i),yt1(i)] = ForwardKimemicDouble(theta1t1(i),theta1t12(i));%二自由度正运动学求解
    [xc1(i),yc1(i)] = ForwardKimemicSingle(theta1c1(i));%单自由度正运动学求解
    [xt2(i),yt2(i)] = ForwardKimemicSingle(theta1t2(i));%单自由度正运动学求解
% %     [xc2(i),yc2(i)] = ForwardKimemicSingle(theta1c2(i));%单自由度正运动学求解
%     [xt3(i),yt3(i)] = ForwardKimemicDouble(theta1t3(i),theta2t3(i));%二自由度正运动学求解
end
figure(2);
plot(xt1+0.12,yt1-0.09,'color',[227/255,118/255,28/255],'linewidth',3);%二自由度末端轨迹
hold on;
plot(xc1+0.12,yc1-0.09,'color',[0/255,176/255,80/255],'linewidth',3);%单自由度末端轨迹
hold on;
plot(xt2+0.12,yt2-0.09,'color',[0/255,176/255,240/255],'linewidth',3);%单自由度末端轨迹


%%
function [coefficients] = quintic_coeff(t0, t1, p0, p1, v0, v1, a0, a1)
    % 增加插入点来消除荣格现象的五次多项式轨迹规划
    % 设置时间段
    dt = t1 - t0;
    if dt >1000
        coefficients = [0 0 0 0 0 0];
        return;
    end
    % 计算多项式系数矩阵
    A1 = [1, t0, t0^2, t0^3, t0^4, t0^5;
         0, 1, 2*t0, 3*t0^2, 4*t0^3, 5*t0^4;
         0, 0, 2, 6*t0, 12*t0^2, 20*t0^3;
         1, t1, t1^2, t1^3, t1^4, t1^5;
         0, 1, 2*t1, 3*t1^2, 4*t1^3, 5*t1^4;
         0, 0, 2, 6*t1, 12*t1^2, 20*t1^3];

    A2 = [t0^5, t0^4, t0^3, t0^2, t0, 1;
         5*t0^4, 4*t0^3, 3*t0^2, 2*t0, 1, 0;
         20*t0^3, 12*t0^2, 6*t0, 2, 0, 0;
         t1^5, t1^4, t1^3, t1^2, t1, 1;
         5*t1^4, 4*t1^3, 3*t1^2, 2*t1, 1, 0;
         20*t1^3, 12*t1^2, 6*t1, 2, 0, 0];
    
    b = [p0; v0; a0; p1; v1; a1];
    
    A = A1;
    % 解线性方程组（考虑奇异性）
    if rank(A) < size(A, 2)
        coefficients = pinv(A) * b;
    else
        coefficients = A \ b;
    end
end
%% 二自由度
function [ Fx,Fy ] = ForwardKimemicDouble( theta1,theta2 )
    l1 = 0.6;%AB杆长度
    l2 = 0.35;%DE杆长度
    l3 = 0.2;%BC杆长度
    l4 = 0.608;%CD杆长度
    l5 = 0.4;%CF杆长度
    l6 = 0.455;%BD杆长度
    l7 = 0.196;%AE杆长度
    
    Ax = 0;%A点横坐标
    Ay = 0;%A点纵坐标
    Ex = 0;%E点横坐标
    Ey = -l7;%E点纵坐标
    
    Bx = l1*cos(theta1);%B点横坐标
    By = l1*sin(theta1);%B点横坐标
    Dx = l2*cos(theta2);%D点横坐标
    Dy = l2*sin(theta2)-l7;%D点横坐标

    n1n2 = sqrt((Bx-Dx).^2+(By-Dy).^2);
    angle_n1 = acos((n1n2.^2+l3.^2-l4.^2)/(2*n1n2*l3));
    Cx = Bx + cos(azimuth(Bx,By,Dx,Dy) + angle_n1) * l3;%C点横坐标
    Cy = By + sin(azimuth(Bx,By,Dx,Dy) + angle_n1) * l3;%C点纵坐标
    
    Angle_BC = azimuth(Bx,By,Cx,Cy);%BC杆方位角
    Angle_CF = Angle_BC;%CF杆方位角
    
    Fx = Cx + cos(Angle_CF) * l5;
    Fy = Cy + sin(Angle_CF) * l5;

end

%% 
function [ Fx,Fy ] = ForwardKimemicSingle( theta1 )
    l1 = 0.6;%AB杆长度
    l2 = 0.35;%DE杆长度
    l3 = 0.2;%BC杆长度
    l4 = 0.608;%CD杆长度
    l5 = 0.4;%CF杆长度
    l6 = 0.455;%BD杆长度
    l7 = 0.196;%AE杆长度
    
    Ax = 0;%A点横坐标
    Ay = 0;%A点纵坐标
    Ex = 0;%E点横坐标
    Ey = -l7;%E点纵坐标
    
    Bx = l1*cos(theta1);%B点横坐标
    By = l1*sin(theta1);%B点横坐标
    
    n1n2 = sqrt((Bx-Ex).^2+(By-Ey).^2);
    angle_n1 = acos((n1n2.^2+l6.^2-l2.^2)/(2*n1n2*l6));
    Dx = Bx + cos(azimuth(Bx,By,Ex,Ey) - angle_n1) * l6;%D点横坐标
    Dy = By + sin(azimuth(Bx,By,Ex,Ey) - angle_n1) * l6;%D点纵坐标
    
    Angle_BD = azimuth(Bx,By,Dx,Dy);%DC杆方位角
    angle_CBD = acos((l3.^2+l6.^2-l4.^2)/(2*l3*l6));
    Angle_BC = Angle_BD + angle_CBD;

    Cx = Bx + cos(Angle_BC) * l3;%C点横坐标
    Cy = By + sin(Angle_BC) * l3;%C点纵坐标
    Angle_CF = Angle_BC;
    
    Fx = Cx + cos(Angle_CF) * l5;%F点纵坐标
    Fy = Cy + sin(Angle_CF) * l5;%F点纵坐标
end
