function AngleVelocityAcc
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
    [time1, position1, velocity1, acceleration1,jerk1] = quintic_trajectory(t0,t(1),theta0,theta1,v0,v1,a0,a1);
    [time12, position12, velocity12, acceleration12,jerk12] = quintic_trajectory(t02,t(1),theta02,theta12,v01,v12,a02,a12);
    
    plot(time1, position1,'LineWidth', 2,"Color","0.00,0.45,0.74");
    hold on;
    plot(time12, position12,'LineWidth', 2,"Color","0.85,0.33,0.10");
    hold on;
    plot(0,0);
    hold on;

    figure(2);
    plot(time1, velocity1,'LineWidth', 2,"Color","0.00,0.45,0.74");
    hold on;
    plot(time12, velocity12,'LineWidth', 2,"Color","0.85,0.33,0.10");
    hold on;
    plot(0,0);
    hold on;
    
    figure(3);
    plot(time1, acceleration1,'LineWidth', 2,"Color","0.00,0.45,0.74");
    hold on;
    plot(time12, acceleration12,'LineWidth', 2,"Color","0.85,0.33,0.10");
    hold on;
    plot(0,0);
    hold on;
    
    figure(4);
    plot(time1, jerk1,'LineWidth', 2,"Color","0.00,0.45,0.74");
    hold on;
    plot(time12, jerk12,'LineWidth', 2,"Color","0.85,0.33,0.10")
    plot(0,0);
    hold on;

    %% 第二段 一次
    vc1 = v1;
    [time2, position2, velocity2, acceleration2] = curveLine(theta1,t(1),t(2),vc1);
    [position22] = coordinate1(position2);
    
    figure(1);
    plot(time2, position2,'LineWidth', 2,"Color","0.00,0.45,0.74");
    plot(time2, position22,'LineWidth', 2,"Color","0.85,0.33,0.10");
    
    plot(0,0);
    hold on;
    plot(0,0);
    hold on;
    
    figure(2);
    plot(time2, velocity2,'LineWidth', 2,"Color","0.00,0.45,0.74");
    hold on;
    plot(0,0);
    hold on;
    plot(0,0);
    hold on;
    
    figure(3);
    plot(time2, acceleration2,'LineWidth', 2,"Color","0.00,0.45,0.74");
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
    
    [time3, position3, velocity3, acceleration3,jerk3] = quintic_trajectory(t(2),t(3),theta2,theta3,vc1,vc2,a0,a1);
    [time33, position33, velocity33, acceleration33,jerk33] = quintic_trajectory(t(2),t(3),theta13,theta23,0,0,0,0);
    figure(1);
    plot(time3, position3,'LineWidth', 2,"Color","0.00,0.45,0.74");
    hold on;
    plot(0,0);
    hold on;
    plot(time33, position33,'LineWidth', 2,"Color","0.93,0.69,0.13");
    hold on;
    
    figure(2);
    plot(time3, velocity3,'LineWidth', 2,"Color","0.00,0.45,0.74");
    hold on;
    plot(0,0);
    hold on;
    plot(time33, velocity33,'LineWidth', 2,"Color","0.93,0.69,0.13");
    hold on;
    
    figure(3);
    plot(time3, acceleration3,'LineWidth', 2,"Color","0.00,0.45,0.74");
    hold on; 
    plot(0,0);
    hold on;
    plot(time33, acceleration33,'LineWidth', 2,"Color","0.93,0.69,0.13");
    hold on; 
    
    figure(4);
    plot(time3, jerk3,'LineWidth', 2,"Color","0.00,0.45,0.74");
    hold on; 
    plot(0,0);
    hold on;
    plot(time33, jerk33,'LineWidth', 2,"Color","0.93,0.69,0.13");
    hold on; 

    %% 第四段：一次
    theta3 = g(3);
    [time4, position4, velocity4, acceleration4] = curveLine(theta3,t(3),t(4),vc2);
    [position42] = coordinate1(position4);
    figure(1);
    plot(time4, position4,'LineWidth', 2,"Color","0.00,0.45,0.74");
    plot(time4, position42,'LineWidth', 2,"Color","0.85,0.33,0.10");
    hold on;
    plot(0,0);
    hold on;
    plot(0,0);
    hold on;
    
    figure(2);
    plot(time4, velocity4,'LineWidth', 2,"Color","0.00,0.45,0.74");
    hold on;
    plot(0,0);
    hold on;
    plot(0,0);
    hold on;
    
    figure(3);
    plot(time4, acceleration4,'LineWidth', 2,"Color","0.00,0.45,0.74");
    hold on; 
    plot(0,0);
    hold on;
    plot(0,0);
    hold on;
    
    %% 第五段：五次
    
    theta4 = g(4);theta42 = g2(4);
    theta5 = thetaD1;theta52 = thetaD2;
    [time5, position5, velocity5, acceleration5,jerk5] = quintic_trajectory(t(4),t(5),theta4,theta5,vc2,0,a0,a1);
    [time52, position52, velocity52, acceleration52,jerk52] = quintic_trajectory(t(4),t(5),theta42,theta52,0,0,0,0);
    
    figure(1);
    
    xlabel('Time(s)');
    ylabel('Angle(rad)');
    
    plot(time5, position5,'LineWidth', 2,"Color","0.00,0.45,0.74");
    hold on;
    plot(time52, position52,'LineWidth', 2,"Color","0.85,0.33,0.10");
    hold on;
    plot(0,0);
    hold on;
    legend({'Axis 1', 'Axis 2','Axis 3'},'Location','southeast');
    hold on;

    figure(2);
    xlabel('Time(s)');
    ylabel('Velocity(rad/s)');

    plot(time5, velocity5,'LineWidth', 2,"Color","0.00,0.45,0.74");
    hold on;
    plot(time52, velocity52,'LineWidth', 2,"Color","0.85,0.33,0.10");
    hold on;
    plot(0,0);
    hold on;
    legend({'Axis 1', 'Axis 2','Axis 3'},'Location','southeast');
    hold on;

    figure(3);
    xlabel('Time(s)');
    ylabel('Acceleration(rad/s^2)');
    
    plot(time5, acceleration5,'LineWidth', 2,"Color","0.00,0.45,0.74");
    hold on; 
    plot(time52, acceleration52,'LineWidth', 2,"Color","0.85,0.33,0.10");
    hold on;
    plot(0,0);
    hold on;
    legend({'Axis 1', 'Axis 2','Axis 3'},'Location','southeast');
    hold on;

    figure(4);
    xlabel('Time(s)');
    ylabel('Jerk(rad/s^3)');
    
    plot(time5, jerk5,'LineWidth', 2,"Color","0.00,0.45,0.74");
    hold on; 
    plot(time52, jerk52,'LineWidth', 2,"Color","0.85,0.33,0.10");
    hold on;
    plot(0,0);
    hold on;
    % 绘制竖直虚线
    figure(1);
    plot([t(1), t(1)], [-2, 2], '--k'); % '--k'表示黑色虚线
    hold on;
    plot([t(2), t(2)], [-2, 2], '--k'); % '--k'表示黑色虚线
    hold on;
    plot([t(3), t(3)], [-2, 2], '--k'); % '--k'表示黑色虚线
    hold on;
    plot([t(4), t(4)], [-2, 2], '--k'); % '--k'表示黑色虚线
    hold on;
    legend({'Axis 1', 'Axis 2','Axis 3'},'Location','southeast'); % Add legend to differentiate the trajectories
end
