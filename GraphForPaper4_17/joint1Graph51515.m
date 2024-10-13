
    thetaA1 = 28.2/180*pi;   % 抓取点  0.4922  
    thetaD1 = 40.7/180*pi;   % 放置点  0.7103
    
    t = [0.235,0.435,1.556,1.756,2.24];
    % 一自由度轨迹点 0.48939	0.5147	0.7679	0.76847
    g = [28.04, 29.49, 44, 44.03]/180*pi; 
    %% 第一段：五次多项式
    t0 = 0;
    theta0 = thetaA1;
    theta1 = g(1);
    v0 = 0;
    v1 = (g(2)-g(1))/0.2;
    a0 = 0;
    a1 = 0;
    figure;
    [time1, position1, velocity1, acceleration1] = quintic_trajectory(t0,t(1),theta0,theta1,v0,v1,a0,a1);
    figure;
    subplot(1, 3,1);
    plot(time1, position1,'LineWidth', 1.5);
    xlabel('Time');
    ylabel('Angel(rad)');
    ax1 = gca; % 获取第一个 subplot 的坐标轴句柄
    ax1.LineWidth = 1.5; % 设置坐标轴线宽为2
    hold on;
    
    subplot(1,3,2);
    plot(time1, velocity1,'LineWidth', 1.5);
    xlabel('Time');
    ylabel('Velocity(rad/s)');
    ax2 = gca; % 获取第一个 subplot 的坐标轴句柄
    ax2.LineWidth = 1.5; % 设置坐标轴线宽为2
    hold on;
    
    subplot(1,3,3);
    plot(time1, acceleration1,'LineWidth', 1.5);
    xlabel('Time');
    ylabel('Acceleration(rad/s^2)');
    ax3 = gca; % 获取第一个 subplot 的坐标轴句柄
    ax3.LineWidth = 1.5; % 设置坐标轴线宽为2
    hold on;
    
    %% 第二段：一次
    vc1 = v1;
    [time2, position2, velocity2, acceleration2] = curveLine(theta1,t(1),t(2),vc1);
    subplot(1, 3, 1);
    plot(time2, position2,'LineWidth', 1.5);
    hold on;
    
    subplot(1, 3, 2);
    plot(time2, velocity2,'LineWidth', 1.5);
    hold on;
    
    subplot(1, 3, 3);
    plot(time2, acceleration2,'LineWidth', 1.5);
    hold on;
    
    %% 第三段：五次
    
    theta2 = g(2);
    theta3 = g(3);
    vc2 = (g(4)-g(3))/0.2;
    [time3, position3, velocity3, acceleration3] = quintic_trajectory(t(2),t(3),theta2,theta3,vc1,vc2,a0,a1);
    subplot(1, 3, 1);
    plot(time3, position3,'LineWidth', 1.5);
    hold on;
    
    subplot(1, 3, 2);
    plot(time3, velocity3,'LineWidth', 1.5);
    hold on;
    
    subplot(1, 3, 3);
    plot(time3, acceleration3,'LineWidth', 1.5);
    hold on; 
    %% 第四段：一次
    theta3 = g(3);
    [time4, position4, velocity4, acceleration4] = curveLine(theta3,t(3),t(4),vc2);
    
    subplot(1,3, 1);
    plot(time4, position4,'LineWidth', 1.5);
    hold on;
    
    subplot(1,3,  2);
    plot(time4, velocity4,'LineWidth', 1.5);
    hold on;
    
    subplot(1,3,  3);
    plot(time4, acceleration4,'LineWidth', 1.5);
    hold on; 
    
    %% 第五段：五次
    theta4 = g(4);
    theta5 = thetaD1;
    [time5, position5, velocity5, acceleration5] = quintic_trajectory(t(4),t(5),theta4,theta5,vc2,0,a0,a1);
    subplot( 1,3, 1);
    plot(time5, position5,'LineWidth', 1.5);
    hold on;


    subplot(1,3,  2);
    plot(time5, velocity5,'LineWidth', 1.5);
    hold on;
    title("Joint 1");

    subplot(1,3,  3);
    plot(time5, acceleration5,'LineWidth', 1.5);
    hold on; 
    legend("P0P1","P1P2","P2P3","P3P4","P4P5")

