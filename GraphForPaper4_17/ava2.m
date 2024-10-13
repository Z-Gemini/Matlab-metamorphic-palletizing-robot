function ava2
%     t = [0.8909 1.3088 1.4069 1.7143 2.30];
    t = [0.235,0.435,1.556-0.14,1.756-0.14,2.380-0.14];
    thetaA1 = deg2rad(28.2);   % 抓取点  0.4922  
    thetaD1 = deg2rad(40.7);   % 放置点  0.7103
%     g = [0.7473,0.7692,0.7773,0.7934]; % 一自由度轨迹点 0.48939	0.5147	0.7679	0.76847
    g = deg2rad([28.04,29.49,44,44.03]);
    thetaA2 = deg2rad(74.7); % 抓取点  0.4922  
    thetaD2 = deg2rad(81.9);   % 放置点  0.7103
%     g2 = [1.3785,1.3859,1.3887,1.3941];
    g2=[1.2893, 1.2982,1.3855,1.3857];
    %% 第一段五次
    t0 = 0; t02 = 0;
    theta0 = thetaA1;theta02 = thetaA2;
    theta1 = g(1);theta12 = g2(1);
    v0 = 0;v01 = 0;
    v1 = (g(2)-g(1))/(t(2)-t(1));v12 = (g2(2)-g2(1))/(t(2)-t(1));
    a0 = 0;a02 = 0;
    a1 = 0;a12 = 0;
    figure(1);
    [time1, position1, velocity1, acceleration1,jerk1] = quintic_trajectory(t0,t(1),theta0,theta1,v0,v1,a0,a1);
    [time12, position12, velocity12, acceleration12,jerk12] = quintic_trajectory(t02,t(1),theta02,theta12,v01,v12,a02,a12);
    
    position1Plot = rad2deg(position1);
    plot(time1, position1Plot,'LineWidth', 2,"Color",[0,139/255,139/255]);
    hold on;
    position12Plot = rad2deg(position12);
    plot(time12, position12Plot,'LineWidth', 2,"Color",[220/255,20/255,60/255]);
    hold on;
    plot(0,0);
    hold on;

    figure(2);
    plotvelocity1=rad2deg(velocity1);
    plot(time1, plotvelocity1,'LineWidth', 2,"Color",[0,139/255,139/255]);
    hold on;
    plotvelocity12 = rad2deg(velocity12);
    plot(time12, plotvelocity12,'LineWidth', 2,"Color",[220/255,20/255,60/255]);
    hold on;
    plot(0,0);
    hold on;
    
    figure(3);
    plotacceleration1 = rad2deg(acceleration1);
    plot(time1, plotacceleration1,'LineWidth', 2,"Color",[0,139/255,139/255]);
    hold on;
    plotacceleration12 = rad2deg(acceleration12);
    plot(time12, plotacceleration12,'LineWidth', 2,"Color",[220/255,20/255,60/255]);
    hold on;
    plot(0,0);
    hold on;
    
    figure(4);
    plotjerk1  =rad2deg(jerk1);
    plot(time1, plotjerk1,'LineWidth', 2,"Color",[0,139/255,139/255]);
    hold on;
    plotjerk12  =rad2deg(jerk12);
    plot(time12, plotjerk12,'LineWidth', 2,"Color",[220/255,20/255,60/255])
    plot(0,0);
    hold on;

    %% 第二段 一次
    vc1 = v1;
    [time2, position2, velocity2, acceleration2] = curveLine(theta1,t(1),t(2),vc1);
    [position22] = coordinate1(position2);
    
    figure(1);
    plotposition2 = rad2deg(position2);
    plot(time2, plotposition2,'LineWidth', 2,"Color",[0,139/255,139/255]);
    plotposition22 = rad2deg(position22);   
    plot(time2, plotposition22,'LineWidth', 2,"Color",[220/255,20/255,60/255]);
    
    plot(0,0);
    hold on;
    plot(0,0);
    hold on;
    
    figure(2);
    plotvelocity2 = rad2deg(velocity2);
    plot(time2, plotvelocity2,'LineWidth', 2,"Color",[0,139/255,139/255]);
    hold on;
    plot(0,0);
    hold on;
    plot(0,0);
    hold on;
    
    figure(3);
    plotacceleration2 = rad2deg(acceleration2);
    plot(time2, plotacceleration2,'LineWidth', 2,"Color",[0,139/255,139/255]);
    hold on;
    plot(0,0);
    hold on;
    plot(0,0);
    hold on;
    
    %% 第三段：五次
    theta2 = g(2);
    theta3 = g(3);
    vc2 = (g(4)-g(3))/(t(4)-t(3));
    
    theta13 = 52/180*pi;   % 抓取点  0.4922  
    theta23 = 0/180*pi;   % 放置点  0.7103
    
    [time3, position3, velocity3, acceleration3,jerk3] = quintic_trajectory(t(2),t(3),theta2,theta3,vc1,vc2,a0,a1);
    [time33, position33, velocity33, acceleration33,jerk33] = quintic_trajectory(t(2),t(3),theta13,theta23,0,0,0,0);
    figure(1);
    plotposition3 = rad2deg(position3);
    plot(time3, plotposition3,'LineWidth', 2,"Color",[0,139/255,139/255]);
    hold on;
    plot(0,0);
    hold on;
    plotposition33 = rad2deg(position33);
    plot(time33, plotposition33,'LineWidth', 2,"Color",[218/255,165/255,32/255]);
    hold on;
    
    figure(2);
    plotvelocity3 = rad2deg(velocity3);
    plot(time3, plotvelocity3,'LineWidth', 2,"Color",[0,139/255,139/255]);
    hold on;
    plot(0,0);
    hold on;
    plotvelocity33 = rad2deg(velocity33);
    plot(time33, plotvelocity33,'LineWidth', 2,"Color",[218/255,165/255,32/255]);
    hold on;
    
    figure(3);
    plotacceleration3 = rad2deg(acceleration3);
    plot(time3, plotacceleration3,'LineWidth', 2,"Color",[0,139/255,139/255]);
    hold on; 
    plot(0,0);
    hold on;
    plotacceleration33 = rad2deg(acceleration33);
    plot(time33, plotacceleration33,'LineWidth', 2,"Color",[218/255,165/255,32/255]);
    hold on; 
    
    figure(4);
    plotjerk3 = rad2deg(jerk3);
    plot(time3, plotjerk3,'LineWidth', 2,"Color",[0,139/255,139/255]);
    hold on; 
    plot(0,0);
    hold on;
    plotjerk33 = rad2deg(jerk33);
    plot(time33, plotjerk33,'LineWidth', 2,"Color",[218/255,165/255,32/255]);
    hold on; 

    %% 第四段：一次
    theta3 = g(3);
    [time4, position4, velocity4, acceleration4] = curveLine(theta3,t(3),t(4),vc2);
    [position42] = coordinate1(position4);
    figure(1);
    plotposition4 = rad2deg(position4);
    plot(time4, plotposition4,'LineWidth', 2,"Color",[0,139/255,139/255]);
    plotposition42 = rad2deg(position42);
    plot(time4, plotposition42,'LineWidth', 2,"Color",[220/255,20/255,60/255]);
    hold on;
    plot(0,0);
    hold on;
    plot(0,0);
    hold on;
    
    figure(2);
    plotvelocity4 = rad2deg(velocity4);
    plot(time4, plotvelocity4,'LineWidth', 2,"Color",[0,139/255,139/255]);
    hold on;
    plot(0,0);
    hold on;
    plot(0,0);
    hold on;
    
    figure(3);
    plotacceleration4 = rad2deg(acceleration4);
    plot(time4, plotacceleration4,'LineWidth', 2,"Color",[0,139/255,139/255]);
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
    ylabel('Angle(°)');
    plotposition5 = rad2deg(position5);
    plot(time5, plotposition5,'LineWidth', 2,"Color",[0,139/255,139/255]);
    hold on;
    plotposition52 = rad2deg(position52);
    plot(time52, plotposition52,'LineWidth', 2,"Color",[220/255,20/255,60/255]);
    hold on;
    plot(0,0);
    hold on;
    legend({'Axis 1', 'Axis 2','Axis 3'},'Location','southeast');
    hold on;

    figure(2);
    xlabel('Time(s)');
    ylabel('Velocity(°/s)');
   
    plotvelocity5 = rad2deg(velocity5);
    plot(time5, plotvelocity5,'LineWidth', 2,"Color",[0,139/255,139/255]);
    hold on;
    plotvelocity52 = rad2deg(velocity52);
    plot(time52, plotvelocity52,'LineWidth', 2,"Color",[220/255,20/255,60/255]);
    hold on;
    plot(0,0);
    hold on;
    legend({'Axis 1', 'Axis 2','Axis 3'},'Location','southeast');
    hold on;

    figure(3);
    xlabel('Time(s)');
    ylabel('Acceleration(°/s^2)');
    
    plotacceleration5 = rad2deg(acceleration5);
    plot(time5, plotacceleration5,'LineWidth', 2,"Color",[0,139/255,139/255]);
    hold on; 
    plotacceleration52 = rad2deg(acceleration52);
    plot(time52, plotacceleration52,'LineWidth', 2,"Color",[220/255,20/255,60/255]);
    hold on;
    plot(0,0);
    hold on;
    legend({'Axis 1', 'Axis 2','Axis 3'},'Location','southeast');
    hold on;

    figure(4);
    xlabel('Time(s)');
    ylabel('Jerk(°/s^3)');
    
    plotjerk5 = rad2deg(jerk5);
    plot(time5, plotjerk5,'LineWidth', 2,"Color",[0,139/255,139/255]);
    hold on; 
    plotjerk52 = rad2deg(jerk52);
    plot(time52, plotjerk52,'LineWidth', 2,"Color",[220/255,20/255,60/255]);
    hold on;
    plot(0,0);
    hold on;
    legend({'Axis 1', 'Axis 2','Axis 3'},'Location','southeast');
    hold on;
    % 绘制竖直虚线
%     figure(1);
%     plot([t(1), t(1)], [-10, 90], '--k'); % '--k'表示黑色虚线
%     hold on;
%     plot([t(2), t(2)], [-10, 90], '--k'); % '--k'表示黑色虚线
%     hold on;
%     plot([t(3), t(3)], [-10, 90], '--k'); % '--k'表示黑色虚线
%     hold on;
%     plot([t(4), t(4)], [-10, 90], '--k'); % '--k'表示黑色虚线
%     hold on;
    legend({'Axis 1', 'Axis 2','Axis 3'},'Location','southeast'); % Add legend to differentiate the trajectories
end
