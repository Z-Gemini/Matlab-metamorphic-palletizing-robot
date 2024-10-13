function FW_Trajectory
%     t = [0.8909 1.3088 1.4069 1.7143 2.30];
    t = [0.235,0.435,1.556-01.4,1.756-0.14,2.380-0.14];
    thetaA1 = deg2rad(28.2);   % 抓取点  0.4922  
    thetaD1 = deg2rad(40.7);   % 放置点  0.7103
    g = deg2rad([28.04,29.49,44,44.03]); % 一自由度轨迹点 0.48939	0.5147	0.7679	0.76847
    
    thetaA2 = deg2rad(74.7); % 抓取点  0.4922  
    thetaD2 = deg2rad(81.9);   % 放置点  0.7103
    g2=[1.2893, 1.2982,1.3855,1.3857];
    
    %% 第一段五次
    t0 = 0; t02 = 0;
    theta0 = thetaA1;theta02 = thetaA2;
    theta1 = g(1);theta12 = g2(1);
    v0 = 0;v01 = 0;
    v1 = (g(2)-g(1))/(t(2)-t(1));v12 = (g2(2)-g2(1))/(t(2)-t(1));
    a0 = 0;a02 = 0;
    a1 = 0;a12 = 0;

    [time1, position1, velocity1, acceleration1,jerk1] = quintic_trajectory(t0,t(1),theta0,theta1,v0,v1,a0,a1);
    [time12, position12, velocity12, acceleration12,jerk12] = quintic_trajectory(t02,t(1),theta02,theta12,v01,v12,a02,a12);
    
    %% 第二段 一次
    vc1 = v1;
    [time2, position2, velocity2, acceleration2] = curveLine(theta1,t(1),t(2),vc1);
    [position22] = coordinate1(position2);
    
    %% 第三段：五次
    theta2 = g(2);
    theta3 = g(3);
    vc2 = (g(4)-g(3))/(t(4)-t(3));
    
    theta13 = 52/180*pi;   % 抓取点  0.4922  
    theta23 = 0/180*pi;   % 放置点  0.7103
    
    [time3, position3, velocity3, acceleration3,jerk3] = quintic_trajectory(t(2),t(3),theta2,theta3,vc1,vc2,a0,a1);
    [time33, position33, velocity33, acceleration33,jerk33] = quintic_trajectory(t(2),t(3),theta13,theta23,0,0,0,0);
    %% 第四段：一次
    theta3 = g(3);
    [time4, position4, velocity4, acceleration4] = curveLine(theta3,t(3),t(4),vc2);
    [position42] = coordinate1(position4);    

    %% 第五段：五次
    theta4 = g(4);theta42 = g2(4);
    theta5 = thetaD1;theta52 = thetaD2;
    [time5, position5, velocity5, acceleration5,jerk5] = quintic_trajectory(t(4),t(5),theta4,theta5,vc2,0,a0,a1);
    [time52, position52, velocity52, acceleration52,jerk52] = quintic_trajectory(t(4),t(5),theta42,theta52,0,0,0,0);

    xt1 = zeros(1,100);yt1 = zeros(1,100);
    xc1 = zeros(1,100);yc1 = zeros(1,100);
    xt2 = zeros(1,100);yt2 = zeros(1,100);
    xc2 = zeros(1,100);yc2 = zeros(1,100);
    xt3 = zeros(1,100);yt3 = zeros(1,100);
    
    for i = 1:100
        [xt1(i),yt1(i)] = ForwardKimemicDouble(position1(i),position12(i));%二自由度正运动学求解
        [xc1(i),yc1(i)] = ForwardKimemicSingle(position2(i));%单自由度正运动学求解
        [xt2(i),yt2(i)] = ForwardKimemicSingle(position3(i));%单自由度正运动学求解
        [xc2(i),yc2(i)] = ForwardKimemicSingle(position4(i));%单自由度正运动学求解
        [xt3(i),yt3(i)] = ForwardKimemicDouble(position5(i),position52(i));%二自由度正运动学求解
    end
    figure(5)
    plot(xt1+0.12,yt1-0.09,'color',[0,139/255,139/255] ,'linewidth',3);%二自由度末端轨迹
    hold on;
    plot(xc1+0.12,yc1-0.09,'color',[220/255,20/255,60/255],'linewidth',3);%单自由度末端轨迹
    hold on;
    plot(xt2+0.12,yt2-0.09,'color',[218/255,165/255,32/255],'linewidth',3);%单自由度末端轨迹
    hold on;
    plot(xc2+0.12,yc2-0.09,'color',[220/255,20/255,60/255],'linewidth',3);%单自由度末端轨迹
    hold on;
    plot(xt3+0.12,yt3-0.09,'color',[0,139/255,139/255],'linewidth',3);%二自由度末端轨迹
    xlabel('x(m)');
    ylabel('y(m)');
    legend('2-DOF','Metamorphic','1-DOF');
end
    
%% fw_double
function [ Fx,Fy ] = ForwardKimemicDouble( theta1,theta2 )
    %  [ Fx,Fy ]为末端点坐标
    %  （theta1,theta2）为两个驱动关节的关节转角
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
%% ForwardKimemicSingle
function [ Fx,Fy ] = ForwardKimemicSingle( theta1 )
    %  [ Fx,Fy ]为末端点坐标
    %  theta1为驱动关节的关节转角
    
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


