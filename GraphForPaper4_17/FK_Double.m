function [ Fx,Fy ] = FK_Double(theta1,theta2)
    %%  二自由度码垛机器人正运动学求解
    %  [ Fx,Fy ]为末端点坐标
    %  （theta1,theta2）为两个驱动关节的关节转角
    %% 机械常数
    l1 = 0.6;%AB杆长度
    l2 = 0.35;%DE杆长度
    l3 = 0.2;%BC杆长度
    l4 = 0.608;%CD杆长度
    l5 = 0.4;%CF杆长度
    l6 = 0.455;%DN杆长度
    l7 = 0.196;%AE杆长度
    
    Ax = 0;%A点横坐标
    Ay = 0;%A点纵坐标
    Ex = 0;%E点横坐标
    Ey = -l7;%E点纵坐标
    
    Bx = l1*cos(theta1);%B点横坐标
    By = l1*sin(theta1);%B点横坐标
    Dx = l2*cos(theta2)+Ex;%D点横坐标
    Dy = l2*sin(theta2)+Ey;%D点横坐标
    %% C点坐标求解
    lBD = sqrt((Bx-Dx).^2+(By-Dy).^2); 
    angle_DBC = acos((lBD.^2+l3.^2-l4.^2)/(2*lBD*l3));
    
    Angle_BD = azimuth([Bx,By;Dx,Dy]);
    Cx = Bx + cos(Angle_BD + angle_DBC) * l3;%C点横坐标
    Cy = By + sin(Angle_BD + angle_DBC) * l3;%C点纵坐标
    Angle_BC = azimuth([Bx,By;Cx,Cy]);%BC杆方位角
    Angle_CF = Angle_BC;%CF杆方位角
    
    %% 最终结果求解
    Fx = Cx + cos(Angle_CF) * l5;
    Fy = Cy + sin(Angle_CF) * l5;
end

