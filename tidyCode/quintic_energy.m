function energy = quintic_energy(T,coefficients1,coefficients2)
    if isequal(coefficients2, [0 0 0 0 0 0])
        energy = quintic_energy2(T,coefficients1);
    else
        energy = quintic_energy1(T,coefficients1,coefficients2);
    end
end

function energy = quintic_energy1(T,coefficients1,coefficients2)
    L = [0.6 0.35 0.2 0.608 0.4 0.455 0.196];

    num_steps = 100;  % Number of steps for interpolation
    energy = 0;  % Initialize total energy
    
    for step = 1:num_steps
        t = T * (step - 0.5) / num_steps;
        p_M1 = coefficients1(1) + coefficients1(2)*t + coefficients1(3)*t^2 + coefficients1(4)*t^3 + coefficients1(5)*t^4 + coefficients1(6)*t^5;
        v_M1 = coefficients1(2) + 2*coefficients1(3)*t + 3*coefficients1(4)*t^2 + 4*coefficients1(5)*t^3 + 5*coefficients1(6)*t^4;
        a_M1 = 2*coefficients1(3) + 6*coefficients1(4)*t + 12*coefficients1(5)*t^2 + 20*coefficients1(6)*t^3;
        
        p_M2 = coefficients2(1) + coefficients2(2)*t + coefficients2(3)*t^2 + coefficients2(4)*t^3 + coefficients2(5)*t^4 + coefficients2(6)*t^5;
        v_M2 = coefficients2(2) + 2*coefficients2(3)*t + 3*coefficients2(4)*t^2 + 4*coefficients2(5)*t^3 + 5*coefficients2(6)*t^4;
        a_M2 = 2*coefficients2(3) + 6*coefficients2(4)*t + 12*coefficients2(5)*t^2 + 20*coefficients2(6)*t^3;
    
        tt = (T/100 + T)/2;
        p_L1 = coefficients1(1) + coefficients1(2)*tt + coefficients1(3)*tt^2 + coefficients1(4)*tt^3 + coefficients1(5)*tt^4 + coefficients1(6)*tt^5;
        p_L2 = coefficients2(1) + coefficients2(2)*tt + coefficients2(3)*tt^2 + coefficients2(4)*tt^3 + coefficients2(5)*tt^4 + coefficients2(6)*tt^5;
    
        ttt = (2*tt+T/100)/2;
        p_R1 = coefficients1(1) + coefficients1(2)*ttt + coefficients1(3)*ttt^2 + coefficients1(4)*ttt^3 + coefficients1(5)*ttt^4 + coefficients1(6)*ttt^5;
        p_R2 = coefficients2(1) + coefficients2(2)*ttt + coefficients2(3)*ttt^2 + coefficients2(4)*ttt^3 + coefficients2(5)*ttt^4 + coefficients2(6)*ttt^5;
    
        % 计算角度值
        angles = [p_M1, p_M2, p_L1, p_L2, p_R1, p_R2];
        cos_values = cos(angles);
        sin_values = sin(angles);
        
        % 计算坐标
        xB = L(1) * cos_values(1);
        yB = L(1) * sin_values(1);
        xD = L(2) * cos_values(2);
        yD = L(2) * sin_values(2);
        rB = L(3);
        rD = L(4);
        [xC, yC] = intersection_of_circles(xD, yD, rD, xB, yB, rB, 0);
        
        xBt = L(1) * cos_values(3);
        yBt = L(1) * sin_values(3);
        xDt = L(2) * cos_values(4);
        yDt = L(2) * sin_values(4);
        [xCt, yCt] = intersection_of_circles(xDt, yDt, rD, xBt, yBt, rB, 0);
        
        xBtt = L(1) * cos_values(5);
        yBtt = L(1) * sin_values(5);
        xDtt = L(2) * cos_values(6);
        yDtt = L(2) * sin_values(6);
        [xCtt, yCtt] = intersection_of_circles(xDtt, yDtt, rD, xBtt, yBtt, rB, 0);
        
        % 三角板CDN转角
        phi3 = azimuth(xD,yD,xC,yC);
        phi3t = azimuth(xDt,yDt,xCt,yCt);
        phi3tt = azimuth(xDtt,yDtt,xCtt,yCtt);
        % 连杆BF转角
        phi4 = azimuth(xB,yB,xC,yC);
        phi4t = azimuth(xBt,yBt,xCt,yCt);
        phi4tt = azimuth(xBtt,yBtt,xCtt,yCtt);
        
        % 计算拉格朗日方程中间变量
        dphi3v_M1 = (phi3t - phi3)/(p_L1 - p_M1);%phi3对p_M1的偏导
        dphi3v_M2 = (phi3t - phi3)/(p_L2 - p_M2);%phi3对p_M2的偏导
        dphi4v_M1 = (phi4t - phi4)/(p_L1 - p_M1);%phi4对p_M1的偏导
        dphi4v_M2 = (phi4t - phi4)/(p_L2 - p_M2);%phi4对p_M2的偏导
        
        dphi3tv_M1t = (phi3tt - phi3t)/(p_R1 - p_L1);%phi3t对p_L1的偏导
        dphi3tv_M2t = (phi3tt - phi3t)/(p_R2 - p_L2);%phi3t对p_L2的偏导
        dphi4tv_M1t = (phi4tt - phi4t)/(p_R1 - p_L1);%phi4t对p_L1的偏导
        dphi4tv_M2t = (phi4tt - phi4t)/(p_R2 - p_L2);%phi4t对p_L2的偏导
        
        j11 = J11(p_M1,phi3,phi4,dphi3v_M1,dphi4v_M1);
        j11t = J11(p_L1,phi3t,phi4t,dphi3tv_M1t,dphi4tv_M1t);
        dJ11v_M1 = (j11t - j11)/(p_L1 - p_M1);
        dJ11v_M2 = (j11t - j11)/(p_L2 - p_M2);
        
        j22 = J22(p_M2,phi3,phi4,dphi3v_M2,dphi4v_M2);
        j22t = J22(p_L2,phi3t,phi4t,dphi3tv_M2t,dphi4tv_M2t);
        dJ22v_M1 = (j22t - j22)/(p_L1 - p_M1);
        dJ22v_M2 = (j22t - j22)/(p_L2 - p_M2);
        
        j12 = J12(p_M1,p_M2,phi3,phi4,dphi3v_M1,dphi4v_M1,dphi3v_M2,dphi4v_M2);
        j12t = J12(p_L1,p_L2,phi3t,phi4t,dphi3tv_M1t,dphi4tv_M1t,dphi3tv_M2t,dphi4tv_M2t);
        dJ12v_M1 = (j12t - j12)/(p_L1 - p_M1);
        dJ12v_M2 = (j12t - j12)/(p_L2 - p_M2);
        
        % 势能计算
        ep = Ep(p_M1,p_M2,phi3,phi4);
        ept = Ep(p_L1,p_L2,phi3t,phi4t);
        dEpv_M1 = (ept - ep)/(p_L1 - p_M1);
        dEpv_M2 = (ept - ep)/(p_L2 - p_M2);
        
        % 拉格朗日方程计算力矩
        tao1 = j11*a_M1 + j12*a_M2 +(1/2)*dJ11v_M1*(v_M1)^2 ...
               + dJ11v_M2*v_M1*v_M2 + (dJ12v_M2 - (1/2)*dJ22v_M1)*(v_M2)^2 + dEpv_M1;
        tao2 = j22*a_M2 + j12*a_M1 + (dJ12v_M1 - (1/2)*dJ11v_M2)*(v_M2)^2 ...
               + dJ22v_M1*v_M1*v_M2 + (1/2)*dJ22v_M2*(v_M2)^2 + dEpv_M2;
        % 计算功率
        energyi = abs(tao1 * v_M1) + abs(tao2 * v_M2);
        if isnan(energyi)
            energy = 0;
            continue;
        end
        energy = energy+energyi;
    end
end

function energy = quintic_energy2(T,coefficients)
    %t:总时间 a:1轴系数 b:2轴系数

    %杆长参数
    L = [0.6 0.35 0.2 0.608 0.4 0.455 0.196];
    %质量及转动惯量参数
    m1 = 9.328;%1轴
    m2 = 6.841;%2轴
    J1 = 42.995/100;%kg*m^2
    J2 = 14.538/100;
    m8 = 2;%负载
    m34 = 19.898;
    J34 = 148.865/100;
    %质心位置参数
    lSA = 0.288;
    lSE = 0.221;
    lSD1 = 0.490;
    ang_BDS = 12.072*pi/180;
    %重力加速度
    grv = 9.81;

    num_steps = 100;  % Number of steps for interpolation
    energy = 0;  % Initialize total energy
    
    for step = 1:num_steps
        t = T * (step - 0.5) / num_steps;
        %计算中间时刻及其附近两个时刻的两个轴的位移、速度、加速度
        theta1 = coefficients(1) + coefficients(2)*t + coefficients(3)*t^2 + coefficients(4)*t^3 + coefficients(5)*t^4 + coefficients(6)*t^5;
        dtheta1 = coefficients(2) + 2*coefficients(3)*t + 3*coefficients(4)*t^2 + 4*coefficients(5)*t^3 + 5*coefficients(6)*t^4;
        ddtheta1 = 2*coefficients(3) + 6*coefficients(4)*t + 12*coefficients(5)*t^2 + 20*coefficients(6)*t^3;
    
        tt = (T/100 + T)/2;
        theta1t = coefficients(1) + coefficients(2)*tt + coefficients(3)*tt^2 + coefficients(4)*tt^3 + coefficients(5)*tt^4 + coefficients(6)*tt^5;
    
        ttt = (2*tt+T/100)/2;
        theta1tt = coefficients(1) + coefficients(2)*ttt + coefficients(3)*ttt^2 + coefficients(4)*ttt^3 + coefficients(5)*ttt^4 + coefficients(6)*ttt^5;
        
        %计算中间时刻及其附近两个时刻的机构简图中C点坐标
        xB = L(1) * cos(theta1);
        yB = L(1) * sin(theta1);
        rB = L(6);
        xE = 0;
        yE = -L(7);
        rE = L(2);
        [xD,yD] = intersection_of_circles(xB,yB,rB,xE,yE,rE,1);
        
        xBt = L(1) * cos(theta1t);
        yBt = L(1) * sin(theta1t);
        [xDt,yDt] = intersection_of_circles(xBt,yBt,rB,xE,yE,rE,1);
        
        xBtt = L(1) * cos(theta1tt);
        yBtt = L(1) * sin(theta1tt);
        [xDtt,yDtt] = intersection_of_circles(xBtt,yBtt,rB,xE,yE,rE,1);
        
        % 三角板CDN合并连杆BF转角
        phi34 = azimuth(xD,yD,xB,yB);
        phi34t = azimuth(xDt,yDt,xBt,yBt);
        phi34tt = azimuth(xDtt,yDtt,xBtt,yBtt);
        
        % 2轴转角
        phi2 = azimuth(xE,yE,xD,yD);
        phi2t = azimuth(xE,yE,xDt,yDt);
        phi2tt = azimuth(xE,yE,xDtt,yDtt);
        
        % BF连杆
        phi4 = phi34;
        phi4t = phi34t;
        phi4tt = phi34tt;
        
        % 计算拉格朗日方程中间变量
        dxDdtheta1 = (xDt - xD)/(theta1t - theta1);%xD对theta1的偏导
        dxDtdtheta1t = (xDtt - xDt)/(theta1tt - theta1t);%xDt对theta1t的偏导
        dyDdtheta1 = (yDt - yD)/(theta1t - theta1);%yD对theta1的偏导
        dyDtdtheta1t = (yDtt - yDt)/(theta1tt - theta1t);%yDt对theta1t的偏导
        
        dphi34dtheta1 = (phi34t - phi34)/(theta1t - theta1);%phi34对theta1的偏导
        dphi34tdtheta1t = (phi34tt - phi34t)/(theta1tt - theta1t);%phi34t对theta1t的偏导
        dphi2dtheta1 = (phi2t - phi2)/(theta1t - theta1);%phi2对theta1的偏导
        dphi2tdtheta1t = (phi2tt - phi2t)/(theta1tt - theta1t);%phi2t对theta1t的偏导
        dphi4dtheta1 = (phi4t - phi4)/(theta1t - theta1);%phi4对theta1的偏导
        dphi4tdtheta1t = (phi4tt - phi4t)/(theta1tt - theta1t);%phi4t对theta1t的偏导
        
        dxS34dtheta1 = dxDdtheta1 + lSD1 * (-sin(phi34 - ang_BDS)) * dphi34dtheta1;
        dyS34dtheta1 = dyDdtheta1 + lSD1 * (cos(phi34 - ang_BDS)) * dphi34dtheta1;
        dxS34tdtheta1t = dxDtdtheta1t + lSD1 * (-sin(phi34t - ang_BDS)) * dphi34tdtheta1t;
        dyS34tdtheta1t = dyDtdtheta1t + lSD1 * (cos(phi34t - ang_BDS)) * dphi34tdtheta1t;
        dxS8dtheta1 = L(1) * (-sin(theta1)) + (L(3) +L(5)) * (-sin(phi4)) * dphi4dtheta1;
        dyS8dtheta1 = L(1) * cos(theta1) + (L(3) +L(5)) * cos(phi4) * dphi4dtheta1;
        dxS8tdtheta1t = L(1) * (-sin(theta1t)) + (L(3) +L(5)) * (-sin(phi4t)) * dphi4tdtheta1t;
        dyS8tdtheta1t = L(1) * cos(theta1t) + (L(3) +L(5)) * cos(phi4t) * dphi4tdtheta1t;
        
        Je34 = m34 * (dxS34dtheta1^2 + dyS34dtheta1^2) + J34 * dphi34dtheta1^2;
        Je34t = m34 * (dxS34tdtheta1t^2 + dyS34tdtheta1t^2) + J34 * dphi34tdtheta1t^2;
        Je2 = J2 * dphi2dtheta1^2;
        Je2t = J2 * dphi2tdtheta1t^2;
        Je1 = J1;
        %负载
        Je8 = m8 * (dxS8dtheta1^2 + dyS8dtheta1^2);
        Je8t = m8 * (dxS8tdtheta1t^2 + dyS8tdtheta1t^2);
        
        J = Je1 + Je2 +Je34 + Je8;
        Jt = Je1 + Je2t + Je34t +Je8t;
        
        dJdtheta1 = (Jt - J)/(theta1t - theta1);
        
        ep1 = m1 * grv * lSA * sin(theta1);
        ep2 = m2 * grv * (lSE * sin(phi2) - L(7));
        ep34 = m34 * grv *(L(2) * sin(phi2) + lSD1 * sin(phi34 - ang_BDS) - L(7));
        ep = ep1 + ep2 + ep34;
        
        ep1t = m1 * grv * lSA * sin(theta1t);
        ep2t = m2 * grv * (lSE * sin(phi2t) - L(7));
        ep34t = m34 * grv *(L(2) * sin(phi2t) + lSD1 * sin(phi34t - ang_BDS) - L(7));
        ept = ep1t + ep2t + ep34t;
        
        dEpdtheta1 = (ept - ep)/(theta1t - theta1);
        
        % 拉格朗日方程计算力矩
        tao = J * ddtheta1 + (1/2) * dJdtheta1 * dtheta1^2 + dEpdtheta1;
        energyi = abs(tao * dtheta1);
        if isnan(energyi)
            energyi = 0;
        end
        % 计算功率
        energy= energy+ energyi;
    end
end
%% 
function [x_target, y_target] = intersection_of_circles(x1, y1, r1, x2, y2, r2, flag)
    % 默认有两个解
    % 标志位 flag = 0默认输出靠右侧解 flag = 1输出靠左侧解
    dx = x2 - x1;
    dy = y2 - y1;
    k1 = dy / dx;
    k2 = -1 / k1;
    L_squared = dx^2 + dy^2;
    dm = (r2^2 - r1^2 + L_squared) / (2 * sqrt(L_squared));
    x0 = x1 + (dm / sqrt(L_squared)) * dx;
    y0 = y1 + (dm / sqrt(L_squared)) * dy;
    cm_squared = r1^2 - (x0 - x1)^2 - (y0 - y1)^2;
    if flag == 0
        x_target = x0 + (sqrt(cm_squared) / sqrt(1 + k2^2));
    else
        x_target = x0 - (sqrt(cm_squared) / sqrt(1 + k2^2));
    end
    y_target = y0 + k2 * (x_target - x0);
end

%% 二自由度运动 J11 不考虑辅助杆件
function J11_sum = J11(theta1, phi3, phi4, dphi3dtheta1, dphi4dtheta1)
    % 杆长参数
    L = [0.6 0.35 0.2 0.608 0.4 0.455 0.196];

    m3 = 11.604;  % 三角板
    m4 = 8.294;  % 连杆
    J1 = 42.995/100;  % kg*m^2
    J3 = 38.373/100;
    J4 = 31.093/100;
    m8 = 2;  % 负载
    lSB = 0.309;  % 连杆
    lSD = 0.319;  % 三角板
    ang_CDS = 5.437*pi/180;
    sin_phi3_CDS = sin(phi3 + ang_CDS);
    cos_phi3_CDS = cos(phi3 + ang_CDS);
    sin_phi4 = sin(phi4);
    cos_phi4 = cos(phi4);
    sin_theta1 = sin(theta1);
    cos_theta1 = cos(theta1);
    
    
    % 求解 J11_
    dxS3dtheta1 = lSD * (-sin_phi3_CDS) * dphi3dtheta1;
    dyS3dtheta1 = lSD * cos_phi3_CDS * dphi3dtheta1;
    J11_3 = m3 * (dxS3dtheta1^2 + dyS3dtheta1^2) + J3 * dphi3dtheta1^2;

    % 求解 J11_4
    dxS4dtheta1 = L(1)*(-sin_theta1) + lSB * (-sin_phi4) * dphi4dtheta1;
    dyS4dtheta1 = L(1)*cos_theta1 + lSB * cos_phi4 * dphi4dtheta1;
    J11_4 = m4 * (dxS4dtheta1^2 + dyS4dtheta1^2) + J4 * dphi4dtheta1^2;

    % J11_8 负载
    dxS8dtheta1 = L(1)*(-sin_theta1) + (L(3) + L(5)) * (-sin_phi4) * dphi4dtheta1;
    dyS8dtheta1 = L(1)*cos_theta1 + (L(3) + L(5)) * cos_phi4 * dphi4dtheta1;
    J11_8 = m8 * (dxS8dtheta1^2 + dyS8dtheta1^2);

    % 求和
    J11_sum = J1 + J11_3 + J11_4 + J11_8;
end

%% 二自由度运动 J12 不考虑辅助杆件
function J12_sum = J12(theta1,theta2,phi3,phi4,dphi3dtheta1,dphi4dtheta1,dphi3dtheta2,dphi4dtheta2)
    %杆长参数
    L = [0.6 0.35 0.2 0.608 0.4 0.455 0.196];

    %质量及转动惯量参数
    m3 = 11.604;%三角板
    m4 = 8.294;%连杆
    J3 = 38.373/100;
    J4 = 31.093/100;
    m8 = 2;%负载
    lSB = 0.309;%连杆
    lSD = 0.319;%三角板
    ang_CDS = 5.437*pi/180;

    J12_1 = 0;
    J12_2 = 0;
    %求解J11_3
    dxS3dtheta1 = lSD * (-sin(phi3 + ang_CDS)) * dphi3dtheta1;
    dxS3dtheta2 = L(2)*(-sin(theta2)) + lSD * (-sin(phi3 + ang_CDS)) * dphi3dtheta2;
    dyS3dtheta1 = lSD * cos(phi3 + ang_CDS) * dphi3dtheta1;
    dyS3dtheta2 = L(2)*cos(theta2) + lSD * cos(phi3 + ang_CDS) * dphi3dtheta2;
    J12_3 = m3 * (dxS3dtheta1 * dxS3dtheta2 + dyS3dtheta1 * dyS3dtheta2) + J3 * dphi3dtheta1 * dphi3dtheta2;
    %求解J11_4
    dxS4dtheta1 = L(1)*(-sin(theta1)) + lSB * (-sin(phi4)) * dphi4dtheta1;
    dxS4dtheta2 = lSB * (-sin(phi4)) * dphi4dtheta2;
    dyS4dtheta1 = L(1)*cos(theta1)+ lSB * cos(phi4) * dphi4dtheta1;
    dyS4dtheta2 = lSB * cos(phi4) * dphi4dtheta2;
    J12_4 = m4 * (dxS4dtheta1 * dxS4dtheta2 + dyS4dtheta1 * dyS4dtheta2) + J4 * dphi4dtheta1 * dphi4dtheta2;
    %求解J11_8 负载
    dxS8dtheta1 = L(1)*(-sin(theta1)) + (L(3) + L(5)) * (-sin(phi4)) * dphi4dtheta1;
    dxS8dtheta2 = (L(3) + L(5)) * (-sin(phi4)) * dphi4dtheta2;
    dyS8dtheta1 = L(1)*cos(theta1)+ (L(3) + L(5)) * cos(phi4) * dphi4dtheta1;
    dyS8dtheta2 = (L(3) + L(5)) * cos(phi4) * dphi4dtheta2;
    J12_8 = m8 * (dxS8dtheta1 * dxS8dtheta2 + dyS8dtheta1 * dyS8dtheta2);
    %求和
    J12_sum = J12_1 + J12_2 + J12_3 + J12_4 + J12_8;
end

%% 二自由度运动 J22 不考虑辅助杆件
function J22_sum = J22(theta2,phi3,phi4,dphi3dtheta2,dphi4dtheta2)
%杆长参数
L = [0.6 0.35 0.2 0.608 0.4 0.455 0.196];

%质量及转动惯量参数

m3 = 11.604;%三角板
m4 = 8.294;%连杆

J2 = 14.538/100;
J3 = 38.373/100;
J4 = 31.093/100;
m8 = 2;%负载
lSB = 0.309;%连杆
lSD = 0.319;%三角板
ang_CDS = 5.437*pi/180;

    J22_1 = 0;
    J22_2 = J2;
    %求解J22_3
    dxS3dtheta2 = L(2)*(-sin(theta2)) + lSD * (-sin(phi3 + ang_CDS)) * dphi3dtheta2;
    dyS3dtheta2 = L(2)*cos(theta2) + lSD * cos(phi3 + ang_CDS) * dphi3dtheta2;
    J22_3 = m3 * (dxS3dtheta2^2 + dyS3dtheta2^2) + J3 * dphi3dtheta2^2;
    %求解J22_4
    dxS4dtheta2 = lSB * (-sin(phi4)) * dphi4dtheta2;
    dyS4dtheta2 = lSB * cos(phi4) * dphi4dtheta2;
    J22_4 = m4 * (dxS4dtheta2^2 + dyS4dtheta2^2) + J4 * dphi4dtheta2^2;
    %求解J22_8 负载
    dxS8dtheta2 = (L(3) + L(5)) * (-sin(phi4)) * dphi4dtheta2;
    dyS8dtheta2 = (L(3) + L(5)) * cos(phi4) * dphi4dtheta2;
    J22_8 = m8 * (dxS8dtheta2^2 + dyS8dtheta2^2);
    %求和
    J22_sum = J22_1 + J22_2 + J22_3 + J22_4 + J22_8;
end

%% 势能计算
function ep = Ep(theta1,theta2,phi3,phi4)
%杆长参数
L = [0.6 0.35 0.2 0.608 0.4 0.455 0.196];
%质量及转动惯量参数
m1 = 9.328;%1轴
m2 = 6.841;%2轴
m3 = 11.604;%三角板
m4 = 8.294;%连杆

lSA = 0.288;
lSE = 0.221;
lSB = 0.309;%连杆
lSD = 0.319;%三角板
ang_CDS = 5.437*pi/180;

%重力加速度
grv = 9.81;

    ep1 = m1 * grv * lSA * sin(theta1);
    ep2 = m2 * grv * (lSE * sin(theta2) - L(7));
    
    ep3 = m3 * grv *(L(2) * sin(theta2) + lSD * sin(phi3 + ang_CDS) - L(7));
    ep4 = m4 * grv *(L(1) * sin(theta1) + lSB * sin(phi4)); 
    
    ep = ep1 + ep2 + ep3 + ep4;
end


