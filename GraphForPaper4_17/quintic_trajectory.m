function [time, position, velocity, acceleration,jerk,coefficients] = quintic_trajectory(t0, t1, p0, p1, v0, v1, a0, a1)
    % 增加插入点来消除荣格现象的五次多项式轨迹规划

    % 设置时间段
    dt = t1 - t0;
    % 计算多项式系数矩阵
    A = [1, t0, t0^2, t0^3, t0^4, t0^5;
         0, 1, 2*t0, 3*t0^2, 4*t0^3, 5*t0^4;
         0, 0, 2, 6*t0, 12*t0^2, 20*t0^3;
         1, t1, t1^2, t1^3, t1^4, t1^5;
         0, 1, 2*t1, 3*t1^2, 4*t1^3, 5*t1^4;
         0, 0, 2, 6*t1, 12*t1^2, 20*t1^3];

    b = [p0; v0; a0; p1; v1; a1];

    % 计算多项式系数
    coefficients = A \ b;

    % 设置插值点数目（可根据需要进行调整）
    num_points = 100;

    % 生成时间向量
    time = linspace(t0, t1, num_points);

    % 计算对应时间点的位置、速度和加速度
    position = zeros(1, num_points);
    velocity = zeros(1, num_points);
    acceleration = zeros(1, num_points);
    jerk = zeros(1,num_points);

    for i = 1:num_points
        t = time(i);
        position(i) = coefficients(1) + coefficients(2)*t + coefficients(3)*t^2 + coefficients(4)*t^3 + coefficients(5)*t^4 + coefficients(6)*t^5;
        velocity(i) = coefficients(2) + 2*coefficients(3)*t + 3*coefficients(4)*t^2 + 4*coefficients(5)*t^3 + 5*coefficients(6)*t^4;
        acceleration(i) = 2*coefficients(3) + 6*coefficients(4)*t + 12*coefficients(5)*t^2 + 20*coefficients(6)*t^3;
        jerk(i) = 6*coefficients(5) + 24*coefficients(4)*t + 60*coefficients(3)*t^2;
    end
    disp("表达式：")
    disp("p(t) = " + string(coefficients(1)) + " + " + ...
    string(coefficients(2)) + " * t + " + ...
    string(coefficients(3)) + " * t^2 + " + ...
    string(coefficients(4)) + " * t^3 + " + ...
    string(coefficients(5)) + " * t^4 + " + ...
    string(coefficients(6)) + " * t^5");
end
