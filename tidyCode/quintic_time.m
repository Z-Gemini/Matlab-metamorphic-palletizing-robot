function t_value = quintic_time(t0, t1, p0, p1, v0, v1, a0, a1)
    vmax = deg2rad(32);
    amax = deg2rad(100);
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
    
    % 解线性方程组（考虑奇异性）
    if rank(A) < size(A, 2)
        coefficients = pinv(A) * b;
    else
        coefficients = A \ b;
    end

    % 设置插值点数目（可根据需要进行调整）
    num_points = 100;

    % 生成时间向量
    time = linspace(t0, t1, num_points);

    % 计算对应时间点的位置、速度和加速度
    velocity = zeros(1, num_points);
    acceleration = zeros(1, num_points);

    for i = 1:num_points
        t = time(i);
        velocity(i) = coefficients(2) + 2*coefficients(3)*t + 3*coefficients(4)*t^2 + 4*coefficients(5)*t^3 + 5*coefficients(6)*t^4;
        acceleration(i) = 2*coefficients(3) + 6*coefficients(4)*t + 12*coefficients(5)*t^2 + 20*coefficients(6)*t^3;
    end
    
    if max(abs(velocity)) > vmax
        k1 = 1;
    else
        k1 = 0;
    end

    if max(abs(acceleration)) > amax
        k2 = 1;
    else
        k2 = 0;
    end

    t_value = time(100) + 1000*k1 + 1000*k2;
end