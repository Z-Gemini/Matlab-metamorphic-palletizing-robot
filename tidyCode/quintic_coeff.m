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
