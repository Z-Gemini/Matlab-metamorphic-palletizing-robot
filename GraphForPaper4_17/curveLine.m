%% 匀速图像
function [time, position, velocity, acceleration] = curveLine(theta0,t0,t1,v)
    % 初始化时间向量
    num_points = 100;
    time = linspace(t0, t1, num_points);

    % 计算位置向量
    position = theta0 + v * (time - t0);

    % 计算速度向量
    velocity = v * ones(size(time));

    % 计算加速度向量
    acceleration = zeros(size(time));
    % Create a string representing the position formula
    position_formula = sprintf('Position = %g + %g * (time - %g)', theta0, v, t0);

    % Display the position formula
    disp('Position Formula:');
    disp(position_formula);

end
