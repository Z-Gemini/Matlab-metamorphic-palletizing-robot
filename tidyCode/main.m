%% main
function main()
    waypoints = deg2rad([28.2, 74.7; 0,0; 0,0; 0,0; 0,0; 40.7, 81.9]);
    
    funNSGAII = @(x) objective_function(x, waypoints);
    A = [1,-1,0,0;0,1,-1,0;0,0,1,-1;0,0,0,0];
    b = [-0.0105;0;-0.0105;0];
    
    lb = [0 0 0 0];
    ub = deg2rad([120 120 120 120]);
    optionsNSGA = optimoptions('gamultiobj', 'PopulationSize', 500, 'MaxGenerations', 80);
    [x, fval] = gamultiobj(funNSGAII, 4, A, b, [], [], lb, ub, optionsNSGA);
    % 确定帕累托最优解
    paretoIdx = paretofront_custom(fval);
    
    % 绘制帕累托前沿
    figure;
    scatter3(fval(:,1), fval(:,2), fval(:,3), 'filled');
    hold on;
    scatter3(fval(paretoIdx,1), fval(paretoIdx,2), fval(paretoIdx,3), 'r', 'filled');
    hold off;
    xlabel('Time');
    ylabel('Energy');
    zlabel('Jerk');
    title('Pareto Front');
    legend('Non-dominant solutions', 'Pareto front');
    for i=1:length(x)
        disp(x(i),fval(i))
    end
end
%%
function paretoIdx = paretofront_custom(fval)
    numSolutions = size(fval, 1);

    dominates = false(numSolutions, numSolutions);
    dominatedByCount = zeros(numSolutions, 1);

    for i = 1:numSolutions
        for j = 1:numSolutions
            if all(fval(j, :) <= fval(i, :)) && any(fval(j, :) < fval(i, :))
                dominates(i, j) = true;
                dominatedByCount(j) = dominatedByCount(j) + 1;
            end
        end
    end
    paretoIdx = dominatedByCount == 0;
end

%% multi-objective_function
function f = objective_function(x, waypoints)
    if x(2) == x(3)
        f = [NaN ,NaN ,NaN];
        return;
    end
    y = coordinate(x);
    waypoints(2:5, :) = [x(1), y(1); x(2), y(2); x(3), y(3); x(4), y(4)];
    
    % 定义多目标优化的目标函数
    [t,coeffs_x,coeffs_y] = coefficients(waypoints);

    time_cost = sum(t);
    if time_cost > 1000
        time_cost = NaN;
    end
    energy_cost = calculate_energy_cost(t,coeffs_x,coeffs_y);
    if energy_cost > 10000
        time_cost = NaN;
    end
    impact_cost = calculate_impact_cost(t,coeffs_x,coeffs_y);
    % 将多个目标组合成一个目标向量
    f = [time_cost, energy_cost, impact_cost];
end
%% coefficients
function [t,coeffs_x,coeffs_y] = coefficients(waypoints)
    coeffs_x = zeros(5, 6);
    coeffs_y = zeros(5, 6); 
    t = zeros(1, 5);
    % 第二段：一次轨迹系数
    [t(2), v2] = compute_time_velocity(waypoints(2,1),waypoints(3,1));
    coeffs_x(2,:) = [waypoints(2,1),v2,0,0,0,0];

    % 第四段：一次轨迹系数
    [ t(4), v4] = compute_time_velocity(waypoints(4,1), waypoints(5,1));
    coeffs_x(4,:) = [waypoints(4,1),v4,0,0,0,0];

    % 第一段：五次轨迹系数
    t1_x = PSO_time(0, waypoints(1,1), waypoints(2,1), 0, v2, 0, 0);
    t1_y = PSO_time(0, waypoints(1,2), waypoints(2,2), 0, 0, 0, 0);
    t(1) = max(t1_x,t1_y);
    coeffs_x(1,:) = quintic_coeff(0, t(1), waypoints(1,1), waypoints(2,1), 0, v2, 0, 0);
    coeffs_y(1,:) = quintic_coeff(0, t(1), waypoints(1,2), waypoints(2,2), 0, 0, 0, 0);

    % 第三段：五次轨迹
    t(3) = PSO_time(0, waypoints(3,1), waypoints(4,1), v2, v4, 0, 0);
    coeffs_x(3,:) = quintic_coeff(0, t(3), waypoints(3,1), waypoints(4,1), v2, v4, 0, 0);

    % 第五段：五次轨迹
    t5_x = PSO_time(0, waypoints(5,1), waypoints(6,1), v4, 0, 0, 0);
    t5_y = PSO_time(0, waypoints(5,2), waypoints(6,2), 0, 0, 0, 0);
    t(5) = max(t5_x,t5_y);
    coeffs_x(5,:) = quintic_coeff(0, t(5), waypoints(5,1), waypoints(6,1), v4, 0, 0, 0);
    coeffs_y(5,:) = quintic_coeff(0, t(5), waypoints(5,2), waypoints(6,2), 0, 0, 0, 0);
end
%% PSO_time
function time = PSO_time(t0, p0, p1, v0, v1, a0, a1)
    % Call particleswarm
    narvs = 1; % 变量个数
    t_lb = 0.01;
    t_ub = 10;
    options = optimoptions('particleswarm', 'SwarmSize', 50, 'MaxIterations', 100,'Display', 'off');
    [~, time] = particleswarm(@(t) quintic_time(t0, t, p0, p1, v0, v1, a0, a1), narvs, t_lb, t_ub, options);
end
%% calculate_energy_cost
function energy_cost = calculate_energy_cost(t,coeffs_x,coeffs_y)
    % 1轴能量消耗
    energy1= quintic_energy(t(1),coeffs_x(1,:),coeffs_y(1,:));
    energy3 = quintic_energy(t(3),coeffs_x(3,:),coeffs_y(3,:));
    energy5 = quintic_energy(t(5),coeffs_x(5,:),coeffs_y(5,:));
    energy_cost = energy1+energy3+energy5;
end
%% calculate_impact_cost
function impact_cost = calculate_impact_cost(t,coeffs_x,coeffs_y)
    % 计算冲击成本，可以根据需要调整具体的计算方法
    impact1 = quintic_jerk(t(1),coeffs_x(1,:),coeffs_y(1,:));
    impact3 = quintic_jerk(t(1),coeffs_x(1,:),coeffs_y(1,:));
    impact5 = quintic_jerk(t(1),coeffs_x(1,:),coeffs_y(1,:));
    impact_cost = impact1+impact3+impact5;
    if impact_cost > 10000
        impact_cost = NaN;
    end
end
