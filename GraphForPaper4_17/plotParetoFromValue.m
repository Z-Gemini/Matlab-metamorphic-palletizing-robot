load("D:\desktop\mydata\桌面\code\matlab\tidyCode\result\文献采用的结果.mat");
    % 确定帕累托最优解
    paretoIdx = paretofront_custom(fval);
    
    % 绘制帕累托前沿
    figure;
    scatter3(fval(:,1), fval(:,2), fval(:,3), 'filled');
    hold on;
    scatter3(fval(paretoIdx,1), fval(paretoIdx,2), fval(paretoIdx,3), 'r', 'filled');
    hold off;
    xlabel('Time(s)');
    ylabel('Energy(J)');
    zlabel('Jerk(rad/s^3)');
    title('Pareto Front');
    legend('Non-dominant solutions', 'Pareto front');
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
