function jerk = quintic_jerk(t, coeffs_x, coeffs_y)
    time = linspace(0, t, 100);
    jerkx = zeros(1, 100);  % Create an array to store jerk values
    jerky = zeros(1, 100);
    for i = 1:100
        ti = time(i);
        jerkx(i) = abs(6*coeffs_x(5) + 24*coeffs_x(4)*ti + 60*coeffs_x(3)*ti^2);
        jerky(i) = abs(6*coeffs_y(5) + 24*coeffs_y(4)*ti + 60*coeffs_y(3)*ti^2);
    end
    jerk = sum(jerkx) + sum(jerky);  % Sum up jerk values for x and y
end

