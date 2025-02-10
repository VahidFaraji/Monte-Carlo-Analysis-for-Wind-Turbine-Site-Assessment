load powercurve_V164.mat
lambda_values = [10.6, 9.7, 9.2, 8.0, 7.8, 8.1, 7.8, 8.1, 9.1, 9.9, 10.6, 10.6];
k_values = [2.0, 2.0, 2.0, 1.9, 1.9, 1.9, 1.9, 1.9, 2.0, 1.9, 2.0, 2.0];

num_samples = 100000;
cut_in_speed = 3.5;
cut_off_speed = 25;

monthly_results = struct();

for month = 1:12
    lambda = lambda_values(month);
    k = k_values(month);

    wind_speeds = wblrnd(lambda, k, num_samples, 1);
    power_outputs_standard = P(wind_speeds);
    
    mean_power_standard = mean(power_outputs_standard);
    std_error_standard = std(power_outputs_standard) / sqrt(num_samples);
    z_value = norminv(0.975);  % 95% confidence interval
    ci_standard = [mean_power_standard - z_value * std_error_standard, ...
                   mean_power_standard + z_value * std_error_standard];

    % Monte Carlo simulation for the truncated version
    truncated_wind_speeds = wind_speeds(wind_speeds >= cut_in_speed & wind_speeds <= cut_off_speed);
    truncated_power_outputs = P(truncated_wind_speeds);

    % Calculate confidence interval for truncated version
    mean_power_truncated = mean(truncated_power_outputs);
    std_error_truncated = std(truncated_power_outputs) / sqrt(length(truncated_wind_speeds));
    ci_truncated = [mean_power_truncated - z_value * std_error_truncated, ...
                    mean_power_truncated + z_value * std_error_truncated];

    % Store results for the month
    monthly_results(month).month = month;
    monthly_results(month).mean_power_standard = mean_power_standard;
    monthly_results(month).ci_standard = ci_standard;
    monthly_results(month).mean_power_truncated = mean_power_truncated;
    monthly_results(month).ci_truncated = ci_truncated;

    fprintf('Month %d:\n', month);
    fprintf('  Standard CI: [%.2f, %.2f]\n', ci_standard(1), ci_standard(2));
    fprintf('  Truncated CI: [%.2f, %.2f]\n\n', ci_truncated(1), ci_truncated(2));
end

