% Define parameters
lambda_values = [10.6, 9.7, 9.2, 8.0, 7.8, 8.1, 7.8, 8.1, 9.1, 9.9, 10.6, 10.6];
k_values = [2.0, 2.0, 2.0, 1.9, 1.9, 1.9, 1.9, 1.9, 2.0, 1.9, 2.0, 2.0];
num_samples = 100000;  % Number of Monte Carlo samples
cut_in_speed = 3.5;    % Minimum operational wind speed
cut_off_speed = 25.0;  % Maximum operational wind speed
optimal_speed_limit = 14.0; % Speed up to which correlation is strong

% Initialize results storage
monthly_results = struct();

for month = 1:12
    % Get Weibull parameters for the current month
    lambda = lambda_values(month);
    k = k_values(month);

    % Generate random wind speeds using the Weibull distribution
    wind_speeds = wblrnd(lambda, k, num_samples, 1);

    % Apply wind speed truncation (within turbine's operational range)
    truncated_wind_speeds = wind_speeds(wind_speeds >= cut_in_speed & wind_speeds <= cut_off_speed);

    % Further limit the data to the range where correlation is strong (up to 14 m/s)
    correlated_wind_speeds = truncated_wind_speeds(truncated_wind_speeds <= optimal_speed_limit);
    
    % Calculate power output using the truncated wind speeds
    correlated_power_outputs = P(correlated_wind_speeds);

    % Calculate optimal beta (control variate coefficient)
    var_wind_speed = var(correlated_wind_speeds);
    cov_output_wind = cov(correlated_power_outputs, correlated_wind_speeds);
    beta = -cov_output_wind(1, 2) / var_wind_speed;

    % Define the corrected estimator
    mean_wind_speed = mean(correlated_wind_speeds);
    corrected_power_outputs = correlated_power_outputs - beta * (correlated_wind_speeds - mean_wind_speed);

    % Calculate the corrected confidence interval
    mean_corrected = mean(corrected_power_outputs);
    std_error_corrected = std(corrected_power_outputs) / sqrt(length(correlated_wind_speeds));
    z_value = norminv(0.975);  % 95% confidence interval
    ci_corrected = [mean_corrected - z_value * std_error_corrected, ...
                    mean_corrected + z_value * std_error_corrected];

    % Store results for the month
    monthly_results(month).month = month;
    monthly_results(month).mean_corrected = mean_corrected;
    monthly_results(month).ci_corrected = ci_corrected;

    % Display results
    fprintf('Month %d:\n', month);
    fprintf('  Corrected CI: [%.2f, %.2f] Watts\n\n', ci_corrected(1), ci_corrected(2));
end
