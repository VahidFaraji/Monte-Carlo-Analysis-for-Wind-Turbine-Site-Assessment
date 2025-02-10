% Constants and parameters
rho = 1.225;             % Air density at sea level (kg/m^3)
d = 164;                 % Rotor diameter (m)
cut_in_speed = 3.5;      % Minimum operational wind speed (m/s)
cut_off_speed = 25;      % Maximum operational wind speed (m/s)

% Monthly Weibull parameters (from the table)
lambda_values = [10.6, 9.7, 9.2, 8.0, 7.8, 8.1, 7.8, 8.1, 9.1, 9.9, 10.6, 10.6];
k_values = [2.0, 2.0, 2.0, 1.9, 1.9, 1.9, 1.9, 1.9, 1.9, 1.9, 2.0, 2.0];

% Initialize results storage
mean_power_coefficients = zeros(1, 12);
ci_lower_bounds = zeros(1, 12);
ci_upper_bounds = zeros(1, 12);
effective_sample_sizes = zeros(1, 12);

% Load the power curve function P(v)
load powercurve_V164.mat  % Assumes P(v) is loaded

% Loop through each month
for month = 1:12
    lambda_w = lambda_values(month);
    k_w = k_values(month);

    % Generate wind speed samples for the month
    num_samples = 500000;
    wind_speeds = wblrnd(lambda_w, k_w, num_samples, 1);

    % Apply truncation to wind speeds within turbine operational limits
    wind_speeds = wind_speeds(wind_speeds >= cut_in_speed & wind_speeds <= cut_off_speed);

    % Calculate total wind power using formula (1)
    P_tot = (1/2) * rho * (pi * d^2 / 4) * wind_speeds.^3;

    % Calculate actual power output from the power curve
    actual_power = P(wind_speeds);

    % Calculate power coefficient samples
    power_coefficient_samples = actual_power ./ P_tot;
    mean_power_coefficients(month) = mean(power_coefficient_samples);

    % Calculate 95% confidence interval
    n_eff = length(wind_speeds);  % Effective sample size after truncation
    std_error = std(power_coefficient_samples) / sqrt(n_eff);
    z_value = norminv(0.975);  % 95% confidence interval z-score
    ci_lower_bounds(month) = mean_power_coefficients(month) - z_value * std_error;
    ci_upper_bounds(month) = mean_power_coefficients(month) + z_value * std_error;
    effective_sample_sizes(month) = n_eff;
end

% Display results
fprintf('\nMonthly Power Coefficients and Confidence Intervals:\n');
for month = 1:12
    fprintf('Month %2d: Mean Power Coefficient = %.4f, 95%% CI = [%.4f, %.4f], Effective Samples = %d\n', ...
            month, mean_power_coefficients(month), ci_lower_bounds(month), ci_upper_bounds(month), effective_sample_sizes(month));
end
