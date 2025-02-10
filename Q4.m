% Antithetic Sampling for Power Curve Estimation (All Months)
% Define parameters
num_samples = 100000;            % Number of samples
cut_in_speed = 3.5;              % Minimum operational wind speed
cut_off_speed = 25;              % Maximum operational wind speed

% Load power curve function (assumed to be in P(v))
load powercurve_V164.mat         % Adjust as needed based on actual data

% Weibull distribution parameters for each month
lambda_values = [10.6, 9.7, 9.2, 8.0, 7.8, 8.1, 7.8, 8.1, 9.1, 9.9, 10.6, 10.6];
k_values = [2.0, 2.0, 2.0, 1.9, 1.9, 1.9, 1.9, 1.9, 2.0, 1.9, 2.0, 2.0];

% Initialize results storage
monthly_results = struct();

% Loop through each month
for month = 1:12
    % Weibull parameters for the current month
    lambda_w = lambda_values(month);
    k_w = k_values(month);

    % Generate samples from Weibull distribution
    U = wblrnd(lambda_w, k_w, num_samples, 1);
    U = U(U >= cut_in_speed & U <= cut_off_speed);  % Truncate within operational range

    % Generate antithetic samples using CDF inversion
    % Compute the CDF values for the original samples
    F_U = wblcdf(U, lambda_w, k_w);
    % Generate antithetic samples by reflecting around 0.5 in the CDF space
    F_U_tilde = 1 - F_U;
    U_tilde = wblinv(F_U_tilde, lambda_w, k_w);
    U_tilde = U_tilde(U_tilde >= cut_in_speed & U_tilde <= cut_off_speed);  % Ensure within operational range

    % Ensure both arrays have the same size by trimming the larger array
    min_length = min(length(U), length(U_tilde));
    U = U(1:min_length);
    U_tilde = U_tilde(1:min_length);

    % Evaluate power curve outputs
    power_U = P(U);
    power_U_tilde = P(U_tilde);

    % Construct antithetic estimator
    W = (power_U + power_U_tilde) / 2;

    % Calculate mean power estimate
    mean_power_estimate = mean(W);

    % Calculate variance reduction
    var_original = var(power_U);
    var_antithetic = var(W);
    reduction_percentage = (1 - (var_antithetic / var_original)) * 100;

    % Compute 95% confidence interval
    std_error = std(W) / sqrt(length(W));
    z_value = norminv(0.975);
    ci_lower = mean_power_estimate - z_value * std_error;
    ci_upper = mean_power_estimate + z_value * std_error;

    % Store results for the month
    monthly_results(month).month = month;
    monthly_results(month).mean_power = mean_power_estimate;
    monthly_results(month).ci = [ci_lower, ci_upper];
    monthly_results(month).var_reduction = reduction_percentage;

    % Display results
    fprintf('Month %d Results:\n', month);
    fprintf('  Mean Power Estimate: %.2f W\n', mean_power_estimate);
    fprintf('  95%% Confidence Interval: [%.2f, %.2f] W\n', ci_lower, ci_upper);
    fprintf('  Variance Reduction: %.2f%%\n\n', reduction_percentage);
end

% Display final summary for all months
disp('Final Summary for All Months:');
for month = 1:12
    fprintf('Month %d: Mean Power = %.2f W, CI = [%.2f, %.2f] W, Variance Reduction = %.2f%%\n', ...
            month, monthly_results(month).mean_power, ...
            monthly_results(month).ci(1), monthly_results(month).ci(2), ...
            monthly_results(month).var_reduction);
end