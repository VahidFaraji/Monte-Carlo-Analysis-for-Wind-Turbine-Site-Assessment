% Load power curve data
load powercurve_V164.mat

% Define parameters for Weibull distribution (monthly values)
lambda_values = [10.6, 9.7, 9.2, 8.0, 7.8, 8.1, 7.8, 8.1, 9.1, 9.9, 10.6, 10.6];
k_values = [2.0, 2.0, 2.0, 1.9, 1.9, 1.9, 1.9, 1.9, 2.0, 1.9, 2.0, 2.0];

% Parameters for Gamma distribution (Instrumental Distribution)
alpha_gamma = 2;      % Shape parameter of Gamma
beta_gamma = 0.23;    % Rate parameter of Gamma

% Number of samples and wind speed limits
num_samples = 500000;
cut_in_speed = 3.5;   % Minimum operational wind speed
cut_off_speed = 25;   % Maximum operational wind speed

% Initialize results storage
monthly_results = struct();

for month = 1:12
    % Get Weibull parameters for the current month
    lambda_w = lambda_values(month);
    k_w = k_values(month);

    % Generate samples from the Gamma distribution
    samples_gamma = gamrnd(alpha_gamma, 1/beta_gamma, num_samples, 1);

    % Apply wind speed truncation to match turbine operational range
    truncated_samples = samples_gamma(samples_gamma >= cut_in_speed & samples_gamma <= cut_off_speed);

    % Compute importance weights
    weights = (k_w / lambda_w) .* (truncated_samples / lambda_w).^(k_w - 1) .* ...
              exp(-(truncated_samples / lambda_w).^k_w) ./ ...
              ((beta_gamma^alpha_gamma / gamma(alpha_gamma)) .* truncated_samples.^(alpha_gamma - 1) .* exp(-beta_gamma * truncated_samples));

    % Calculate power output for truncated samples
    power_outputs = P(truncated_samples);

    % Compute weighted mean and variance
    weighted_mean_power = sum(weights .* power_outputs) / sum(weights);
    weighted_variance_power = sum(weights .* (power_outputs - weighted_mean_power).^2) / sum(weights);

    % Calculate effective sample size for proper scaling of the standard error
    n_eff = (sum(weights)^2) / sum(weights.^2);

    % Calculate the 95% confidence interval
    std_error_weighted = sqrt(weighted_variance_power / n_eff);
    z_value = norminv(0.975);  % 95% confidence interval
    ci_weighted = [weighted_mean_power - z_value * std_error_weighted, ...
                   weighted_mean_power + z_value * std_error_weighted];

    % Store results for the month
    monthly_results(month).month = month;
    monthly_results(month).mean_weighted = weighted_mean_power;
    monthly_results(month).ci_weighted = ci_weighted;

    % Display results
    fprintf('N_eff %d:\n', n_eff);
    fprintf('Month %d:\n', month);
    fprintf('  Weighted Mean Power: %.2f Watts\n', weighted_mean_power);
    fprintf('  Weighted CI: [%.2f, %.2f] Watts\n\n', ci_weighted(1), ci_weighted(2));
end

% Display final results for all months
disp('Final Results for Each Month:');
for month = 1:12
    fprintf('Month %d: Mean Power = %.2f W, 95%% CI = [%.2f, %.2f] W\n', ...
            month, monthly_results(month).mean_weighted, ...
            monthly_results(month).ci_weighted(1), monthly_results(month).ci_weighted(2));
end
