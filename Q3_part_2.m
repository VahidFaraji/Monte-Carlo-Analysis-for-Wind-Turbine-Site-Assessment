%% Question 3 part b 

% Define parameters
k = 1.96;               % Shape parameter for Weibull distribution
lambda = 9.13;          % Scale parameter for Weibull distribution
alpha = 0.638;          % Copula parameter
p = 3;                  % Copula parameter
q = 1.5;                % Copula parameter
num_samples = 100000;   % Number of samples for Monte Carlo simulation
v_min = 3.5;                  % Minimum wind speed (operational area)
v_max = 25;                   % Maximum wind speed (operational area)


load('powercurve_V164.mat');  
P = @(v) P(v);               

% Generate independent Weibull samples
v1_samples = wblrnd(lambda, k, num_samples, 1);
v2_samples = wblrnd(lambda, k, num_samples, 1);

valid_indices = (v1_samples >= v_min & v1_samples <= v_max) & ...
                (v2_samples >= v_min & v2_samples <= v_max);

% Compute CDF values for the samples using the built-in function
F_v1 = wblcdf(v1_samples, lambda, k);
F_v2 = wblcdf(v2_samples, lambda, k);

% Apply the copula correction to obtain joint density
joint_density = wblpdf(v1_samples, lambda, k) .* wblpdf(v2_samples, lambda, k) .* ...
    (1 + alpha * (1 - F_v1.^p).^(q - 1) .* (1 - F_v2.^p).^(q - 1) .* ...
    (F_v1.^p .* (1 + p * q) - 1) .* (F_v2.^p .* (1 + p * q) - 1));

% Define the instrumental density (independent Weibull assumption)
instrumental_density = wblpdf(v1_samples, lambda, k) .* wblpdf(v2_samples, lambda, k);

% Compute importance sampling weights
weights = joint_density ./ instrumental_density;
mean_weight = mean(weights);
max_weight = max(weights);
min_weight = min(weights);

power_v1 = P(v1_samples)/1e6;
power_v2 = P(v2_samples)/1e6;

mean_power_v1 = sum(weights .* power_v1) / sum(weights);
mean_power_v2 = sum(weights .* power_v2) / sum(weights);

% Estimate the covariance between power outputs
covariance_estimate = sum(weights .* (power_v1 - mean_power_v1) .* (power_v2 - mean_power_v2)) / sum(weights);

%fprintf('Sample power outputs: %.2f kW, %.2f kW, %.2f MW\n', power_v1(1), power_v1(2), power_v1(3));
fprintf('Estimated covariance between power outputs: %.4f MW^2\n', covariance_estimate);


% Number of bootstrap samples
num_bootstrap_samples = 1000;
bootstrap_covariances = zeros(num_bootstrap_samples, 1);

% Generate bootstrap samples and compute covariance for each
for i = 1:num_bootstrap_samples
    bootstrap_indices = randi(num_samples, num_samples, 1);
    bootstrap_weights = weights(bootstrap_indices);
    bootstrap_power_v1 = power_v1(bootstrap_indices);
    bootstrap_power_v2 = power_v2(bootstrap_indices);

    % Recalculate covariance for the bootstrap sample
    mean_bootstrap_v1 = sum(bootstrap_weights .* bootstrap_power_v1) / sum(bootstrap_weights);
    mean_bootstrap_v2 = sum(bootstrap_weights .* bootstrap_power_v2) / sum(bootstrap_weights);
    bootstrap_covariances(i) = sum(bootstrap_weights .* (bootstrap_power_v1 - mean_bootstrap_v1) .* (bootstrap_power_v2 - mean_bootstrap_v2)) / sum(bootstrap_weights);
end

% Calculate the confidence interval (percentile method)
ci_lower = prctile(bootstrap_covariances, 2.5);
ci_upper = prctile(bootstrap_covariances, 97.5);

fprintf('95%% Confidence Interval for covariance (Bootstrap): [%.4f, %.4f] MW^2\n', ci_lower, ci_upper);

mean_P1P2 = mean(weights .* power_v1 .* power_v2);  % E[P(V1)P(V2)]
variance_P1P2 = var(weights .* power_v1 .* power_v2);  % Variance of weighted samples

std_P1P2 = sqrt(variance_P1P2 / num_samples);
conf_width = 1.96 * std_P1P2;    
ci_lower = covariance_estimate - conf_width;
ci_upper = covariance_estimate + conf_width;

fprintf('95%% Confidence Interval for covariance (Analytical): [%.4f, %.4f] MW^2\n', ci_lower, ci_upper);

%% Quesstion 3 part c
var_power_v1 = sum(weights .* (power_v1 - mean_power_v1).^2) / sum(weights);
var_power_v2 = sum(weights .* (power_v2 - mean_power_v2).^2) / sum(weights);
var_combined_power = var_power_v1 + var_power_v2 + 2 * covariance_estimate;
std_combined_power = sqrt(var_combined_power);
fprintf('Variance of combined power output: %.4f MW^2\n', var_combined_power);
fprintf('Standard deviation of combined power output: %.4f MW\n', std_combined_power);

%% Question 3 part d

% Condition 1: Combined power > 9.5 MW
condition1 = (power_v1 + power_v2) > 9.5;  % Check condition for each sample
probability1 = sum(weights .* condition1) / sum(weights);

% Calculate 95% confidence interval for probability1
std_prob1 = sqrt(probability1 * (1 - probability1) / num_samples);
ci_prob1 = [probability1 - 1.96 * std_prob1, probability1 + 1.96 * std_prob1];

fprintf('Probability (Combined Power > 9.5 MW): %.4f\n', probability1);
fprintf('95%% Confidence Interval: [%.4f, %.4f]\n', ci_prob1(1), ci_prob1(2));

% Condition 2: Combined power < 9.5 MW
condition2 = (power_v1 + power_v2) < 9.5; 
probability2 = sum(weights .* condition2) / sum(weights);

% Calculate 95% confidence interval for probability2
std_prob2 = sqrt(probability2 * (1 - probability2) / num_samples);
ci_prob2 = [probability2 - 1.96 * std_prob2, probability2 + 1.96 * std_prob2];

fprintf('Probability (Combined Power < 9.5 MW): %.4f\n', probability2);
fprintf('95%% Confidence Interval: [%.4f, %.4f]\n', ci_prob2(1), ci_prob2(2));

condition3 = (power_v1 + power_v2) == 9.5;  % Equality condition
probability3 = sum(weights .* condition3) / sum(weights);
fprintf('Probability (Combined Power = 9.5 MW): %.4f\n', probability3);

% Sum of probabilities
sum_of_probabilities = probability1 + probability2 + probability3;
fprintf('Sum of probabilities: %.4f\n', sum_of_probabilities);

