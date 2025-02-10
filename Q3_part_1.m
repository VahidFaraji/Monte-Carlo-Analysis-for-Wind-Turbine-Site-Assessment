clear; clc; rng(0);

load('powercurve_V164.mat'); 

% --- Step 2: Define Weibull distribution parameters ---
lambda_weibull = 9.13;  % Scale parameter (m/s)
k_weibull = 1.96;       % Shape parameter

% --- Step 3: Define truncation limits ---
v_min = 3.5;  % Cut-in speed (m/s)
v_max = 25.0; % Cut-off speed (m/s)

% Calculate CDF values at truncation points
F_low = wblcdf(v_min, lambda_weibull, k_weibull);
F_high = wblcdf(v_max, lambda_weibull, k_weibull);

% Define truncated Weibull PDF
f_truncated_weibull = @(v) wblpdf(v, lambda_weibull, k_weibull) / (F_high - F_low);

% --- Step 4: Define the Gamma instrumental distribution ---
alpha_gamma = 3.85;  % Shape parameter
beta_gamma = 0.475;   % Rate parameter

% Define Gamma PDF and sampling function
g_gamma = @(v) (beta_gamma^alpha_gamma / gamma(alpha_gamma)) * v.^(alpha_gamma - 1) .* exp(-beta_gamma * v);
sample_gamma = @(n) gamrnd(alpha_gamma, 1 / beta_gamma, [n, 1]);

% --- Step 5: Perform importance sampling with truncation ---
n_samples = 500000;     % Number of samples
v_samples = sample_gamma(n_samples);  % Generate samples from Gamma distribution

% Truncate samples to the range [v_min, v_max]
v_samples = v_samples(v_samples >= v_min & v_samples <= v_max);

% Recalculate sample size after truncation
n_valid_samples = length(v_samples);

% Calculate importance weights and power output
weights = f_truncated_weibull(v_samples) ./ g_gamma(v_samples);
power_outputs = P(v_samples);  % Evaluate power output using the power curve
weighted_outputs = power_outputs .* weights;

% --- Step 6: Estimate expected power output ---
expected_power = mean(weighted_outputs);
variance_estimate = var(weighted_outputs) / n_valid_samples;
confidence_interval = [expected_power - 1.96 * sqrt(variance_estimate), ...
                       expected_power + 1.96 * sqrt(variance_estimate)];
% Assuming 'weights' is a vector of importance weights
ess = (sum(weights)^2) / sum(weights.^2);
ess_percentage = (ess / n_valid_samples) * 100;

fprintf('Effective Sample Size: %.2f (%.2f%% of total samples)\n', ess, ess_percentage);

% --- Step 7: Display results ---
fprintf('Expected Power Output: %.2f W\n', expected_power);
fprintf('95%% Confidence Interval: [%.2f W, %.2f W]\n', confidence_interval(1), confidence_interval(2));

% --- Step 8: Plot distributions for comparison ---
figure;
v_range = linspace(0, 30, 1000);
plot(v_range, f_truncated_weibull(v_range), 'b-', 'LineWidth', 2); hold on;
plot(v_range, g_gamma(v_range), 'r--', 'LineWidth', 2);
xlabel('Wind Speed (m/s)');
ylabel('PDF');
legend('Truncated Weibull Distribution', 'Gamma Distribution');
title('Comparison of Target and Instrumental Distributions');
grid on;
