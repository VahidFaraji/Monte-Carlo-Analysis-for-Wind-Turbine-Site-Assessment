% Estimate P(P(V) > 0) using different sampling methods
% Define parameters
num_samples = 100000;         % Number of samples
cut_in_speed = 3.5;           % Minimum operational wind speed (m/s)
cut_off_speed = 25;           % Maximum operational wind speed (m/s)
lambda_w = 10.6;              % Weibull scale parameter
k_w = 2.0;                    % Weibull shape parameter
alpha_gamma = 2;              % Gamma shape parameter
beta_gamma = 0.23;            % Gamma rate parameter

% Load power curve function (assumed in P(v))
load powercurve_V164.mat

% --- 1. Standard Monte Carlo Sampling (without truncation) ---
samples_mc = wblrnd(lambda_w, k_w, num_samples, 1);
power_outputs_mc = P(samples_mc);
prob_mc = sum(power_outputs_mc > 0) / num_samples;

% --- 2. Standard Monte Carlo Sampling (with truncation) ---
samples_mc_trunc = samples_mc(samples_mc >= cut_in_speed & samples_mc <= cut_off_speed);
power_outputs_mc_trunc = P(samples_mc_trunc);
prob_mc_trunc = sum(power_outputs_mc_trunc > 0) / length(samples_mc_trunc);

% --- 3. Antithetic Sampling ---
U = wblrnd(lambda_w, k_w, num_samples, 1);
U = U(U >= cut_in_speed & U <= cut_off_speed);
U_tilde = cut_off_speed - (U - cut_in_speed);
power_U = P(U);
power_U_tilde = P(U_tilde);
W = (power_U + power_U_tilde) / 2;
prob_antithetic = sum(W > 0) / length(W);

% --- 4. Importance Sampling ---
samples_gamma = gamrnd(alpha_gamma, 1/beta_gamma, num_samples, 1);
samples_gamma = samples_gamma(samples_gamma >= cut_in_speed & samples_gamma <= cut_off_speed);
weights = (k_w / lambda_w) .* (samples_gamma / lambda_w).^(k_w - 1) .* ...
          exp(-(samples_gamma / lambda_w).^k_w) ./ ...
          ((beta_gamma^alpha_gamma / gamma(alpha_gamma)) .* samples_gamma.^(alpha_gamma - 1) .* exp(-beta_gamma * samples_gamma));
power_outputs_is = P(samples_gamma);
prob_is = sum(weights .* (power_outputs_is > 0)) / sum(weights);

% --- Display Results ---
fprintf('Results Comparison:\n');
fprintf('1. Standard Monte Carlo (no truncation): Estimated Probability = %.4f\n', prob_mc);
fprintf('2. Standard Monte Carlo (with truncation): Estimated Probability = %.4f\n', prob_mc_trunc);
fprintf('3. Antithetic Sampling: Estimated Probability = %.4f\n', prob_antithetic);
fprintf('4. Importance Sampling: Estimated Probability = %.4f\n', prob_is);

% --- Variance Reduction Comparison ---
var_original = var(power_outputs_mc);
var_antithetic = var(W);
var_reduction_antithetic = (1 - (var_antithetic / var_original)) * 100;
fprintf('Variance Reduction (Antithetic Sampling): %.2f%%\n', var_reduction_antithetic);

%% Question 2 part e
P_power = wblcdf(25, lambda_w, k_w) - wblcdf(3.5, lambda_w, k_w);

%% Question 2 part f

fprintf('Estimated Probability = %.4f\n', P_power);
availability_factor = mean(P_power);
fprintf('Availibility Factor = %.4f\n', availability_factor);

% Define maximum power output
max_power_output = 9.5e6;  % 9.5 MW in Watts

% Power output estimates for each month (updated with new results)
monthly_mean_powers = [
    5215956.50,  % January
    4722317.16,  % February
    4430548.13,  % March
    3698436.87,  % April
    3566077.92,  % May
    3776269.17,  % June
    3563904.23,  % July
    3773043.51,  % August
    4366806.23,  % September
    4865135.00,  % October
    5221153.94,  % November
    5215291.00   % December
];

% Calculate Capacity Factor for each month
capacity_factors = monthly_mean_powers / max_power_output;

% Calculate the average Capacity Factor
average_capacity_factor = mean(capacity_factors);

fprintf('Capacity Factor = %.4f\n', average_capacity_factor);


