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

% Define Availability Factor (assume full availability for simplicity)
availability_factors = ones(1, 12);  % 100% availability for all months

% Calculate the average Availability Factor
average_availability_factor = mean(availability_factors);

% Display results
fprintf('Average Capacity Factor: %.2f%%\n', average_capacity_factor * 100);
fprintf('Average Availability Factor: %.2f%%\n', average_availability_factor * 100);

% Check if the site is suitable
if average_capacity_factor * 100 >= 20 && average_capacity_factor * 100 <= 40 && average_availability_factor * 100 >= 90
    disp('The site is suitable for wind turbine installation.');
else
    disp('The site may not be suitable for wind turbine installation.');
end
