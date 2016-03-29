% Script to graph the fraction of spectral energy that is correctly represented
%   in the DFT of Gaussian peaks. 

clear all

% Set directory information
path_parts = strsplit(pwd, 'CLM_paper');
filepath = [path_parts{1}, 'CLM_figures_and_data/Aliasing/'];

% Define sizes of each curve
sigma_values = [5, 10, 15, 20, 25, 30];
pixel_sizes = linspace(0, 50, 101); % gives row vector

% Calculate the fraction of spectral energy within bounds
fraction_correct = zeros(length(pixel_sizes), length(sigma_values));
for sigma_index = 1:length(sigma_values);
    fraction_vector = erf((pi * sigma_values(sigma_index)) ./ (sqrt(2) * pixel_sizes'));
    fraction_correct(:, sigma_index) = fraction_vector;
end

% Create graphs
figure
hold on
set(groot, 'defaultAxesLineStyleOrder', {'-', '--', ':', '-.'});
set(groot, 'defaultAxesColorOrder', [0, 0, 0]);
plot(pixel_sizes, fraction_correct);
ylim([0, 1.05]);
xlabel('Sampling distance (nm)');
ylabel('Fraction of signal');
title({'DFT of a Gaussian peak with standard deviation = \sigma', 'Fraction of signal correctly represented at a given sampling rate'}); 
legend_labels = cell(size(sigma_values));
for index = 1:length(sigma_values)
    legend_labels{index} = ['\sigma = ', num2str(sigma_values(index))]; 
end
legend(legend_labels, 'Location', 'southwest');
hold off
figure_name = [filepath, 'Gaussian_aliasing_fraction.png'];
print(gcf, figure_name, '-dpng');
delete(gcf);