% Script to create informal improvement characteristic curves


% Load data
load('fractional_discrepency_curves_doublecheck_ssq_SN10.mat')
SN_ratio = 10;
method_name = 'OLD_SSQ';

% Get mean and standard deviation of the replicates
[approx_FD_mean, approx_FD_stdev] = calc_mean_stdev_FD_matrix(approx_SSQ_results);
[ideal_FD_mean, ideal_FD_stdev] = calc_mean_stdev_FD_matrix(ideal_SSQ_results);

% Calculate information improvement curves
% Get max and min discrepency values
approx_max_discrepency = approx_FD_mean(end, :);
approx_min_discrepency = approx_FD_mean(1, :);

ideal_max_discrepency = ideal_FD_mean(end, :);
ideal_min_discrepency = ideal_FD_mean(1, :);

% Transform to fractional improvement of the discrepency 
approx_delta_discrepency = approx_FD_mean - repmat(approx_min_discrepency, size(approx_FD_mean, 1), 1);
approx_max_delta = approx_max_discrepency - approx_min_discrepency;
approx_info_improvement = 1-(approx_delta_discrepency ./ repmat(approx_max_delta, size(approx_FD_mean, 1), 1));

ideal_delta_discrepency = ideal_FD_mean - repmat(ideal_min_discrepency, size(ideal_FD_mean, 1), 1);
ideal_max_delta = ideal_max_discrepency - ideal_min_discrepency;
ideal_info_improvement = 1-(ideal_delta_discrepency ./ repmat(ideal_max_delta, size(ideal_FD_mean, 1), 1));

% Calculate AUC and TCI/ICI/ECI
dFrac = event_fractions(1:end-1) - event_fractions(2:end);
midpoint_II_approx = (approx_info_improvement(2:end, :) + approx_info_improvement(1:end-1, :))/2;
AUC_approx = sum(repmat(dFrac, 1, size(midpoint_II_approx, 2)) .* midpoint_II_approx, 1);
ECI = (2*AUC_approx - 1)';

midpoint_II_ideal = (ideal_info_improvement(2:end, :) + ideal_info_improvement(1:end-1, :))/2;
AUC_ideal = sum(repmat(dFrac, 1, size(midpoint_II_ideal, 2)) .* midpoint_II_ideal, 1);
ICI = (2*AUC_ideal - 1)';

TCI = (1 - (ideal_min_discrepency ./ ideal_max_discrepency))';

% Get figure path
figure_path_parts = strsplit(pwd, 'CLM_paper');
figure_path = [figure_path_parts{1}, 'CLM_figures_and_data/FD/'];

% plot curves
figure
hold on
set(gcf, 'Position', [100, 100, 800, 600]);
set(gcf, 'PaperPositionMode', 'auto');
plot(repmat(event_fractions, 1, size(ideal_FD_mean, 2)), ideal_info_improvement);
xlabel('Fraction of events included')
ylabel([method_name, ' fractional information improvement'])
title([method_name, ' Ideal Information Improvement Charactistic curve']);
hold off
print([figure_path, 'New_Ideal_IIC_SN', num2str(SN_ratio), '_', method_name, '.png'], '-dpng');
delete(gcf)

figure
hold on
set(gcf, 'Position', [100, 100, 800, 600]);
set(gcf, 'PaperPositionMode', 'auto');
plot(repmat(event_fractions, 1, size(approx_FD_mean, 2)), approx_info_improvement);
xlabel('Fraction of events included')
ylabel([method_name, ' fractional information improvement'])
title([method_name, ' Approximated Information Improvement Charactistic curve']);
hold off
print([figure_path, 'New_Approx_IIC_SN', num2str(SN_ratio), '_', method_name, '.png'], '-dpng');
delete(gcf)

figure
hold on
set(gcf, 'Position', [100, 100, 800, 600]);
set(gcf, 'PaperPositionMode', 'auto');
plot(repmat(true_event_numbers, 1, 3), [TCI, ICI, ECI]);
xlabel('Number of events in image')
ylabel([method_name, ' Ideal or Estimated Completeness Index'])
title([method_name, ' Ideal or Estimated Completeness Index vs. event number']);
legend('TCI', 'ICI', 'ECI', 'Location', 'southeast');
set(gca,'XScale','log');
hold off
print([figure_path, 'New_CI_SN', num2str(SN_ratio), '_', method_name, '.png'], '-dpng');
delete(gcf)

figure
hold on
set(gcf, 'Position', [100, 100, 800, 600]);
set(gcf, 'PaperPositionMode', 'auto');
plot(TCI, ECI);
xlabel('Theoretical Completeness Index')
ylabel('Estimated Completeness Index')
title([method_name, ' Theoretical vs Estimated Completeness Index']);
hold off
print([figure_path, 'New_TCI_vs_ECI', num2str(SN_ratio), '_', method_name, '.png'], '-dpng');
delete(gcf)
