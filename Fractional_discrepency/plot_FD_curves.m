
% Load data
load('fractional_discrepency_curves_ssq_SN10.mat')
SN_ratio = 10;
method_name = 'SSQ';

% Get mean and standard deviation of the replicates
[approx_FD_mean, approx_FD_stdev] = calc_mean_stdev_FD_matrix(approx_SSQ_results);
[ideal_FD_mean, ideal_FD_stdev] = calc_mean_stdev_FD_matrix(ideal_SSQ_results);

% Get figure path
figure_path_parts = strsplit(pwd, 'CLM_paper');
figure_path = [figure_path_parts{1}, 'CLM_figures_and_data/FD/'];

% plot curves
figure
hold on
set(gcf, 'Position', [100, 100, 800, 600]);
set(gcf, 'PaperPositionMode', 'auto');
ideal_zeros = ideal_FD_mean(end, :);
plot(repmat(event_fractions, 1, size(ideal_FD_mean, 2)), ideal_FD_mean./repmat(ideal_zeros, size(ideal_FD_mean, 1), 1));
xlabel('Fraction of events included')
ylabel('Relative Sum of Squares')
title('Old SSQ ideal relative fractional discrepency');
hold off
print([figure_path, 'Old_Ideal_relative_FD_SN', num2str(SN_ratio), '_', method_name, '.png'], '-dpng');
delete(gcf)

figure
hold on
set(gcf, 'Position', [100, 100, 800, 600]);
set(gcf, 'PaperPositionMode', 'auto');
approx_zeros = approx_FD_mean(end, :);
plot(repmat(event_fractions, 1, size(approx_FD_mean, 2)), approx_FD_mean./repmat(approx_zeros, size(approx_FD_mean, 1), 1));
xlabel('Fraction of events included')
ylabel('Relative Sum of Squares')
title('Old SSQ approximated relative fractional discrepency');
hold off
print([figure_path, 'Old_Approx_relative_FD_SN', num2str(SN_ratio), '_', method_name, '.png'], '-dpng');
delete(gcf)

figure
hold on
set(gcf, 'Position', [100, 100, 800, 600]);
set(gcf, 'PaperPositionMode', 'auto');
plot(repmat(event_fractions, 1, size(ideal_FD_mean, 2)), ideal_FD_mean);
xlabel('Fraction of events included')
ylabel('Sum of Squares')
title('Old SSQ ideal fractional discrepency');
hold off
print([figure_path, 'Old_Ideal_raw_FD_SN', num2str(SN_ratio), '_', method_name, '.png'], '-dpng');
delete(gcf)

figure
hold on
set(gcf, 'Position', [100, 100, 800, 600]);
set(gcf, 'PaperPositionMode', 'auto');
plot(repmat(event_fractions, 1, size(approx_FD_mean, 2)), approx_FD_mean);
xlabel('Fraction of events included')
ylabel('Sum of Squares')
title('Old SSQ approximated fractional discrepency');
hold off
print([figure_path, 'Old_Approx_raw_FD_SN', num2str(SN_ratio), '_', method_name, '.png'], '-dpng');
delete(gcf)
