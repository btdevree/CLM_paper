% Make a figure from the completeness index calculations made in
% make_part_C_data

% For only one SN ratio only

% Load the CI data
SN_ratio = 10;
load(['NPIF_part_C_CIdata_ssq99_SN', num2str(SN_ratio), '.mat']);

% Pick the method used

% approx_CI = approx_CI_method1;
% ideal_CI = ideal_CI_method1;
% method_name = 'NVI-100';

% approx_CI = approx_CI_method2;
% ideal_CI = ideal_CI_method2;
% method_name = 'SSQ';

% approx_CI = approx_CI_method1;
% ideal_CI = ideal_CI_method1;
% method_name = 'SSQ-100';

approx_CI = approx_CI_method2;
ideal_CI = ideal_CI_method2;
method_name = 'l2-99';


% Collapse the replicates and pseudoreplicates for the approximate CI data
approx_CI_values = permute(reshape(permute(approx_CI, [4, 3, 2, 1]),...
    [size(approx_CI, 3) * size(approx_CI, 4), size(approx_CI, 2), size(approx_CI, 1)]), [3, 2, 1]);

% Get average and stdev of the approximate CI values
approx_CI_mean = mean(approx_CI_values, 3);
approx_CI_median = median(approx_CI_values, 3);
approx_CI_stdev = std(approx_CI_values, 0, 3);
approx_CI_mad = mad(approx_CI_values, 1, 3);

% Get average and stdev of the ideal CI values
ideal_CI_mean = mean(ideal_CI, 3);
ideal_CI_stdev = std(ideal_CI, 0, 3);

% Make figure
hfig = figure;
set(gcf, 'Position', [100, 100, 800, 600]);
set(gcf, 'PaperPositionMode', 'auto');
hold on
errorbar(true_event_numbers, ideal_CI_mean, ideal_CI_stdev, 'ko-', 'markerfacecolor', 'k', 'markersize', 3);
%errorbar(true_event_numbers, approx_CI_mean, approx_CI_stdev, 'ko-', 'markerfacecolor', 'k', 'markersize', 3);
errorbar(true_event_numbers, approx_CI_median, approx_CI_mad, 'ro-', 'markerfacecolor', 'k', 'markersize', 3);
xlabel('Number of events');
ylabel('Completeness index');
set(gca,'XScale','log');
ylim([0, 1]);
xlim([1e3, 1e7]);
title({'Ideal vs Approximated Completeness Index', ['S/N ratio = ', num2str(SN_ratio), '; Method = ', method_name]});
hold off

figure_path_parts = strsplit(pwd, 'CLM_paper');
figure_path = [figure_path_parts{1}, 'CLM_figures_and_data/'];
print([figure_path, 'NPIF_SN', num2str(SN_ratio), '_', method_name, '.png'], '-dpng');