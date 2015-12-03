
% Load data
load('fractional_discrepency_curves_SN10.mat')

% Get mean and standard deviation of the replicates
[approx_FD_mean, approx_FD_stdev] = calc_mean_stdev_FD_matrix(approx_SSQ_results);
[ideal_FD_mean, ideal_FD_stdev] = calc_mean_stdev_FD_matrix(ideal_SSQ_results);

% plot curves
figure
hold on
%errorbar(repmat(event_fractions, 1, size(approx_FD_mean, 2)), approx_FD_mean, approx_FD_stdev, 'k-d');
%approx_zeros = approx_FD_mean(end, :);
ideal_zeros = ideal_FD_mean(end, :);
%plot(repmat(event_fractions, 1, size(approx_FD_mean, 2)), approx_FD_mean./repmat(approx_zeros, size(approx_FD_mean, 1), 1));
plot(repmat(event_fractions, 1, size(ideal_FD_mean, 2)), ideal_FD_mean./repmat(ideal_zeros, size(ideal_FD_mean, 1), 1));
xlabel('Fraction of events included')
ylabel('Sum of Squares')
hold off
% 
% % plot curves
% figure
% hold on
% approx_zero_mean = mean(approx_SSQ_results(end, :));
% approx_relative_SSQ = approx_SSQ_results/approx_zero_mean;
% ideal_zero_mean = mean(ideal_SSQ_results(end, :));
% ideal_relative_SSQ = ideal_SSQ_results/ideal_zero_mean;
% errorbar(event_fractions, mean(approx_relative_SSQ, 2), std(approx_relative_SSQ, 0, 2), 'r-*');
% errorbar(event_fractions, mean(ideal_relative_SSQ, 2), std(ideal_relative_SSQ, 0, 2), 'k-d');
% xlabel('Fraction of events included')
% ylabel('Relative Sum of Squares')
% hold off
