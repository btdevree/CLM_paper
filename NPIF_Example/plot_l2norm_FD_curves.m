
% Load data
load('fractional_discrepency_curves_l2norm_SN10.mat')

% Get mean and standard deviation of the replicates
[approx_FD_mean, approx_FD_stdev] = calc_mean_stdev_FD_matrix(approx_SSQ_results);
[ideal_FD_mean, ideal_FD_stdev] = calc_mean_stdev_FD_matrix(ideal_SSQ_results);

% plot curves
figure
hold on
ideal_zeros = ideal_FD_mean(end, :);
plot(repmat(event_fractions, 1, size(ideal_FD_mean, 2)), ideal_FD_mean./repmat(ideal_zeros, size(ideal_FD_mean, 1), 1));
xlabel('Fraction of events included')
ylabel('Relative Sum of Squares')
title('L2 norm ideal relative fractional discrepency');
hold off

figure
hold on
approx_zeros = approx_FD_mean(end, :);
plot(repmat(event_fractions, 1, size(approx_FD_mean, 2)), approx_FD_mean./repmat(approx_zeros, size(approx_FD_mean, 1), 1));
xlabel('Fraction of events included')
ylabel('Relative Sum of Squares')
title('L2 norm approximated relative fractional discrepency');
hold off

figure
hold on
plot(repmat(event_fractions, 1, size(ideal_FD_mean, 2)), ideal_FD_mean);
xlabel('Fraction of events included')
ylabel('Sum of Squares')
title('L2 norm ideal fractional discrepency');
hold off

figure
hold on
plot(repmat(event_fractions, 1, size(approx_FD_mean, 2)), approx_FD_mean);
xlabel('Fraction of events included')
ylabel('Sum of Squares')
title('L2 norm approximated fractional discrepency');
hold off

