
% Load data
load('test_SSQ_data_v2');

% plot curves
figure
hold on
errorbar(event_fractions, mean(approx_SSQ_results, 2), std(approx_SSQ_results, 0, 2), 'r-*');
errorbar(event_fractions, mean(ideal_SSQ_results, 2), std(ideal_SSQ_results, 0, 2), 'k-d');
xlabel('Fraction of events included')
ylabel('Sum of Squares')
hold off

% plot curves
figure
hold on
approx_zero_mean = mean(approx_SSQ_results(end, :));
approx_relative_SSQ = approx_SSQ_results/approx_zero_mean;
ideal_zero_mean = mean(ideal_SSQ_results(end, :));
ideal_relative_SSQ = ideal_SSQ_results/ideal_zero_mean;
errorbar(event_fractions, mean(approx_relative_SSQ, 2), std(approx_relative_SSQ, 0, 2), 'r-*');
errorbar(event_fractions, mean(ideal_relative_SSQ, 2), std(ideal_relative_SSQ, 0, 2), 'k-d');
xlabel('Fraction of events included')
ylabel('Relative Sum of Squares')
hold off
