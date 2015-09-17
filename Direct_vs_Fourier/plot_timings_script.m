% Script for graphing the results from get_timings_script

% Expects:
%   number_points_vector, number_points_direct, number_points_memory:
%       The number of points evaluated in the correlation, independent
%       variable, column vector of doubles.
%   binning_Fourier_times, pdf_Fourier_times, pdist2_direct_times,
%       custom_memory_direct_times: Time of evaluation for each method
%       given the specified number of points, dependent variable, column
%       vector of doubles.

% Make new figure
hfig = figure;

% Make axes and set options
haxes = axes('XScale', 'log', 'XLim', [1e2, 1e7], 'YScale', 'log', 'YLim', [1e-4, 1100]);
xlabel('Number of points');
ylabel('Time (sec)');
title('Execution time of autocorrelation calculation algorithms')

% Plot each dataset
h1 = line(number_points_vector, binning_Fourier_times, 'Marker', 'o',...
    'Color', 'black', 'LineStyle', '-');
h2 = line(number_points_vector, pdf_Fourier_times, 'Marker', '+',...
    'Color', 'black', 'LineStyle', '-');
h3 = line(number_points_direct, pdist2_direct_times, 'Marker', 'x',...
    'Color', 'black', 'LineStyle', '-');
h4 = line(number_points_memory, custom_memory_direct_times, 'Marker', '+',...
    'Color', 'black', 'LineStyle', '-');

% Add Legend



