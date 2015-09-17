% Script for graphing the results from get_timings_script

% Expects:
%   number_points_vector, number_points_direct, number_points_memory:
%       The number of points evaluated in the correlation, independent
%       variable, column vector of doubles.
%   binning_Fourier_times, pdf_Fourier_times, pdist2_direct_times,
%       custom_memory_direct_times: Time of evaluation for each method
%       given the specified number of points, dependent variable, column
%       vector of doubles.

% testdata
number_points_vector = [1e2; 3e2; 1e3; 3e3; 1e4; 3e4; 1e5; 3e5; 1e6; 3e6; 1e7];
binning_Fourier_times = [12; 13; 14; 15; 16; 17; 18; 20; 22; 25; 30];


% Make new figure
hfig = figure;

% Make axes and set options
haxes = axes('XScale', 'log', 'XLim', [1e2, 1e7], 'YLim', [0, 30],...
    'XLabel', 'Number of points', 'YLabel', 'Time (sec)');

% Plot each dataset
h1 = plot(haxes, number_points_vector, binning_Fourier_times, 'Marker', 'o',...
    'Color', 'black', 'LineStyle', '-');
h2 = plot(haxes, number_points_vector, binning_Fourier_times-12, 'Marker', '+',...
    'Color', 'black', 'LineStyle', '-');
h3 = plot(haxes, number_points_vector, binning_Fourier_times-5, 'Marker', 'x',...
    'Color', 'black', 'LineStyle', '-');
h4 = plot(haxes, number_points_vector, binning_Fourier_times-10, 'Marker', '+',...
    'Color', 'black', 'LineStyle', '-');

% Add Legend

print(gcf, 'average_xcor', '-dpng');
clf

