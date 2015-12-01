% Script for graphing the results from get_timings_script

% Expects:
%   number_points_binning, number_points_pdf, number_points_parllel, 
%       number_points_direct, number_points_memory: The number of points 
%       evaluated in the correlation, independent variable, column vector 
%       of doubles.
%   binning_Fourier_times, pdf_Fourier_times, pdf_parallel_Fourier_times, 
%       pdist2_direct_times, custom_memory_direct_times: Time of evaluation 
%       for each method given the specified number of points, dependent 
%       variable, column vector of doubles.

% Make new figure
hfig = figure;

% Make axes and set options
haxes = axes('XScale', 'log', 'XLim', [1e2, 3e7], 'YScale', 'log', 'YLim', [2e-4, 2e3]);
xlabel('Number of points');
ylabel('Time (sec)');
title('Execution time of autocorrelation calculation algorithms')

% Plot each dataset
h1 = line(number_points_binning, binning_Fourier_times, 'Marker', 'o',...
    'Color', 'black', 'LineStyle', '-');
h2 = line(number_points_pdf, pdf_Fourier_times, 'Marker', '*',...
    'Color', 'black', 'LineStyle', '-');
h3 = line(number_points_pdf_parallel, pdf_parallel_Fourier_times, 'Marker', '.',...
    'Color', 'black', 'LineStyle', '-');
% h4 = line(number_points_pdf_MEX, pdf_MEX_Fourier_times, 'Marker', 's',...
%     'Color', 'black', 'LineStyle', '-');
% h5 = line(number_points_pdf_parallel_MEX, pdf_parallel_MEX_Fourier_times, 'Marker', 'd',...
%     'Color', 'black', 'LineStyle', '-');
h4 = line(number_points_binning, pdf_MEX_Fourier_times, 'Marker', 's',...
    'Color', 'black', 'LineStyle', '-');
h5 = line(number_points_binning, pdf_parallel_MEX_Fourier_times, 'Marker', 'd',...
    'Color', 'black', 'LineStyle', '-');
h6 = line(number_points_direct, pdist2_direct_times, 'Marker', 'x',...
    'Color', 'black', 'LineStyle', '-');
h7 = line(number_points_memory, custom_memory_direct_times, 'Marker', '+',...
    'Color', 'black', 'LineStyle', '-');

% Add Legend
legend('Fourier - binning','Fourier - pdf','Fourier - parallel pdf',...
    'Fourier - MEX pdf', 'Fourier - parallel MEX pdf', ...
    'Direct - MATLAB', 'Direct - memory saving','Location','southeast');

% Save figure
figure_path_parts = strsplit(pwd, 'CLM_paper');
figure_path = [figure_path_parts{1}, 'CLM_figures_and_data/'];
print(gcf, [figure_path, 'Figure1_v2.png'], '-dpng');

