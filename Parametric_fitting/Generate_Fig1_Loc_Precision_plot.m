% Script to generate the localization precision graph used in figure 1 & 2 
% Run in Parametric_fitting folder

% Get path for output
figure_path_parts = strsplit(pwd, 'CLM_paper');
figure_path = [figure_path_parts{1}, 'CLM_figures_and_data/Fig1/'];

% Load in parameters structure
load('parameters_Fig1.mat'); % loads 'params' into local namespace

% Get localizaiton precision
precision = params.STORM_precision;

% Generate distance vector
dist_vector = [-100:1:100];

% Generate the normal pdf
pdf_vector = normpdf(dist_vector, 0, 25);

% Create new figure
hfig = figure('Units', 'pixels', 'Position', [100, 100, 700, 500], 'PaperPositionMode', 'auto', 'InvertHardcopy', 'off', 'Color', [1, 1, 1]);

% Create clear axes for plotting shapes
haxes = axes('Units', 'pixels', 'Position', [70, 50, 600, 400]);
plot(haxes, dist_vector, pdf_vector, 'LineWidth', 4, 'Color', [0, 0, 0]);
set(haxes, 'Xlim', [dist_vector(1), dist_vector(end)], 'Ylim', [0, 0.02], 'Color', 'none');
title('Localization Precision', 'FontSize', 16);
xlabel('Distance (nanometers)','FontSize', 12);
ylabel('Probability Density','FontSize', 12);

% Save image
filename = [figure_path, 'localization_precision'];
print(hfig, filename, '-dpng');

% Cleanup
delete(hfig);