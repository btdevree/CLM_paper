function make_ECI_TCI_plot(filename, event_numbers, ECI_mean, TCI_mean, ECI_stdev, TCI_stdev)
%MAKE_NONPARAMETRIC_PLOT Creates a .png image of a information improvement
% characteristic curve. 
%
% Inputs:
%   filename: string, full filename and path of image to be created
%   event_numbers: column vector of the total number of events in the
%       dataset used to generate each curve.
%   ECI/TCI_mean: column vectors of experimental and theoretical
%       completeness index mean values, a seperate cuve for each 
%   ECI/TCI_stdev: column vector of experimental and theoretical
%       completeness index standard deviation values. Optional, default =
%       [] (no error bars).
%   IIC_stdev: standard deviation matrix of the IIC_mean values, with each
%       curve arranged as a column. Optional, default = [], no error bars.
% Output:
%   Image saved with given filename, no return arguments.

% Set defaults
if nargin < 5; ECI_stdev = []; end;
if nargin < 6; TCI_stdev = []; end;

% Create new figure
hfig = figure('Units', 'pixels', 'Position', [150, 150, 650, 400], 'PaperPositionMode', 'auto', 'InvertHardcopy', 'off', 'Color', [1, 1, 1]);

% Define linespec options to cycle through
linespec_options = {'-+k', '-ok', '-*k', '-sk', '-xk', '-dk', ':+k', ':ok', ':*k', ':sk', ':xk', ':dk', '--+k', '--ok', '--*k', '--sk', '--xk', '--dk'};

% Create clear axes for plotting shapes
haxes = axes('Units', 'pixels', 'Position', [70, 50, 500, 300]);

% Plot curves
plot_handles = gobjects(2, 1);
legend_labels = cell(2, 1);
hold on

% ECI
plot_handles(1) = plot(event_numbers, ECI_mean, linespec_options{1});
if ~isempty(ECI_stdev)
    errorbar(event_numbers, ECI_mean, ECI_stdev, linespec_options{1});
end
legend_labels{1} = num2str('Experimental Completeness Index');

% TCI
plot_handles(2) = plot(event_numbers, TCI_mean, linespec_options{2});
if ~isempty(ECI_stdev)
    errorbar(event_numbers, TCI_mean, TCI_stdev, linespec_options{2});
end
legend_labels{2} = num2str('Theoretical Completeness Index');

% Finish graph
set(haxes, 'Xscale', 'log', 'Ylim', [0, 1]);
legend(plot_handles, legend_labels, 'Location', 'southeast'); % Errorbar legend markers don't show up, use Plot curves for legend 
title('Completeness Index vs. Localization Event Number', 'FontSize', 16);
xlabel('Number of Localization Events', 'FontSize', 12);
ylabel('Completeness Index, Sum of Squares method', 'FontSize', 12);
hold off

% Save image
print(hfig, filename, '-dpng');

% Delete figure
delete(hfig);
end

