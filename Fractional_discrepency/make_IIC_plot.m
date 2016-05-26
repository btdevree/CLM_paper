function make_IIC_plot(filename, fraction_vector, event_numbers, IIC_mean, IIC_stdev)
%MAKE_NONPARAMETRIC_PLOT Creates a .png image of a information improvement
% characteristic curve. 
%
% Inputs:
%   filename: string, full filename and path of image to be created
%   fraction_vector: column vector of event fractions
%   event_numbers: column vector of the total number of events in the
%       dataset used to generate each curve.
%   IIC_mean: information improvement curve matrix with each curve arranged
%       as a column.
%   IIC_stdev: standard deviation matrix of the IIC_mean values, with each
%       curve arranged as a column. Optional, default = [], no error bars.
% Output:
%   Image saved with given filename, no return arguments.

% Set defaults
if nargin < 5; IIC_stdev = []; end;

% Create new figure
hfig = figure('Units', 'pixels', 'Position', [150, 150, 650, 500], 'PaperPositionMode', 'auto', 'InvertHardcopy', 'off', 'Color', [1, 1, 1]);

% Define linespec options to cycle through
linespec_options = {'-+k', '-ok', '-*k', '-sk', '-xk', '-dk', ':+k', ':ok', ':*k', ':sk', ':xk', ':dk', '--+k', '--ok', '--*k', '--sk', '--xk', '--dk'};

% Create clear axes for plotting shapes
haxes = axes('Units', 'pixels', 'Position', [70, 50, 500, 400]);

% Loop through each curve
hold on
plot_handles = gobjects(size(IIC_mean, 2), 1);
legend_labels = cell(size(IIC_mean, 2), 1);
for curve_index = 1:size(IIC_mean, 2);
    
    % Plot the curve
    current_linespec = linespec_options{mod(curve_index - 1, length(linespec_options))+1};
    plot_handles(curve_index) = plot(fraction_vector, IIC_mean(:, curve_index), current_linespec);
    if ~isempty(IIC_stdev)
        errorbar(fraction_vector, IIC_mean(:, curve_index), IIC_stdev(:, curve_index), current_linespec);
    end
    
    % Add to legend
    legend_labels{curve_index} = num2str(event_numbers(curve_index),'%.1E'); 
end
hold off

% Finish graph
set(haxes, 'Xlim', [0, 1], 'Ylim', [0, 1]);
hlegend = legend(plot_handles, legend_labels, 'Location', 'southeast'); % Errorbar legend markers don't show up, use Plot curves for legend 
if verLessThan('matlab','8.4') % Stupid MATLAB changed lots of stuff about graphics
    hlt = get(hlegend, 'title');
    set(hlt, 'String', 'Total Number of Events');
else
    hlt = text('Parent', hlegend.DecorationContainer, 'String', 'Total Number of Events', 'HorizontalAlignment', 'center', ...
        'VerticalAlignment', 'bottom', 'Position', [0.5, 1.05, 0], 'Units', 'normalized', 'FontSize', 10);
end
title('Information Improvement Characteristic Curves', 'FontSize', 16);
xlabel('Fraction of Localization Events Included','FontSize', 12);
ylabel('Normalized Information Improvement','FontSize', 12);

% Save image
print(hfig, filename, '-dpng');

% Delete figure
delete(hfig);
end

