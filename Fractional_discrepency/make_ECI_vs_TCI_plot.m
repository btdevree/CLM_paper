function make_ECI_vs_TCI_plot(filename, ECI_mean, TCI_mean, method_list, legend_method_list, ECI_stdev, TCI_stdev)
%MAKE_ECI_VS_TCI_PLOT Creates a .png image of the theoretical vs
% experimental completeness indices
%
% Inputs:
%   filename: string, full filename and path of image to be created
%   ECI/TCI_mean: column vectors of experimental and theoretical
%       completeness index mean values, a seperate cuve for each
%   method_list: cell array of strings representing the discrepency method 
%       used to calculate each pair of TCI and ECI curves
%   legend_method_list: cell array of strings with the names of the 
%       discrepency method for the plot legend 
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
hfig = figure('Units', 'pixels', 'Position', [150, 150, 540, 500], 'PaperPositionMode', 'auto', 'InvertHardcopy', 'off', 'Color', [1, 1, 1]);

% Define linespec options to cycle through
linespec_options = {'-xk', '-ok', '-*k', '-sk', '-+k', '-dk', ':xk', ':ok', ':*k', ':sk', ':+k', ':dk', '--xk', '--ok', '--*k', '--sk', '--+k', '--dk'};

% Create clear axes for plotting shapes
haxes = axes('Units', 'pixels', 'Position', [70, 50, 400, 400]);

% Plot curves
plot_handles = gobjects(length(method_list), 1);
hold on

% Loop through each method
for method_cell_index = 1:length(method_list)
    plot_handles(method_cell_index) = plot(TCI_mean{method_cell_index}, ECI_mean{method_cell_index}, linespec_options{method_cell_index});
    if ~isempty(ECI_stdev) && ~isempty(TCI_stdev)
        errorbar(TCI_mean{method_cell_index}, ECI_mean{method_cell_index}, ECI_stdev{method_cell_index}, linespec_options{method_cell_index});
        herrorbar(TCI_mean{method_cell_index}, ECI_mean{method_cell_index}, TCI_stdev{method_cell_index}, linespec_options{method_cell_index});
    end
end
  
% Finish graph
set(haxes, 'Xlim', [0, 1], 'Ylim', [0, 1]);
legend(plot_handles, legend_method_list, 'Location', 'northwest'); % Errorbar legend markers don't show up, use Plot curves for legend 
title('Theoretical vs. Experimental Completeness Index', 'FontSize', 16);
xlabel('Theoretical Completeness Index', 'FontSize', 12);
ylabel('Experimental Completeness Index', 'FontSize', 12);
hold off

% Save image
print(hfig, filename, '-dpng');

% Delete figure
delete(hfig);
end

