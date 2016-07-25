function graph_all_tests(test_directory, output_directory)
%GRAPH_all_tests Graph results from all test_summary_*.mat files in the
%test_directory. 
%
% Requires functions from Common folder.
%   
% Input:
%   test_directory: string, filepath for the test folder where the test
%       archive and response are kept. Optional, default = current 
%       directory.
%   output_directory: string, filepath for saving the test results.
%       Optional, default = test_directory.
% Output:
%   Graphs are saved to the output directory as .png files.

% Set defaults
if nargin < 1; test_directory = ''; end;
if nargin < 2; output_directory = test_directory; end;

% File files
if strcmp(test_directory, '')
    summary_paths = dir('test_summary_*.mat');
else
    summary_paths = find([test_directory, '/test_summary_*.mat']);
end

% Special case if only one file
if size(summary_paths, 1) == 1
    load([test_directory, '/', summary_paths.name]);
    summary_struct = test_summary;

% If more than 1 file, add them together    
elseif size(summary_paths, 1) >= 2
    
    % Load first file into overall summary structure
    load([test_directory, '/', summary_paths(1).name]);
    summary_struct = test_summary;
    
    % Concatenate with other summaries
    for summary_index = 2:size(summary_paths, 1);
        load([test_directory, '/', summary_paths(summary_index).name]); % Loads test_summary

        % Add summaries together into one big summary
        summary_struct = add_summaries(summary_struct, test_summary, 1);
    end
    
% No files, give error
else
    error('No test_summary files found');
end

% ---- Regions -----

% Rename for convenience
s = summary_struct.region;

% Found accuracy image
% Get vectors for the ECI
ECI_cr_0_5 = select_summary_subset(s.found_regions, 'ECI', 'contrast_ratios', 0.5);
ECI_cr_1_5 = select_summary_subset(s.found_regions, 'ECI', 'contrast_ratios', 1.5);
ECI_cr_4 = select_summary_subset(s.found_regions, 'ECI', 'contrast_ratios', 4);

% Get area vectors and calculate relative absolute area
abs_area_cr_0_5 = select_summary_subset(s.found_regions, 'abs_area', 'contrast_ratios', 0.5);
abs_area_cr_1_5 = select_summary_subset(s.found_regions, 'abs_area', 'contrast_ratios', 1.5);
abs_area_cr_4 = select_summary_subset(s.found_regions, 'abs_area', 'contrast_ratios', 4);
true_area_cr_0_5 = select_summary_subset(s.found_regions, 'true_area', 'contrast_ratios', 0.5);
true_area_cr_1_5 = select_summary_subset(s.found_regions, 'true_area', 'contrast_ratios', 1.5);
true_area_cr_4 = select_summary_subset(s.found_regions, 'true_area', 'contrast_ratios', 4);
relative_area_cr_0_5 = abs_area_cr_0_5 ./ true_area_cr_0_5;
relative_area_cr_1_5 = abs_area_cr_1_5 ./ true_area_cr_1_5;
relative_area_cr_4 = abs_area_cr_4 ./ true_area_cr_4;

% Create new figure
hfig = figure('Units', 'pixels', 'Position', [100, 100, 700, 500], 'PaperPositionMode', 'auto', 'InvertHardcopy', 'off', 'Color', [1, 1, 1]);

% Plot scatter plot of each contrast ratio
haxes = axes('Units', 'pixels', 'Position', [70, 50, 600, 400]);
scatter(haxes, ECI_cr_0_5, relative_area_cr_0_5, 'Marker', 'o', 'MarkerEdgeColor', 'k');
hold on
scatter(haxes, ECI_cr_1_5, relative_area_cr_1_5, 'Marker', '^', 'MarkerEdgeColor', 'k');
scatter(haxes, ECI_cr_4, relative_area_cr_4, 'Marker', 's', 'MarkerEdgeColor', 'k');
set(haxes, 'Xlim', [0, 1]);
legend_labels = {'1:1.5', '1:2.5', '1:5'};
hlegend = legend(legend_labels, 'Location', 'best');
if verLessThan('matlab','8.4') % Stupid MATLAB changed lots of stuff about graphics
    hlt = get(hlegend, 'title');
    set(hlt, 'String', 'Contrast Ratio');
else
    hlt = text('Parent', hlegend.DecorationContainer, 'String', 'Contrast Ratio', 'HorizontalAlignment', 'center', ...
        'VerticalAlignment', 'bottom', 'Position', [0.5, 1.05], 'Units', 'normalized', 'FontSize', 10);
end
title('Relative Area Between Response and Ground Truth Regions', 'FontSize', 14);
xlabel('Estimated Completeness Index','FontSize', 12);
ylabel('Relative Area','FontSize', 12);
hold off

% Save image
if ~strcmp(output_directory, '');
    filepath = [output_directory, '/region_abs_area.png'];
else
    filepath = 'region_abs_area.png';
end
print(hfig, filepath, '-dpng');

% Cleanup
delete(hfig);

% Missing answer histogram

% Create new figure
hfig = figure('Units', 'pixels', 'Position', [100, 100, 700, 500], 'PaperPositionMode', 'auto', 'InvertHardcopy', 'off', 'Color', [1, 1, 1]);

% Get vectors for the ECI
ECI_cr_0_5 = select_summary_subset(s.missing_regions, 'ECI', 'contrast_ratios', 0.5);
ECI_cr_1_5 = select_summary_subset(s.missing_regions, 'ECI', 'contrast_ratios', 1.5);
ECI_cr_4 = select_summary_subset(s.missing_regions, 'ECI', 'contrast_ratios', 4);

% Calcaulate histogram
max_ECI = 0.3;
histogram_edges = linspace(0, max_ECI, 11)';
missing_cr_0_5 = histc(ECI_cr_0_5, histogram_edges, 1);
if isempty(missing_cr_0_5); missing_cr_0_5 = zeros(size(histogram_edges)); end;
missing_cr_1_5 = histc(ECI_cr_1_5, histogram_edges, 1);
if isempty(missing_cr_1_5); missing_cr_1_5 = zeros(size(histogram_edges)); end;
missing_cr_4 = histc(ECI_cr_4, histogram_edges, 1);
if isempty(missing_cr_4); missing_cr_4 = zeros(size(histogram_edges)); end;


% Plot histogram for each contrast ratio
haxes = axes('Units', 'pixels', 'Position', [70, 50, 600, 400]);
    hold on
    hbar = bar(haxes, histogram_edges, [missing_cr_0_5, missing_cr_1_5, missing_cr_4], 'EdgeColor', 'k', 'BarLayout', 'grouped');
    set(hbar(1), 'FaceColor', [0, 0, 0]);
    set(hbar(2), 'FaceColor', [.4, .4, .4]);
    set(hbar(3), 'FaceColor', [.8, .8, .8]);
    set(haxes, 'Xlim', [-0.05, max_ECI + 0.05]);legend_labels = {'1:1.5', '1:2.5', '1:5'};
hlegend = legend(legend_labels, 'Location', 'best');
if verLessThan('matlab','8.4') % Stupid MATLAB changed lots of stuff about graphics
    hlt = get(hlegend, 'title');
    set(hlt, 'String', 'Contrast Ratio');
else
    hlt = text('Parent', hlegend.DecorationContainer, 'String', 'Contrast Ratio', 'HorizontalAlignment', 'center', ...
        'VerticalAlignment', 'bottom', 'Position', [0.5, 1.05], 'Units', 'normalized', 'FontSize', 10);
end
title('Missing Responses Histogram of Regions', 'FontSize', 14);
xlabel('Estimated Completeness Index','FontSize', 12);
ylabel('Number Missing Responses','FontSize', 12);
hold off

% Save image
if ~strcmp(output_directory, '');
    filepath = [output_directory, '/region_missing.png'];
else
    filepath = 'region_missing.png';
end
print(hfig, filepath, '-dpng');

% Cleanup
delete(hfig);

% ---- Dots -----

% Rename for convenience
s = summary_struct.dots;

% Found accuracy images

% Repeat images for each dot size
dot_sizes = [20, 50, 100, 200, 500]; % nanometers
for dot_radius = dot_sizes


    % Get vectors for the ECI
    ECI_cr_1= select_summary_subset(s.found_points, 'ECI', {'contrast_ratios'; 'radius'}, {1; dot_radius});
    ECI_cr_4 = select_summary_subset(s.found_points, 'ECI', {'contrast_ratios'; 'radius'}, {4; dot_radius});
    ECI_cr_10 = select_summary_subset(s.found_points, 'ECI', {'contrast_ratios'; 'radius'}, {10; dot_radius});

    % Get distance vectors and calculate relative absolute area
    dist_cr_1 = select_summary_subset(s.found_points, 'distances', {'contrast_ratios'; 'radius'}, {1; dot_radius});
    dist_cr_4 = select_summary_subset(s.found_points, 'distances', {'contrast_ratios'; 'radius'}, {4; dot_radius});
    dist_cr_10 = select_summary_subset(s.found_points, 'distances', {'contrast_ratios'; 'radius'}, {10; dot_radius});
    
    % Create new figure
    hfig = figure('Units', 'pixels', 'Position', [100, 100, 700, 500], 'PaperPositionMode', 'auto', 'InvertHardcopy', 'off', 'Color', [1, 1, 1]);

    % Plot scatter plot of each contrast ratio
    haxes = axes('Units', 'pixels', 'Position', [70, 50, 600, 400]);
    scatter(haxes, ECI_cr_1, dist_cr_1, 'Marker', 'o', 'MarkerEdgeColor', 'k');
    hold on
    scatter(haxes, ECI_cr_4, dist_cr_4, 'Marker', '^', 'MarkerEdgeColor', 'k');
    scatter(haxes, ECI_cr_10, dist_cr_10, 'Marker', 's', 'MarkerEdgeColor', 'k');
    set(haxes, 'Xlim', [0, 1]);
    legend_labels = {'1:1.5', '1:2.5', '1:5'};
    hlegend = legend(legend_labels, 'Location', 'best');
    if verLessThan('matlab','8.4') % Stupid MATLAB changed lots of stuff about graphics
        hlt = get(hlegend, 'title');
        set(hlt, 'String', 'Contrast Ratio');
    else
        hlt = text('Parent', hlegend.DecorationContainer, 'String', 'Contrast Ratio', 'HorizontalAlignment', 'center', ...
            'VerticalAlignment', 'bottom', 'Position', [0.5, 1.05], 'Units', 'normalized', 'FontSize', 10);
    end
    title(['Distance Between Response and Ground Truth ', num2str(dot_radius), ' nm Dots'], 'FontSize', 14);
    xlabel('Estimated Completeness Index','FontSize', 12);
    ylabel('Distance (nm)','FontSize', 12);
    hold off

    % Save image
    name = ['dots_', num2str(dot_radius), 'nm_distance.png'];
    if ~strcmp(output_directory, '');
        filepath = [output_directory, '/', name];
    else
        filepath = name;
    end
    print(hfig, filepath, '-dpng');

    % Cleanup
    delete(hfig);

    % Missing answer histogram
    % Create new figure
    hfig = figure('Units', 'pixels', 'Position', [100, 100, 700, 500], 'PaperPositionMode', 'auto', 'InvertHardcopy', 'off', 'Color', [1, 1, 1]);

    % Get vectors for the ECI
    ECI_cr_1= select_summary_subset(s.missing_points, 'ECI', {'contrast_ratios'; 'radius'}, {1; dot_radius});
    ECI_cr_4 = select_summary_subset(s.missing_points, 'ECI', {'contrast_ratios'; 'radius'}, {4; dot_radius});
    ECI_cr_10 = select_summary_subset(s.missing_points, 'ECI', {'contrast_ratios'; 'radius'}, {10; dot_radius});

    % Calcaulate histogram
    max_ECI = 1;
    histogram_edges = linspace(0, max_ECI, 11)';
    missing_cr_1 = histc(ECI_cr_1, histogram_edges, 1);
    if isempty(missing_cr_1); missing_cr_1 = zeros(size(histogram_edges)); end;
    missing_cr_4 = histc(ECI_cr_4, histogram_edges, 1);
    if isempty(missing_cr_4); missing_cr_4 = zeros(size(histogram_edges)); end;
    missing_cr_10 = histc(ECI_cr_10, histogram_edges, 1);
    if isempty(missing_cr_10); missing_cr_10 = zeros(size(histogram_edges)); end;


    % Plot histogram for each contrast ratio
    haxes = axes('Units', 'pixels', 'Position', [70, 50, 600, 400]);
    hold on
    hbar = bar(haxes, histogram_edges, [missing_cr_1, missing_cr_4, missing_cr_10], 'EdgeColor', 'k', 'BarLayout', 'grouped');
    set(hbar(1), 'FaceColor', [0, 0, 0]);
    set(hbar(2), 'FaceColor', [.4, .4, .4]);
    set(hbar(3), 'FaceColor', [.8, .8, .8]);
    set(haxes, 'Xlim', [-0.05, max_ECI + 0.05]);
    legend_labels = {'1:1.5', '1:2.5', '1:5'};
    hlegend = legend(legend_labels, 'Location', 'best');
    if verLessThan('matlab','8.4') % Stupid MATLAB changed lots of stuff about graphics
        hlt = get(hlegend, 'title');
        set(hlt, 'String', 'Contrast Ratio');
    else
        hlt = text('Parent', hlegend.DecorationContainer, 'String', 'Contrast Ratio', 'HorizontalAlignment', 'center', ...
            'VerticalAlignment', 'bottom', 'Position', [0.5, 1.05], 'Units', 'normalized', 'FontSize', 10);
    end
    title(['Missing Responses Histogram of ', num2str(dot_radius), ' nm Dots'], 'FontSize', 14);
    xlabel('Estimated Completeness Index','FontSize', 12);
    ylabel('Number Missing Responses','FontSize', 12);
    hold off

    % Save image
    name = ['dots_', num2str(dot_radius), 'nm_missing.png'];
    if ~strcmp(output_directory, '');
        filepath = [output_directory, '/', name];
    else
        filepath = name;
    end
    print(hfig, filepath, '-dpng');

    % Cleanup
    delete(hfig);
    
    % Extra answer histogram
    % Create new figure
    hfig = figure('Units', 'pixels', 'Position', [100, 100, 700, 500], 'PaperPositionMode', 'auto', 'InvertHardcopy', 'off', 'Color', [1, 1, 1]);

    % Get vectors for the ECI
    ECI_cr_1= select_summary_subset(s.extra_points, 'ECI', {'contrast_ratios'; 'radius'}, {1; dot_radius});
    ECI_cr_4 = select_summary_subset(s.extra_points, 'ECI', {'contrast_ratios'; 'radius'}, {4; dot_radius});
    ECI_cr_10 = select_summary_subset(s.extra_points, 'ECI', {'contrast_ratios'; 'radius'}, {10; dot_radius});

    % Calcaulate histogram
    max_ECI = 1;
    histogram_edges = linspace(0, max_ECI, 11)';
    extra_cr_1 = histc(ECI_cr_1, histogram_edges, 1);
    if isempty(extra_cr_1); extra_cr_1 = zeros(size(histogram_edges)); end;
    extra_cr_4 = histc(ECI_cr_4, histogram_edges, 1);
    if isempty(extra_cr_4); extra_cr_4 = zeros(size(histogram_edges)); end;
    extra_cr_10 = histc(ECI_cr_10, histogram_edges, 1);
    if isempty(extra_cr_10); extra_cr_10 = zeros(size(histogram_edges)); end;


    % Plot histogram for each contrast ratio
    haxes = axes('Units', 'pixels', 'Position', [70, 50, 600, 400]);
    hold on
    hbar = bar(haxes, histogram_edges, [extra_cr_1, extra_cr_4, extra_cr_10], 'EdgeColor', 'k', 'BarLayout', 'grouped');
    set(hbar(1), 'FaceColor', [0, 0, 0]);
    set(hbar(2), 'FaceColor', [.4, .4, .4]);
    set(hbar(3), 'FaceColor', [.8, .8, .8]);
    set(haxes, 'Xlim', [-0.05, max_ECI + 0.05]);
    legend_labels = {'1:1.5', '1:2.5', '1:5'};
    hlegend = legend(legend_labels, 'Location', 'best');
    if verLessThan('matlab','8.4') % Stupid MATLAB changed lots of stuff about graphics
        hlt = get(hlegend, 'title');
        set(hlt, 'String', 'Contrast Ratio');
    else
        hlt = text('Parent', hlegend.DecorationContainer, 'String', 'Contrast Ratio', 'HorizontalAlignment', 'center', ...
            'VerticalAlignment', 'bottom', 'Position', [0.5, 1.05], 'Units', 'normalized', 'FontSize', 10);
    end
    title(['Extra Responses Histogram of ', num2str(dot_radius), ' nm Dots'], 'FontSize', 14);
    xlabel('Estimated Completeness Index','FontSize', 12);
    ylabel('Number Missing Responses','FontSize', 12);
    hold off

    % Save image
    name = ['dots_', num2str(dot_radius), 'nm_extra.png'];
    if ~strcmp(output_directory, '');
        filepath = [output_directory, '/', name];
    else
        filepath = name;
    end
    print(hfig, filepath, '-dpng');

    % Cleanup
    delete(hfig);

end % End looping through each dot radius

% ---- Actin -----

% Rename for convenience
s = summary_struct.actin;

% Found accuracy images

% Repeat images for each dot size
line_widths = [26, 9]; % nanometers
for width = line_widths;

    % Get vectors for the ECI
    ECI_cr_1= select_summary_subset(s.found_lines, 'ECI', {'contrast_ratios'; 'width'}, {1; width});
    ECI_cr_4 = select_summary_subset(s.found_lines, 'ECI', {'contrast_ratios'; 'width'}, {4; width});
    ECI_cr_10 = select_summary_subset(s.found_lines, 'ECI', {'contrast_ratios'; 'width'}, {10; width});

    % Get area vectors and calculate area/length
    area_cr_1 = select_summary_subset(s.found_lines, 'area', {'contrast_ratios'; 'width'}, {1; width});
    area_cr_4 = select_summary_subset(s.found_lines, 'area', {'contrast_ratios'; 'width'}, {4; width});
    area_cr_10 = select_summary_subset(s.found_lines, 'area', {'contrast_ratios'; 'width'}, {10; width});
    length_cr_1 = select_summary_subset(s.found_lines, 'true_length', {'contrast_ratios'; 'width'}, {1; width});
    length_cr_4 = select_summary_subset(s.found_lines, 'true_length', {'contrast_ratios'; 'width'}, {4; width});
    length_cr_10 = select_summary_subset(s.found_lines, 'true_length', {'contrast_ratios'; 'width'}, {10; width});
    norm_area_cr_1 = area_cr_1 ./ length_cr_1;
    norm_area_cr_4 = area_cr_4 ./ length_cr_4;
    norm_area_cr_10 = area_cr_10 ./ length_cr_10;
    
    % Create new figure
    hfig = figure('Units', 'pixels', 'Position', [100, 100, 700, 500], 'PaperPositionMode', 'auto', 'InvertHardcopy', 'off', 'Color', [1, 1, 1]);

    % Plot scatter plot of each contrast ratio
    haxes = axes('Units', 'pixels', 'Position', [70, 50, 600, 400]);
    scatter(haxes, ECI_cr_1, norm_area_cr_1, 'Marker', 'o', 'MarkerEdgeColor', 'k');
    hold on
    scatter(haxes, ECI_cr_4, norm_area_cr_4, 'Marker', '^', 'MarkerEdgeColor', 'k');
    scatter(haxes, ECI_cr_10, norm_area_cr_10, 'Marker', 's', 'MarkerEdgeColor', 'k');
    set(haxes, 'Xlim', [0, 1]);
    legend_labels = {'1:1.5', '1:2.5', '1:5'};
    hlegend = legend(legend_labels, 'Location', 'best');
    if verLessThan('matlab','8.4') % Stupid MATLAB changed lots of stuff about graphics
        hlt = get(hlegend, 'title');
        set(hlt, 'String', 'Contrast Ratio');
    else
        hlt = text('Parent', hlegend.DecorationContainer, 'String', 'Contrast Ratio', 'HorizontalAlignment', 'center', ...
            'VerticalAlignment', 'bottom', 'Position', [0.5, 1.05], 'Units', 'normalized', 'FontSize', 10);
    end
    title(['Area Between Response and Ground Truth Curves ', num2str(width), ' nm Wide'], 'FontSize', 14);
    xlabel('Estimated Completeness Index','FontSize', 12);
    ylabel('Area/Curve Length (nm^2 / nm)','FontSize', 12);
    hold off

    % Save image
    name = ['actin_', num2str(width), 'nm_area.png'];
    if ~strcmp(output_directory, '');
        filepath = [output_directory, '/', name];
    else
        filepath = name;
    end
    print(hfig, filepath, '-dpng');

    % Cleanup
    delete(hfig);

    % Missing answer histogram
    % Create new figure
    hfig = figure('Units', 'pixels', 'Position', [100, 100, 700, 500], 'PaperPositionMode', 'auto', 'InvertHardcopy', 'off', 'Color', [1, 1, 1]);

    % Get vectors for the ECI
    ECI_cr_1= select_summary_subset(s.missing_lines, 'ECI', {'contrast_ratios'; 'width'}, {1; width});
    ECI_cr_4 = select_summary_subset(s.missing_lines, 'ECI', {'contrast_ratios'; 'width'}, {4; width});
    ECI_cr_10 = select_summary_subset(s.missing_lines, 'ECI', {'contrast_ratios'; 'width'}, {10; width});

    % Calcaulate histogram
    max_ECI = 1;
    histogram_edges = linspace(0, max_ECI, 11)';
    missing_cr_1 = histc(ECI_cr_1, histogram_edges, 1);
    if isempty(missing_cr_1); missing_cr_1 = zeros(size(histogram_edges)); end;
    missing_cr_4 = histc(ECI_cr_4, histogram_edges, 1);
    if isempty(missing_cr_4); missing_cr_4 = zeros(size(histogram_edges)); end;
    missing_cr_10 = histc(ECI_cr_10, histogram_edges, 1);
    if isempty(missing_cr_10); missing_cr_10 = zeros(size(histogram_edges)); end;


    % Plot histogram for each contrast ratio
    haxes = axes('Units', 'pixels', 'Position', [70, 50, 600, 400]);
    hold on
    hbar = bar(haxes, histogram_edges, [missing_cr_1, missing_cr_4, missing_cr_10], 'EdgeColor', 'k', 'BarLayout', 'grouped');
    set(hbar(1), 'FaceColor', [0, 0, 0]);
    set(hbar(2), 'FaceColor', [.4, .4, .4]);
    set(hbar(3), 'FaceColor', [.8, .8, .8]);
    set(haxes, 'Xlim', [-0.05, max_ECI + 0.05]);
    legend_labels = {'1:1.5', '1:2.5', '1:5'};
    hlegend = legend(legend_labels, 'Location', 'best');
    if verLessThan('matlab','8.4') % Stupid MATLAB changed lots of stuff about graphics
        hlt = get(hlegend, 'title');
        set(hlt, 'String', 'Contrast Ratio');
    else
        hlt = text('Parent', hlegend.DecorationContainer, 'String', 'Contrast Ratio', 'HorizontalAlignment', 'center', ...
            'VerticalAlignment', 'bottom', 'Position', [0.5, 1.05], 'Units', 'normalized', 'FontSize', 10);
    end
    title(['Missing Responses Histogram of ', num2str(width), ' nm Curves'], 'FontSize', 14);
    xlabel('Estimated Completeness Index','FontSize', 12);
    ylabel('Number Missing Responses','FontSize', 12);
    hold off

    % Save image
    name = ['actin_', num2str(width), 'nm_missing.png'];
    if ~strcmp(output_directory, '');
        filepath = [output_directory, '/', name];
    else
        filepath = name;
    end
    print(hfig, filepath, '-dpng');

    % Cleanup
    delete(hfig);
    
    % Extra answer histogram
    % Create new figure
    hfig = figure('Units', 'pixels', 'Position', [100, 100, 700, 500], 'PaperPositionMode', 'auto', 'InvertHardcopy', 'off', 'Color', [1, 1, 1]);

    % Get vectors for the ECI
    ECI_cr_1= select_summary_subset(s.extra_lines, 'ECI', {'contrast_ratios'; 'width'}, {1; width});
    ECI_cr_4 = select_summary_subset(s.extra_lines, 'ECI', {'contrast_ratios'; 'width'}, {4; width});
    ECI_cr_10 = select_summary_subset(s.extra_lines, 'ECI', {'contrast_ratios'; 'width'}, {10; width});

    % Calcaulate histogram
    max_ECI = 1;
    histogram_edges = linspace(0, max_ECI, 11)';
    extra_cr_1 = histc(ECI_cr_1, histogram_edges, 1);
    if isempty(extra_cr_1); extra_cr_1 = zeros(size(histogram_edges)); end;
    extra_cr_4 = histc(ECI_cr_4, histogram_edges, 1);
    if isempty(extra_cr_4); extra_cr_4 = zeros(size(histogram_edges)); end;
    extra_cr_10 = histc(ECI_cr_10, histogram_edges, 1);
    if isempty(extra_cr_10); extra_cr_10 = zeros(size(histogram_edges)); end;


    % Plot histogram for each contrast ratio
    haxes = axes('Units', 'pixels', 'Position', [70, 50, 600, 400]);
    hold on
    hbar = bar(haxes, histogram_edges, [extra_cr_1, extra_cr_4, extra_cr_10], 'EdgeColor', 'k', 'BarLayout', 'grouped');
    set(hbar(1), 'FaceColor', [0, 0, 0]);
    set(hbar(2), 'FaceColor', [.4, .4, .4]);
    set(hbar(3), 'FaceColor', [.8, .8, .8]);
    set(haxes, 'Xlim', [-0.05, max_ECI + 0.05]);
    legend_labels = {'1:1.5', '1:2.5', '1:5'};
    hlegend = legend(legend_labels, 'Location', 'best');
    if verLessThan('matlab','8.4') % Stupid MATLAB changed lots of stuff about graphics
        hlt = get(hlegend, 'title');
        set(hlt, 'String', 'Contrast Ratio');
    else
        hlt = text('Parent', hlegend.DecorationContainer, 'String', 'Contrast Ratio', 'HorizontalAlignment', 'center', ...
            'VerticalAlignment', 'bottom', 'Position', [0.5, 1.05], 'Units', 'normalized', 'FontSize', 10);
    end
    title(['Extra Responses Histogram of ', num2str(width), ' nm Curves'], 'FontSize', 14);
    xlabel('Estimated Completeness Index','FontSize', 12);
    ylabel('Number Missing Responses','FontSize', 12);
    hold off

    % Save image
    name = ['actin_', num2str(width), 'nm_extra.png'];
    if ~strcmp(output_directory, '');
        filepath = [output_directory, '/', name];
    else
        filepath = name;
    end
    print(hfig, filepath, '-dpng');

    % Cleanup
    delete(hfig);

end % End looping through each line width

% ---- Borders -----

% Rename for convenience
s = summary_struct.borders;

% Accuracy images
% Get vectors for the ECI
ECI_cr_0_5 = select_summary_subset(s.all_borders, 'ECI', 'contrast_ratios', 0.5);
ECI_cr_1_5 = select_summary_subset(s.all_borders, 'ECI', 'contrast_ratios', 1.5);
ECI_cr_4 = select_summary_subset(s.all_borders, 'ECI', 'contrast_ratios', 4);

% Get area vectors and calculate absolute area in square micrometers
abs_area_cr_0_5 = select_summary_subset(s.all_borders, 'abs_area', 'contrast_ratios', 0.5) / 1e6; %square micrometers
abs_area_cr_1_5 = select_summary_subset(s.all_borders, 'abs_area', 'contrast_ratios', 1.5) / 1e6;
abs_area_cr_4 = select_summary_subset(s.all_borders, 'abs_area', 'contrast_ratios', 4) / 1e6;

% Create new figure
hfig = figure('Units', 'pixels', 'Position', [100, 100, 700, 500], 'PaperPositionMode', 'auto', 'InvertHardcopy', 'off', 'Color', [1, 1, 1]);

% Plot scatter plot of each contrast ratio
haxes = axes('Units', 'pixels', 'Position', [70, 50, 600, 400]);
scatter(haxes, ECI_cr_0_5, abs_area_cr_0_5, 'Marker', 'o', 'MarkerEdgeColor', 'k');
hold on
scatter(haxes, ECI_cr_1_5, abs_area_cr_1_5, 'Marker', '^', 'MarkerEdgeColor', 'k');
scatter(haxes, ECI_cr_4, abs_area_cr_4, 'Marker', 's', 'MarkerEdgeColor', 'k');
set(haxes, 'Xlim', [0, 1]);
legend_labels = {'1:1.5', '1:2.5', '1:5'};
hlegend = legend(legend_labels, 'Location', 'best');
if verLessThan('matlab','8.4') % Stupid MATLAB changed lots of stuff about graphics
    hlt = get(hlegend, 'title');
    set(hlt, 'String', 'Contrast Ratio');
else
    hlt = text('Parent', hlegend.DecorationContainer, 'String', 'Contrast Ratio', 'HorizontalAlignment', 'center', ...
        'VerticalAlignment', 'bottom', 'Position', [0.5, 1.05], 'Units', 'normalized', 'FontSize', 10);
end
title('Absolute Area Between Response and Ground Truth Border', 'FontSize', 14);
xlabel('Estimated Completeness Index','FontSize', 12);
ylabel('Absolute Area (\mum^2)','FontSize', 12);
hold off

% Save image
if ~strcmp(output_directory, '');
    filepath = [output_directory, '/border_abs_area.png'];
else
    filepath = 'border_abs_area.png';
end
print(hfig, filepath, '-dpng');

% Cleanup
delete(hfig);

end

function [new_struct] = add_summaries(s1, s2, cat_dimension)
% Add together summary structures and returns a larger structure
%
% Input:
%   s1/2: summary structures with identical fields that contain matricies 
%       to be stacked together.
%   cat_dimension: Dimension along which to stack the results. Optional, 
%       default = 1.
% Output:
%   new_struct: structures with the same fields and stacked entries.

% Preallocate a structure.
new_struct = struct();

% Get first level fieldnames and loop through each name (region, dots, etc.)
fn_lvl1 = fieldnames(s1);
for field_index_lvl1 = 1:size(fn_lvl1, 1)
    lvl1_name = fn_lvl1{field_index_lvl1};

    % Get second level fieldnames and loop through each name (found_region, missing_region, etc.)
    fn_lvl2 = fieldnames(s1.(lvl1_name));
    for field_index_lvl2 = 1:size(fn_lvl2, 1)
        lvl2_name = fn_lvl2{field_index_lvl2};

        % Add together the two structures with matrices/cell arrays
        new_substruct = add_struct_fields_together(s1.(lvl1_name).(lvl2_name), s2.(lvl1_name).(lvl2_name), cat_dimension);
        
        % Put the new substructure in the new structure
        new_struct.(lvl1_name).(lvl2_name) = new_substruct;
    end
end
end

function [values] = select_summary_subset(summary_substructure, value_fieldnames, selector_fieldnames, selector_values)
% Get the specified values in a summary structure that match the selector
%   criteria.
%
% Input:
%   summary_substructure: substructure from a summary, i.e. s.region.found_regions
%   value_fieldnames: String or cell array of strings of fieldnames for the
%       desired values to be returned. 
%   selector_fieldnames: String or cell array of strings of fieldnames for 
%       the desired values to be returned.
%   selector_valeus: String/numeric value or cell array of strings/numeric 
%       values  of fieldnames for the desired values to be returned.
% Output:
%   values: matrix or cell array of matricies of the same shape as
%       value_fieldnames. 

% Default behavior
remove_values_cell_array = false;

% Rename for convenience
s = summary_substructure;

% If any of the arguments are single values, put them in a cell array.
if ~iscell(value_fieldnames)
    value_fieldnames = {value_fieldnames};
    remove_values_cell_array = true;
end
if ~iscell(selector_fieldnames)
    selector_fieldnames = {selector_fieldnames};
end
if ~iscell(selector_values)
    selector_values = {selector_values};
end

% Get size of the values
number_entries = size(s.(value_fieldnames{1}), 1);

% Get binary selection vectors for each selector fieldname
binary_selection_vector = true(number_entries, 1);
for selector_index = 1:length(selector_fieldnames(:))
    
    % For string values
    if iscell(s.(selector_fieldnames{selector_index}))
        new_binary_vector = strcmp(s.(selector_fieldnames{selector_index}), selector_values{selector_index});
    
    % For double values
    elseif isa(s.(selector_fieldnames{selector_index}), 'double')
        new_binary_vector = s.(selector_fieldnames{selector_index}) == selector_values{selector_index};
    end

    % Get intersection of existing and new binary vectors
    binary_selection_vector = binary_selection_vector & new_binary_vector;
end

% Create output cell array
values = cell(size(value_fieldnames));

% Put values into the value cell array
for value_index = 1:length(value_fieldnames(:))
    values{value_index} = s.(value_fieldnames{value_index})(binary_selection_vector);
end

% If there's only one value field, just return the values instead of a cell array
if remove_values_cell_array
    values = values{1};
end
end
