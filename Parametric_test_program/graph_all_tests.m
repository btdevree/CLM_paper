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
% Define CR values and legend
cr_values = [0.5, 1, 4];
legend_labels = {'1:1.5', '1:2', '1:5'};

% Get vectors for the ECI
ECI_cr1 = select_summary_subset(s.found_regions, 'ECI', 'contrast_ratios', cr_values(1));
ECI_cr2 = select_summary_subset(s.found_regions, 'ECI', 'contrast_ratios', cr_values(2));
ECI_cr3 = select_summary_subset(s.found_regions, 'ECI', 'contrast_ratios', cr_values(3));

% Get area vectors and calculate relative absolute area
abs_area_cr1 = select_summary_subset(s.found_regions, 'abs_area', 'contrast_ratios', cr_values(1));
abs_area_cr2 = select_summary_subset(s.found_regions, 'abs_area', 'contrast_ratios', cr_values(2));
abs_area_cr3 = select_summary_subset(s.found_regions, 'abs_area', 'contrast_ratios', cr_values(3));
true_area_cr1 = select_summary_subset(s.found_regions, 'true_area', 'contrast_ratios', cr_values(1));
true_area_cr2 = select_summary_subset(s.found_regions, 'true_area', 'contrast_ratios', cr_values(2));
true_area_cr3 = select_summary_subset(s.found_regions, 'true_area', 'contrast_ratios', cr_values(3));
relative_area_cr1 = abs_area_cr1 ./ true_area_cr1;
relative_area_cr2 = abs_area_cr2 ./ true_area_cr2;
relative_area_cr3 = abs_area_cr3 ./ true_area_cr3;

% Create new figure
hfig = figure('Units', 'pixels', 'Position', [100, 100, 600, 375], 'PaperPositionMode', 'auto', 'InvertHardcopy', 'off', 'Color', [1, 1, 1]);

% Plot scatter plot of each contrast ratio
haxes = axes('Units', 'pixels', 'Position', [70, 50, 500, 275]);
scatter(haxes, ECI_cr1, relative_area_cr1, 'Marker', 'o', 'MarkerEdgeColor', 'k', 'MarkerFaceColor', [1,1,1]);
hold on
scatter(haxes, ECI_cr2, relative_area_cr2, 'Marker', 'o', 'MarkerEdgeColor','k', 'MarkerFaceColor', [0.7,0.7,0.7]);
scatter(haxes, ECI_cr3, relative_area_cr3, 'Marker', 'o', 'MarkerEdgeColor', 'k', 'MarkerFaceColor', [0.4,0.4,0.4]);
set(haxes, 'Xlim', [0, 0.5]);
hlegend = legend(legend_labels, 'Position', [450, 250, 90, 40]);
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
refresh
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
ECI_cr1 = select_summary_subset(s.missing_regions, 'ECI', 'contrast_ratios', cr_values(1));
ECI_cr2 = select_summary_subset(s.missing_regions, 'ECI', 'contrast_ratios', cr_values(2));
ECI_cr3 = select_summary_subset(s.missing_regions, 'ECI', 'contrast_ratios', cr_values(3));

% Calcaulate histogram
max_ECI = 0.3;
histogram_edges = linspace(0, max_ECI, 11)';
missing_cr1 = histc(ECI_cr1, histogram_edges, 1);
if isempty(missing_cr1); missing_cr1 = zeros(size(histogram_edges)); end;
missing_cr2 = histc(ECI_cr2, histogram_edges, 1);
if isempty(missing_cr2); missing_cr2 = zeros(size(histogram_edges)); end;
missing_cr3 = histc(ECI_cr3, histogram_edges, 1);
if isempty(missing_cr3); missing_cr3 = zeros(size(histogram_edges)); end;


% Plot histogram for each contrast ratio
haxes = axes('Units', 'pixels', 'Position', [70, 50, 600, 400]);
    hold on
    hbar = bar(haxes, histogram_edges, [missing_cr1, missing_cr2, missing_cr3], 'EdgeColor', 'k', 'BarLayout', 'grouped');
    set(hbar(1), 'FaceColor', [1, 1, 1]);
    set(hbar(2), 'FaceColor', [.7, .7, .7]);
    set(hbar(3), 'FaceColor', [.4, .4, .4]);
    set(haxes, 'Xlim', [-0.05, max_ECI + 0.05]);
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
refresh
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
dot_sizes = [50; 100; 200]; % nanometers
cr_values = [4, 19; 2, 9; 1, 4];
legend_labels_distance = {{'1:5', '1:20'}, {'1:3', '1:10'}, {'1:2', '1:5'}};
legend_labels_errors = {{'1:5 extra', '1:5 missing', '1:20 extra', '1:20 missing'},...
                          {'1:3 extra', '1:3 missing', '1:10 extra', '1:10 missing'},...
                          {'1:2 extra', '1:2 missing', '1:5 extra', '1:5 missing'}};
for size_index = 1:size(dot_sizes, 1)
    dot_radius = dot_sizes(size_index);
    
    % Get vectors for the ECI
    ECI_cr1= select_summary_subset(s.found_points, 'ECI', {'contrast_ratios'; 'radius'}, {cr_values(size_index, 1); dot_radius});
    ECI_cr2 = select_summary_subset(s.found_points, 'ECI', {'contrast_ratios'; 'radius'}, {cr_values(size_index, 2); dot_radius});
    
    % Get distance vectors and calculate relative absolute area
    dist_cr1 = select_summary_subset(s.found_points, 'distances', {'contrast_ratios'; 'radius'}, {cr_values(size_index, 1); dot_radius});
    dist_cr2 = select_summary_subset(s.found_points, 'distances', {'contrast_ratios'; 'radius'}, {cr_values(size_index, 2); dot_radius});
    
    % Create new figure
    hfig = figure('Units', 'pixels', 'Position', [100, 100, 700, 450], 'PaperPositionMode', 'auto', 'InvertHardcopy', 'off', 'Color', [1, 1, 1]);

    % Plot scatter plot of each contrast ratio
    haxes = axes('Units', 'pixels', 'Position', [70, 50, 600, 350]);
    scatter(haxes, ECI_cr1, dist_cr1, 'Marker', 'o', 'MarkerEdgeColor', 'k', 'MarkerFaceColor', [1, 1, 1]);
    hold on
    scatter(haxes, ECI_cr2, dist_cr2, 'Marker', 'o', 'MarkerEdgeColor', 'k', 'MarkerFaceColor', [.5, .5, .5]);
    set(haxes, 'Xlim', [0, 1]);
    hlegend = legend(legend_labels_distance{size_index}, 'Position', [580, 350, 90, 30]);
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
    refresh
    name = ['dots_', num2str(dot_radius), 'nm_distance.png'];
    if ~strcmp(output_directory, '');
        filepath = [output_directory, '/', name];
    else
        filepath = name;
    end
    print(hfig, filepath, '-dpng');

    % Cleanup
    delete(hfig);

    % Missing and extra answer histogram
    % Create new figure
    hfig = figure('Units', 'pixels', 'Position', [100, 100, 500, 500], 'PaperPositionMode', 'auto', 'InvertHardcopy', 'off', 'Color', [1, 1, 1]);

    % Get vectors for the ECI
    ECI_cr1_extra = select_summary_subset(s.extra_points, 'ECI', {'contrast_ratios'; 'radius'}, {cr_values(size_index, 1); dot_radius});
    ECI_cr2_extra = select_summary_subset(s.extra_points, 'ECI', {'contrast_ratios'; 'radius'}, {cr_values(size_index, 2); dot_radius});
    ECI_cr1_missing = select_summary_subset(s.missing_points, 'ECI', {'contrast_ratios'; 'radius'}, {cr_values(size_index, 1); dot_radius});
    ECI_cr2_missing = select_summary_subset(s.missing_points, 'ECI', {'contrast_ratios'; 'radius'}, {cr_values(size_index, 2); dot_radius});     
    ECI_cr1_found = select_summary_subset(s.found_points, 'ECI', {'contrast_ratios'; 'radius'}, {cr_values(size_index, 1); dot_radius});
    ECI_cr2_found = select_summary_subset(s.found_points, 'ECI', {'contrast_ratios'; 'radius'}, {cr_values(size_index, 2); dot_radius});
    
    % Calcaulate histogram
    max_ECI = 1;
    histogram_edges = linspace(0, max_ECI, 11)';
    extra_cr1 = histc(ECI_cr1_extra, histogram_edges, 1);
    if isempty(extra_cr1); extra_cr1 = zeros(size(histogram_edges)); end;
    extra_cr2 = histc(ECI_cr2_extra, histogram_edges, 1);
    if isempty(extra_cr2); extra_cr2 = zeros(size(histogram_edges)); end;
    missing_cr1 = histc(ECI_cr1_missing, histogram_edges, 1);
    if isempty(missing_cr1); missing_cr1 = zeros(size(histogram_edges)); end;
    missing_cr2 = histc(ECI_cr2_missing, histogram_edges, 1);
    if isempty(missing_cr2); missing_cr2 = zeros(size(histogram_edges)); end;
    found_cr1 = histc(ECI_cr1_found, histogram_edges, 1);
    if isempty(found_cr1); found_cr1 = zeros(size(histogram_edges)); end;
    found_cr2 = histc(ECI_cr2_found, histogram_edges, 1);
    if isempty(found_cr2); found_cr2 = zeros(size(histogram_edges)); end;
    
    % Normalize histogram counts
    extra_cr1_normalized = extra_cr1 ./ (missing_cr1 + found_cr1);
    extra_cr2_normalized = extra_cr2 ./ (missing_cr2 + found_cr2);
    missing_cr1_normalized = missing_cr1 ./ (missing_cr1 + found_cr1);
    missing_cr2_normalized = missing_cr2 ./ (missing_cr2 + found_cr2);
    
    % Plot histogram for each contrast ratio
    haxes = axes('Units', 'pixels', 'Position', [70, 50, 400, 400], 'LineWidth', 1);
    hold on
    hbar = bar(haxes, histogram_edges, [extra_cr1_normalized, missing_cr1_normalized, extra_cr2_normalized, missing_cr2_normalized],...
        'EdgeColor', 'k', 'BarLayout', 'grouped', 'BarWidth', .8, 'LineWidth', 1);
    set(hbar(1), 'FaceColor', [.6, .6, .6]);
    set(hbar(2), 'FaceColor', 'none'); % Add .5 grey hatches
    set(hbar(3), 'FaceColor', [0, 0, 0]);
    set(hbar(4), 'FaceColor','none'); % Add black hatches   
    set(haxes, 'Xlim', [-0.05, max_ECI + 0.05]);
    set(haxes, 'Ylim', [0, 1.05]);
    set(haxes, 'Xtick', [0:.1:1]);
    hlegend = legend(legend_labels_errors{size_index}, 'Position', [320, 330, 150, 80]);
    
    % Add hatchs and legend title - must happen after legend is made
    if verLessThan('matlab','8.4') % Stupid MATLAB changed lots of stuff about graphics
        hlt = get(hlegend, 'title');
        set(hlt, 'String', 'Contrast Ratio');
        hatchfill2(hbar(2), 'single', 'HatchLineWidth', 2, 'HatchColor', [.6, .6, .6], 'HatchAngle', 0);
        hatchfill2(hbar(4), 'single', 'HatchLineWidth', 2, 'HatchColor', [0, 0, 0], 'HatchAngle', 0);
        hlegendpatch = findobj(hlegend, 'type', 'patch');
        hatchfill2(hlegendpatch(1), 'single', 'HatchLineWidth',2, 'HatchColor', [0, 0, 0], 'HatchAngle', 0);
        hatchfill2(hlegendpatch(3), 'single', 'HatchLineWidth', 2, 'HatchColor', [.6, .6, .6], 'HatchAngle', 0);
        hhatch = findobj(gcf,'Type','line');
        for hatch_index = 1:length(hhatch)
            uistack(hhatch(hatch_index), 'bottom'); % Must be done one-at-a-time
        end
    else
        hlt = text('Parent', hlegend.DecorationContainer, 'String', 'Contrast Ratio', 'HorizontalAlignment', 'center', ...
            'VerticalAlignment', 'bottom', 'Position', [0.5, 1.05], 'Units', 'normalized', 'FontSize', 10);
    end
    title(['Extra and Missing Responses for ', num2str(dot_radius), ' nm Dots'], 'FontSize', 14);
    xlabel('Estimated Completeness Index','FontSize', 12);
    ylabel('Fraction of Ground Truth Dots','FontSize', 12);
    hold off

    % Save image
    refresh
    name = ['dots_', num2str(dot_radius), 'nm_detection_errors.png'];
    if ~strcmp(output_directory, '');
        filepath = [output_directory, '/', name];
    else
        filepath = name;
    end
    print(hfig, filepath, '-dpng');

    % Cleanup
    T = evalc('close(hfig)'); % Something is wrong with the hatching listeners, but we'll supress output and ignore it because the saved images are fine.

end % End looping through each dot radius

% ---- Actin -----

% Rename for convenience
s = summary_struct.actin;

% Found accuracy images

% Repeat images for each dot size
line_widths = [9; 26]; % nanometers
cr_values = [5, 19; 3, 9];
legend_labels_distance = {{'1:6', '1:20'}, {'1:4', '1:10'}};
legend_labels_errors = {{'1:6 extra', '1:6 missing', '1:20 extra', '1:20 missing'},...
                          {'1:4 extra', '1:4 missing', '1:10 extra', '1:10 missing'}};
for width_index = 1:size(line_widths, 1);
    width = line_widths(width_index);

    % Get vectors for the ECI
    ECI_cr1= select_summary_subset(s.found_lines, 'ECI', {'contrast_ratios'; 'width'}, {cr_values(width_index, 1); width});
    ECI_cr2 = select_summary_subset(s.found_lines, 'ECI', {'contrast_ratios'; 'width'}, {cr_values(width_index, 2); width});

    % Get area vectors and calculate area/length
    area_cr1 = select_summary_subset(s.found_lines, 'area', {'contrast_ratios'; 'width'}, {cr_values(width_index, 1); width});
    area_cr2 = select_summary_subset(s.found_lines, 'area', {'contrast_ratios'; 'width'}, {cr_values(width_index, 2); width});
    length_cr1 = select_summary_subset(s.found_lines, 'true_length', {'contrast_ratios'; 'width'}, {cr_values(width_index, 1); width});
    length_cr2 = select_summary_subset(s.found_lines, 'true_length', {'contrast_ratios'; 'width'}, {cr_values(width_index, 2); width});
    norm_area_cr1 = area_cr1 ./ length_cr1;
    norm_area_cr2 = area_cr2 ./ length_cr2;
    
    % Create new figure
    hfig = figure('Units', 'pixels', 'Position', [100, 100, 700, 500], 'PaperPositionMode', 'auto', 'InvertHardcopy', 'off', 'Color', [1, 1, 1]);

    % Plot scatter plot of each contrast ratio
    haxes = axes('Units', 'pixels', 'Position', [70, 50, 600, 400]);
    scatter(haxes, ECI_cr1, norm_area_cr1, 'Marker', 'o', 'MarkerEdgeColor', 'k', 'MarkerFaceColor', [1, 1, 1]);
    hold on
    scatter(haxes, ECI_cr2, norm_area_cr2, 'Marker', 'o', 'MarkerEdgeColor', 'k', 'MarkerFaceColor', [.5, .5, .5]);
    set(haxes, 'Xlim', [0, 1]);
    hlegend = legend(legend_labels_distance{size_index}, 'Position', [570, 370, 90, 30]);
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
    refresh
    name = ['actin_', num2str(width), 'nm_area.png'];
    if ~strcmp(output_directory, '');
        filepath = [output_directory, '/', name];
    else
        filepath = name;
    end
    print(hfig, filepath, '-dpng');

    % Cleanup
    delete(hfig);
    
    % Missing and extra answer histogram
    % Create new figure
    hfig = figure('Units', 'pixels', 'Position', [100, 100, 500, 500], 'PaperPositionMode', 'auto', 'InvertHardcopy', 'off', 'Color', [1, 1, 1]);

    % Get vectors for the ECI
    ECI_cr1_extra = select_summary_subset(s.extra_lines, 'ECI', {'contrast_ratios'; 'width'}, {cr_values(width_index, 1); width});
    ECI_cr2_extra = select_summary_subset(s.extra_lines, 'ECI', {'contrast_ratios'; 'width'}, {cr_values(width_index, 2); width});
    ECI_cr1_missing = select_summary_subset(s.missing_lines, 'ECI', {'contrast_ratios'; 'width'}, {cr_values(width_index, 1); width});
    ECI_cr2_missing = select_summary_subset(s.missing_lines, 'ECI', {'contrast_ratios'; 'width'}, {cr_values(width_index, 2); width});     
    ECI_cr1_found = select_summary_subset(s.found_lines, 'ECI', {'contrast_ratios'; 'width'}, {cr_values(width_index, 1); width});
    ECI_cr2_found = select_summary_subset(s.found_lines, 'ECI', {'contrast_ratios'; 'width'}, {cr_values(width_index, 2); width});
    
    % Calcaulate histogram
    max_ECI = 1;
    histogram_edges = linspace(0, max_ECI, 11)';
    extra_cr1 = histc(ECI_cr1_extra, histogram_edges, 1);
    if isempty(extra_cr1); extra_cr1 = zeros(size(histogram_edges)); end;
    extra_cr2 = histc(ECI_cr2_extra, histogram_edges, 1);
    if isempty(extra_cr2); extra_cr2 = zeros(size(histogram_edges)); end;
    missing_cr1 = histc(ECI_cr1_missing, histogram_edges, 1);
    if isempty(missing_cr1); missing_cr1 = zeros(size(histogram_edges)); end;
    missing_cr2 = histc(ECI_cr2_missing, histogram_edges, 1);
    if isempty(missing_cr2); missing_cr2 = zeros(size(histogram_edges)); end;
    found_cr1 = histc(ECI_cr1_found, histogram_edges, 1);
    if isempty(found_cr1); found_cr1 = zeros(size(histogram_edges)); end;
    found_cr2 = histc(ECI_cr2_found, histogram_edges, 1);
    if isempty(found_cr2); found_cr2 = zeros(size(histogram_edges)); end;
    
    % Normalize histogram counts
    extra_cr1_normalized = extra_cr1 ./ (missing_cr1 + found_cr1);
    extra_cr2_normalized = extra_cr2 ./ (missing_cr2 + found_cr2);
    missing_cr1_normalized = missing_cr1 ./ (missing_cr1 + found_cr1);
    missing_cr2_normalized = missing_cr2 ./ (missing_cr2 + found_cr2);
    
    % Plot histogram for each contrast ratio
    haxes = axes('Units', 'pixels', 'Position', [70, 50, 400, 400], 'LineWidth', 1);
    hold on
    hbar = bar(haxes, histogram_edges, [extra_cr1_normalized, missing_cr1_normalized, extra_cr2_normalized, missing_cr2_normalized],...
        'EdgeColor', 'k', 'BarLayout', 'grouped', 'BarWidth', .8, 'LineWidth', 1);
    set(hbar(1), 'FaceColor', [.6, .6, .6]);
    set(hbar(2), 'FaceColor', 'none'); % Add .5 grey hatches
    set(hbar(3), 'FaceColor', [0, 0, 0]);
    set(hbar(4), 'FaceColor','none'); % Add black hatches   
    set(haxes, 'Xlim', [-0.05, max_ECI + 0.05]);
    set(haxes, 'Ylim', [0, 1.05]);
    set(haxes, 'XTick', [0:.1:1]);
    hlegend = legend(legend_labels_errors{width_index}, 'Position', [300, 330, 150, 80]);
    
    % Add hatchs and legend title - must happen after legend is made
    if verLessThan('matlab','8.4') % Stupid MATLAB changed lots of stuff about graphics
        hlt = get(hlegend, 'title');
        set(hlt, 'String', 'Contrast Ratio');
        hatchfill2(hbar(2), 'single', 'HatchLineWidth', 2, 'HatchColor', [.6, .6, .6], 'HatchAngle', 0);
        hatchfill2(hbar(4), 'single', 'HatchLineWidth', 2, 'HatchColor', [0, 0, 0], 'HatchAngle', 0);
        hlegendpatch = findobj(hlegend, 'type', 'patch');
        hatchfill2(hlegendpatch(1), 'single', 'HatchLineWidth',2, 'HatchColor', [0, 0, 0], 'HatchAngle', 0);
        hatchfill2(hlegendpatch(3), 'single', 'HatchLineWidth', 2, 'HatchColor', [.6, .6, .6], 'HatchAngle', 0);
        hhatch = findobj(gcf,'Type','line');
        for hatch_index = 1:length(hhatch)
            uistack(hhatch(hatch_index), 'bottom'); % Must be done one-at-a-time
        end
    else
        hlt = text('Parent', hlegend.DecorationContainer, 'String', 'Contrast Ratio', 'HorizontalAlignment', 'center', ...
            'VerticalAlignment', 'bottom', 'Position', [0.5, 1.05], 'Units', 'normalized', 'FontSize', 10);
    end
    title(['Extra and Missing Responses for ', num2str(dot_radius), ' nm Wide Curves'], 'FontSize', 14);
    xlabel('Estimated Completeness Index','FontSize', 12);
    ylabel('Fraction of Ground Truth Lines','FontSize', 12);
    hold off

    % Save image
    refresh
    name = ['actin_', num2str(dot_radius), 'nm_detection_errors.png'];
    if ~strcmp(output_directory, '');
        filepath = [output_directory, '/', name];
    else
        filepath = name;
    end
    print(hfig, filepath, '-dpng');

    % Cleanup
    T = evalc('close(hfig)'); % Something is wrong with the hatching listeners, but we'll supress output and ignore it because the saved images are fine.
    
end % End looping through each line width

% ---- Borders -----

% Rename for convenience
s = summary_struct.borders;

% Roughness classes
roughness_values = [.45; .6; .75];

% Get average ECI values for each event number group
ECI1_vector = select_summary_subset(s.all_borders, 'ECI', 'roughness', roughness_values(1));
ECI2_vector = select_summary_subset(s.all_borders, 'ECI', 'roughness', roughness_values(2));
ECI3_vector = select_summary_subset(s.all_borders, 'ECI', 'roughness', roughness_values(3));

% Get vectors for the fractal dim
FD1_vector = select_summary_subset(s.all_borders, 'fractal_dim', 'roughness', roughness_values(1));
FD2_vector = select_summary_subset(s.all_borders, 'fractal_dim', 'roughness', roughness_values(2));
FD3_vector = select_summary_subset(s.all_borders, 'fractal_dim', 'roughness', roughness_values(3));
FD1_mean = mean(FD1_vector, 1);
FD2_mean = mean(FD2_vector, 1);
FD3_mean = mean(FD3_vector, 1);
FD1_stdev = std(FD1_vector, 0, 1);
FD2_stdev = std(FD2_vector, 0, 1);
FD3_stdev = std(FD3_vector, 0, 1);
legend_labels = {num2str(FD1_mean, 3), num2str(FD2_mean, 3), num2str(FD3_mean, 3)};


% Get area vectors and calculate absolute area in square micrometers
abs_area_ECI1 = select_summary_subset(s.all_borders, 'abs_area', 'roughness', roughness_values(1)) / 1e6; %square micrometers
abs_area_ECI2 = select_summary_subset(s.all_borders, 'abs_area', 'roughness', roughness_values(2)) / 1e6;
abs_area_ECI3 = select_summary_subset(s.all_borders, 'abs_area', 'roughness', roughness_values(3)) / 1e6;

% Create new figure
hfig = figure('Units', 'pixels', 'Position', [100, 100, 600, 375], 'PaperPositionMode', 'auto', 'InvertHardcopy', 'off', 'Color', [1, 1, 1]);

% Plot scatter plot of each contrast ratio
haxes = axes('Units', 'pixels', 'Position', [70, 50, 500, 275]);
scatter(haxes, ECI1_vector, abs_area_ECI1, 'Marker', 'o', 'MarkerEdgeColor', 'k', 'MarkerFaceColor', [1,1,1]);
hold on
scatter(haxes, ECI2_vector, abs_area_ECI2, 'Marker', 'o', 'MarkerEdgeColor','k', 'MarkerFaceColor', [0.7,0.7,0.7]);
scatter(haxes, ECI3_vector, abs_area_ECI3, 'Marker', 'o', 'MarkerEdgeColor', 'k', 'MarkerFaceColor', [0.4,0.4,0.4]);
set(haxes, 'Xlim', [0, 0.5]);
hlegend = legend(legend_labels, 'Position', [450, 260, 90, 35]);
if verLessThan('matlab','8.4') % Stupid MATLAB changed lots of stuff about graphics
    hlt = get(hlegend, 'title');
    set(hlt, 'String', 'Fractal Dimension of Border');
else
    hlt = text('Parent', hlegend.DecorationContainer, 'String', 'Fractal Dimension of Border', 'HorizontalAlignment', 'center', ...
        'VerticalAlignment', 'bottom', 'Position', [0.5, 1.05], 'Units', 'normalized', 'FontSize', 10);
end
title('Area Between Response and Ground Truth Region Border', 'FontSize', 14);
xlabel('Estimated Completeness Index','FontSize', 12);
ylabel('Area (micrometers squared)','FontSize', 12);
hold off

% Save image
refresh
if ~strcmp(output_directory, '');
    filepath = [output_directory, '/border_abs_area.png'];
else
    filepath = 'border_abs_area.png';
end
print(hfig, filepath, '-dpng');

% Save border group uncertantity info
fileID = fopen('border_group_data.txt','w');
fprintf(fileID,'Mean fractal dimension = %4.3f with standard deviation = %6.5f\n', FD1_mean, FD1_stdev);
fprintf(fileID,'Mean fractal dimension = %4.3f with standard deviation = %6.5f\n', FD2_mean, FD2_stdev);
fprintf(fileID,'Mean fractal dimension = %4.3f with standard deviation = %6.5f\n', FD3_mean, FD3_stdev);
fclose(fileID);

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
%   selector_values: String/numeric value or cell array of strings/numeric 
%       values of fieldnames for the desired values to be returned.
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
    binary_selection_vector = binary_selection_vector & new_binary_vector(:);
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
