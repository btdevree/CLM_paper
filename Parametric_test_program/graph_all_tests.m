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
legend_labels = {'1:1.5'; '1:2.5'; '1:5'};
hlegend = legend(legend_labels, 'Location', 'best');
if verLessThan('matlab','8.4') % Stupid MATLAB changed lots of stuff about graphics
    hlt = get(hlegend, 'title');
    set(hlt, 'String', 'Contrast Ratio');
else
    hlt = text('Parent', hlegend.DecorationContainer, 'String', 'Contrast Ratio', 'HorizontalAlignment', 'center', ...
        'VerticalAlignment', 'bottom', 'Position', [0.5, 1.05, 0], 'Units', 'normalized', 'FontSize', 10);
end
title('Relative Area Between Response and Ground Truth Regions', 'FontSize', 16);
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
histogram_edges = linspace(0, max_ECI, 11);
missing_cr_0_5 = histc(ECI_cr_0_5, histogram_edges);
missing_cr_1_5 = histc(ECI_cr_1_5, histogram_edges);
missing_cr_4 = histc(ECI_cr_4, histogram_edges);

% Plot histogram for each contrast ratio
haxes = axes('Units', 'pixels', 'Position', [70, 50, 600, 400]);
bar(haxes, histogram_edges, missing_cr_0_5, 'FaceColor', [0, 0, 0], 'EdgeColor', 'k');
hold on
bar(haxes, histogram_edges, missing_cr_1_5, 'FaceColor', [.4, .4, .4], 'EdgeColor', 'k');
bar(haxes, histogram_edges, missing_cr_4, 'FaceColor', [.8, .8, .8], 'EdgeColor', 'k');
set(haxes, 'Xlim', [0, max_ECI]);
legend_labels = {'1:1.5'; '1:2.5'; '1:5'};
hlegend = legend(legend_labels, 'Location', 'best');
if verLessThan('matlab','8.4') % Stupid MATLAB changed lots of stuff about graphics
    hlt = get(hlegend, 'title');
    set(hlt, 'String', 'Contrast Ratio');
else
    hlt = text('Parent', hlegend.DecorationContainer, 'String', 'Contrast Ratio', 'HorizontalAlignment', 'center', ...
        'VerticalAlignment', 'bottom', 'Position', [0.5, 1.05, 0], 'Units', 'normalized', 'FontSize', 10);
end
title('Missing Responses Histogram of Regions', 'FontSize', 16);
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
    value_fieldnames ={value_fieldnames};
    remove_values_cell_array = true;
end
if ~iscell(selector_fieldnames)
    selector_fieldnames ={selector_fieldnames};
end
if ~iscell(selector_values)
    selector_values ={selector_values};
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
    values{value_index} = s.(value_fieldnames{selector_index})(binary_selection_vector);
end

% If there's only one value field, just return the values instead of a cell array
if remove_values_cell_array
    values = values{1};
end
end
