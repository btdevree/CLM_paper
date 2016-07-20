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
        summary_struct = add_summaries(summary_struct, test_summary);
    end
    
% No files, give error
else
    error('No test_summary files found');
end

disp('Hi')

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

% preallocate a structure.
new_struct = struct('region', [], 'dots', [], 'actin', [], 'borders', []);
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