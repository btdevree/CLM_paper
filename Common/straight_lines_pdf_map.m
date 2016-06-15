function [pdf_map, line_start_coords, line_end_coords] = straight_lines_pdf_map(parameter_struct, number_of_lines, line_width, line_to_background_ratio, line_min_length, line_max_length)
%STRAIGHT_LINES_PDF_MAP Makes a pdf map of line structures in the field of
% view.
%
% Each line is randomly placed within the image, assuming a Carteisain
% coordinate system. The origin is assumed to be in the lower left corner
% of the lower left pixel unless the bounds indicate otherwise in the 
% parameter structure.
%
% Inputs:
%   parameter_struct: parameter structure made with
%       test_movie_parameters_dv.
%   number_of_lines: integer, number of lines to put in the pdf map.
%   line_width: width of the lines, in nanometers.
%   line_to_background_ratio: ratio of the density of events inside the
%       line region to the general background. Use Inf for a image with no 
%       background.
%   line_min_length: minimum length of the lines in the image, given in
%       nanometers or as a string in the format 'xxx%' representing the 
%       percentage of the smallest image dimenstion. Optional, default =
%       '5%'.
%   line_max_length: maximum length of the lines in the image, given in
%       nanometers or as a string in the format 'xxx%' representing the 
%       percentage of the smallest image dimenstion. Optional, default =
%       '95%'.
% Outputs:
%   pdf_map: array of floating-point doubles, normalized sampling of the
%       analytical pdf.
%   line_start_coords/line_end_cords: coordinates of the start and end of 
%       the lines in the pdf map, given as an n by 2 matrix of (x, y)
%       values.

% Set defaults
if nargin < 5; line_min_length = '5%'; end;
if nargin < 6; line_max_length = '95%'; end;

% Rename parameter structure for convenience
params = parameter_struct;

% Determine coordinate bounds
min_x_bound = params.bounds(1);
min_y_bound = params.bounds(2);
max_x_bound = params.bounds(3);
max_y_bound = params.bounds(4);
x_length = max_x_bound - min_x_bound;
y_length = max_y_bound - min_y_bound;

% Determine line minimum length
if ~isnumeric(line_min_length)
    tokens = regexp(params.cell_radius, '(\d+)([%])', 'tokens'); % Look for a percentage value
    if ~isempty(tokens)
        base_length = min([x_length, y_length]);
        line_min_length = (str2num(tokens{1}{1})/100) * base_length;
    end
end

% Determine line max lengths
if ~isnumeric(line_max_length)
    tokens = regexp(params.cell_radius, '(\d+)([%])', 'tokens'); % Look for a percentage value
    if ~isempty(tokens)
        base_length = min([x_length, y_length]);
        line_max_length = (str2num(tokens{1}{1})/100) * base_length;
    end
end

% Double check the min and max are not mixed up
if line_min_length > line_max_length
    error('Minimum line length must be smaller or equal to the maximum length');
end

% Generate line coordinates, work in coordinates relative to lower left corner
% Initialize variables
line_start_coords = zeros(number_of_lines, 2);
line_end_coords = zeros(number_of_lines, 2);
line_counter = 0;

% Repeat loop until enough valid lines are generated
while line_counter < number_of_lines 
    line_start = rand(1,2) .* [x_length, y_length];
    line_end = rand(1,2) .* [x_length, y_length];
    line_length = sqrt(sum((line_end - line_start).^2);
    if line_length >= line_min_length && line_length <= line_max_length
       line_counter = line_counter + 1;
       line_start_coords(line_counter, :) = line_start;
       line_end_coords(line_counter, :) = line_end;
    end
end

%Calc needed number of pixels
map_resolution = params.ch1_distribution_params{2};
num_pixels_x = ceil(x_length/map_resolution);
num_pixels_y = ceil(y_length/map_resolution);

% Calc meshgrid
x_vec = [0.5: 1: num_pixels_x - 0.5] .* map_resolution;
y_vec = [num_pixels_y - 0.5: -1: 0.5] .* map_resolution;
[xmesh, ymesh] = meshgrid(x_vec, y_vec);

% Initalize pdf map image
pdf_map = zeros(num_pixels_y, num_pixels_x);

% Calculate each line seperately and add to total map
for line_index = 1:size(line_start_coords, 1)
    start_x = line_start_coords(line_index, 1);
    start_y = line_start_coords(line_index, 2);
    end_x = line_end_coords(line_index, 1);
    end_y = line_end_coords(line_index, 2);
    
    % Select a rectangular region of the meshgrid that must contain the entire line
    window_left_index = floor((min([start_x, end_x]) - line_width) / map_resolution) + 1;
    window_right_index = floor((max([start_x, end_x]) - line_width) / map_resolution);
    window_top_index = floor((min([start_y, end_y]) - line_width) / map_resolution) + 1;
    window_bottom_index = floor((max([start_y, end_y]) - line_width) / map_resolution);
    window_xmesh = xmesh(window_top_index:window_bottom_index, window_left_index:window_right_index);
    window_ymesh = ymesh(window_top_index:window_bottom_index, window_left_index:window_right_index);
    
    % Calc distances from line
    window_distance = abs((end_y - start_y) * window_xmesh - (end_x - start_x) * window_ymesh + end_x * start_y - end_y * start_x)...
        / sqrt((end_y - start_y).^2 + (end_x - start_x).^2);
    
    % Add 1 to the pdf map for the pixels that are within the specified width
    pdf_map(window_top_index:window_bottom_index, window_left_index:window_right_index) =...
        pdf_map((window_top_index:window_bottom_index, window_left_index:window_right_index) +...
        double(window_distance <= line_width);
end

% Normalize the pdf map
if ~isinf(line_to_background_ratio)
    pdf_map = pdf_map * line_to_background_ratio; % Multiply by ratio
    pdf_map = pdf_map + ones(size(pdf_map)); % Add one for background
end
map_sum = sum(pdf_map(:));
pdf_map = pdf_map / map_sum;

end


