function [pdf_map, control_points_x, control_points_y] = straight_lines_pdf_map(parameter_struct, number_of_lines, line_width, line_to_background_ratio, line_type, line_min_length, line_max_length)
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
%   line_type: string, type of lines to draw. Options are:
%       'line_segment' - straight line segments
%       'quadratic' - quadratic bezier (3 control points)
%       'cubic' - cubic bezier (4 control points)
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
    tokens = regexp(line_min_length, '(\d+)([%])', 'tokens'); % Look for a percentage value
    if ~isempty(tokens)
        base_length = min([x_length, y_length]);
        line_min_length = (str2num(tokens{1}{1})/100) * base_length;
    end
end

% Determine line max lengths
if ~isnumeric(line_max_length)
    tokens = regexp(line_max_length, '(\d+)([%])', 'tokens'); % Look for a percentage value
    if ~isempty(tokens)
        base_length = min([x_length, y_length]);
        line_max_length = (str2num(tokens{1}{1})/100) * base_length;
    end
end

% Double check that min and max are not mixed up
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
    line_length = sqrt(sum((line_end - line_start).^2, 2));
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
pdf_map = zeros(size(xmesh));

% Calculate each line seperately and add to total map
for line_index = 1:size(line_start_coords, 1)
    % Report
    disp(['Working on line number ', num2str(line_index)]);
    
    % Calc distances from line
    endpoints = [line_start_coords(line_index, :); line_end_coords(line_index, :)];
    line_length = sqrt(sum((endpoints(2, :) - endpoints(1, :)).^2 , 2));
    number_curve_points = 5 * (line_length / map_resolution); % 1/10th pixel maximum error
    curve_points = calc_bezier_line(endpoints, number_curve_points);
    
    % Select pixels next to the line, too inaccurate to use as distance measurements but avoids calculating many unused distances
    [coord_linear_indices, x_coords, y_coords] = select_bezier_region(curve_points, xmesh, ymesh, line_width/2);
    
    % Find the minimum distance for each coordinate to any point on the curve
    distances = distance_to_bezier(curve_points, x_coords, y_coords, false);
    
    % Add 1 to the pdf map for the pixels that are within the specified width
    pdf_map(coord_linear_indices) = pdf_map(coord_linear_indices) + double(distances <= line_width/2);
end

% Normalize the pdf map
if ~isinf(line_to_background_ratio) % Any non-perfect image
    pdf_map = pdf_map * line_to_background_ratio; % Multiply by ratio
    pdf_map = pdf_map + ones(size(pdf_map)); % Add one for background
end
map_sum = sum(pdf_map(:));
pdf_map = pdf_map / map_sum;
end


