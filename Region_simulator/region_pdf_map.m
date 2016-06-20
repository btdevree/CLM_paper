function [pdf_map, polygon_points] = region_pdf_map(parameter_struct, number_of_vertices, region_to_background_ratio)
%REGION_PDF_MAP Makes a pdf map of a randomized polygon in then center of
% the image.
%
% A polygon with randomized verterx positions. Assumes a Carteisain
% coordinate system. The origin is assumed to be in the lower left corner
% of the lower left pixel unless the bounds indicate otherwise in the 
% parameter structure.
%
% Inputs:
%   parameter_struct: parameter structure made with
%       test_movie_parameters_dv.
%   number_of_vertices: integer, number of vertices to create the polygon
%       with.
%   region_to_background_ratio: ratio of the density of events inside the
%       polygon region to the general background. Use Inf for a image with 
%       no background.
% Outputs:
%   pdf_map: array of floating-point doubles, normalized sampling of the
%       analytical pdf.
%   polygon_points: n by 2 matrix of (x, y) coordinates of the vertex 
%       points of the polygon region.

% Rename parameter structure for convenience
params = parameter_struct;

% Determine coordinate bounds
min_x_bound = params.bounds(1);
min_y_bound = params.bounds(2);
max_x_bound = params.bounds(3);
max_y_bound = params.bounds(4);
x_length = max_x_bound - min_x_bound;
y_length = max_y_bound - min_y_bound;

% Determine the angle for each vertex and size of polygon
base_angle = 2 * pi / number_of_vertices;
max_dist_from_center = min(0.9 * x_length, 0.9 * y_length);
max_base_radius = max_dist_from_center / 1.5; % 

% Generate line coordinates, work in coordinates relative to lower left corner
% Initialize variables
start_coords = zeros(number_of_lines, 2);
end_coords = zeros(number_of_lines, 2);
if strcmp(line_type, 'quadratic') || strcmp(line_type, 'cubic')
    cp1_coords = zeros(number_of_lines, 2);
end
if strcmp(line_type, 'cubic')
    cp2_coords = zeros(number_of_lines, 2);
end
line_counter = 0;

% Repeat loop until enough valid lines are generated
while line_counter < number_of_lines 
    line_start = rand(1,2) .* [x_length, y_length];
    line_end = rand(1,2) .* [x_length, y_length];
    chord_length = sqrt(sum((line_end - line_start).^2, 2));
    if chord_length >= chord_min_length && chord_length <= chord_max_length
        line_counter = line_counter + 1;
        start_coords(line_counter, :) = line_start;
        end_coords(line_counter, :) = line_end;
        if strcmp(line_type, 'quadratic') || strcmp(line_type, 'cubic')
            cp1_coords(line_counter, :) = rand(1,2) .* [x_length, y_length];
        end
        if strcmp(line_type, 'cubic')
            cp2_coords(line_counter, :) = rand(1,2) .* [x_length, y_length];
        end
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

% Initalize pdf_map and control point matricies
pdf_map = zeros(size(xmesh));
if strcmp(line_type, 'line_segment')
    control_points_x = zeros(number_of_lines, 2);
    control_points_y = zeros(number_of_lines, 2);
elseif strcmp(line_type, 'quadratic')
    control_points_x = zeros(number_of_lines, 3);
    control_points_y = zeros(number_of_lines, 3);
elseif strcmp(line_type, 'cubic')
    control_points_x = zeros(number_of_lines, 4);
    control_points_y = zeros(number_of_lines, 4);
end

% Calculate each line seperately and add to total map
for line_index = 1:size(start_coords, 1)
    % Report
    disp(['Working on line number ', num2str(line_index)]);
    
    % Prepare the control points for each type of line
    if strcmp(line_type, 'line_segment')
        cp = [start_coords(line_index, :); end_coords(line_index, :)];
    elseif strcmp(line_type, 'quadratic')
        cp = [start_coords(line_index, :); cp1_coords(line_index, :); end_coords(line_index, :)];
    elseif strcmp(line_type, 'cubic')
        cp = [start_coords(line_index, :); cp1_coords(line_index, :); cp2_coords(line_index, :); end_coords(line_index, :)];
    end
    
    % Record control points
    control_points_x(line_index, :) = cp(:, 1)';
    control_points_y(line_index, :) = cp(:, 2)';
    
    % Approximate the arc length
    testpoints = calc_bezier_line(cp, 25); % 25 points should most always us a length within 1% of the true value for cubic beziers
    approx_length = sum(sqrt(sum((testpoints(2:end, :) - testpoints(1:end-1, :)).^2 , 2)), 1);
    
    % Get approreate number of points along the bezier curve 
    number_curve_points = 5 * (approx_length / map_resolution); % about 1/10th pixel maximum errors on distance measurements
    curve_points = calc_bezier_line(cp, number_curve_points);
    
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


