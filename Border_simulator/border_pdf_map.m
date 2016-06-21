function [pdf_map, border_coords] = border_pdf_map(parameter_struct, number_of_vertices, region_to_background_ratio)
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
%   polygon_coords: n by 2 matrix of (x, y) coordinates of the vertex 
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

% Determine the angle for each vertex and size of polygon and randomization radius
base_angle = 2 * pi / number_of_vertices;
max_dist_from_center = min(0.9 * x_length, 0.9 * y_length) / 2;
polygon_radius = max_dist_from_center / (1 + (sqrt(2 * (1 - cos(base_angle))) / 2)); % Max radius = polygon radius + 1/2 side length
randomization_radius = polygon_radius * (sqrt(2 * (1 - cos(base_angle))) / 2); % Randomization radius = 1/2 side length

% Calculate basic coordinates
center_coords = [x_length / 2, y_length / 2];
angles = base_angle * [0:number_of_vertices-1];
relative_x_coords = cos(angles)' * polygon_radius;
relative_y_coords = sin(angles)' * polygon_radius;
base_polygon_coords = [relative_x_coords, relative_y_coords] + repmat(center_coords, number_of_vertices, 1);

% Randomize polygon - random polar coords, so small displacements are preferred
rand_angles = 2 * pi * rand(number_of_vertices, 1);
rand_radii = randomization_radius * rand(number_of_vertices, 1);
rand_offsets = [cos(rand_angles), sin(rand_angles)] .* repmat(rand_radii, 1, 2);
polygon_coords = base_polygon_coords + rand_offsets;

%Calc needed number of pixels
map_resolution = params.ch1_distribution_params{2};
num_pixels_x = ceil(x_length/map_resolution);
num_pixels_y = ceil(y_length/map_resolution);

% Calc meshgrid
x_vec = [0.5: 1: num_pixels_x - 0.5] .* map_resolution;
y_vec = [num_pixels_y - 0.5: -1: 0.5] .* map_resolution;
[xmesh, ymesh] = meshgrid(x_vec, y_vec);

% Create pdf_map
pdf_map = zeros(size(xmesh));
pdf_map = pdf_map + double(inpolygon(xmesh, ymesh, polygon_coords(:, 1), polygon_coords(:, 2)));

% Normalize the pdf map
if ~isinf(region_to_background_ratio) % Any non-perfect image
    pdf_map = pdf_map * region_to_background_ratio; % Multiply by ratio
    pdf_map = pdf_map + ones(size(pdf_map)); % Add one for background
end
map_sum = sum(pdf_map(:));
pdf_map = pdf_map / map_sum;
end


