function [ pdf_map, center_coords] = dots_pdf_map(parameter_struct, number_of_dots, dot_radius, dot_to_background_ratio)
%DOTS_PDF_MAP Makes a pdf map of circular regions with uniform density 
%   randomly distributed in an image.
%
% Inputs:
%   number_of_dots: number of regions to put in the pdf map.
%   dot_radius: radius of each circular region.
%   dot_to_background_ratio: ratio of the density of events inside the
%       dot region to the general background. Use Inf for a image with no 
%       background.
%   parameter_struct: parameter structure made with
%       test_movie_parameters_dv.
% Outputs:
%   pdf_map: array of floating-point doubles, normalized sampling of the
%       analytical pdf.
%   center_coords: coordinates of the centers of the circular regions

% Rename parameter structure for convenience
params = parameter_struct;

% Determine coordinate bounds
min_x_bound = params.bounds(1);
min_y_bound = params.bounds(2);
max_x_bound = params.bounds(3);
max_y_bound = params.bounds(4);
x_length = max_x_bound - min_x_bound;
y_length = max_y_bound - min_y_bound;

%Calc needed number of pixels
STORM_resolution = params.STORM_pixel_size;
map_resolution = params.ch1_distribution_params{2};
num_STORM_pixels_x = ceil(x_length/STORM_resolution);
num_STORM_pixels_y = ceil(y_length/STORM_resolution);
num_pixels_x = ceil(num_STORM_pixels_x * STORM_resolution / map_resolution);
num_pixels_y = ceil(num_STORM_pixels_y * STORM_resolution / map_resolution);

% Calc meshgrid
x_vec = [0.5: 1: num_pixels_x - 0.5] .* map_resolution;
y_vec = [num_pixels_y - 0.5: -1: 0.5] .* map_resolution;
[xmesh, ymesh] = meshgrid(x_vec, y_vec);

% Initalize pdf_map and center point matricies
pdf_map = zeros(size(xmesh));
center_coords_x = rand(number_of_dots, 1) * (x_length - 2 * dot_radius) + dot_radius;
center_coords_y = rand(number_of_dots, 1) * (y_length - 2 * dot_radius) + dot_radius;
center_coords = [center_coords_x, center_coords_y];

% Multiply all pixels inside spot radius by deisred dot to background ratio
for dot_index = 1:number_of_dots
    center_x = center_coords(dot_index, 1);
    center_y = center_coords(dot_index, 2);
    
    % Calc delta distances from center
    delta_x = xmesh - center_x;
    delta_y = ymesh - center_y;
    delta_dist = sqrt(delta_x.^2 + delta_y.^2);
    
    % Add 1 to the pdf map for the pixels that are within the radius
    pdf_map = pdf_map + double(delta_dist <= dot_radius); 
end

% Normalize the pdf map
if ~isinf(dot_to_background_ratio) % Any non-perfect image
    pdf_map = pdf_map * dot_to_background_ratio; % Multiply by ratio
    pdf_map = pdf_map + ones(size(pdf_map)); % Add one for background
end
map_sum = sum(pdf_map(:));
pdf_map = pdf_map / map_sum;
end