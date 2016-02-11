function [ mask ] = create_test_cell_STORM_mask(params, passed_variables)
%MAKE_TEST_CELL_STORM_MASK Creates a non-smoothend circular window for 
%   masking a test cell. 
%
%   Function that returns a double matrix with a circular region drawn in 
%       the cell region.

% Unpack required variables
x_length = passed_variables.max_x_bound - passed_variables.min_x_bound;
y_length = passed_variables.max_y_bound - passed_variables.min_y_bound;
cell_center_x = passed_variables.cell_center_x;
cell_center_y = passed_variables.cell_center_y;
cell_radius = passed_variables.cell_radius;
map_resolution = params.STORM_pixel_size;

%Calc needed number of pixels
num_pixels_x = ceil(x_length/map_resolution);
num_pixels_y = ceil(y_length/map_resolution);

% Calc meshgrid
x_vec = [0.5: 1: num_pixels_x - 0.5] .* map_resolution;
y_vec = [num_pixels_y - 0.5: -1: 0.5] .* map_resolution;
[xmesh, ymesh] = meshgrid(x_vec, y_vec);

% Calc delta distances from center
delta_x = xmesh - cell_center_x;
delta_y = ymesh - cell_center_y;
delta_dist = sqrt(delta_x.^2 + delta_y.^2);

% Returun image with 1 as for the pixels that are within the cell radius and 0 for those outside.
dist_bools = delta_dist <= cell_radius;
mask = double(dist_bools); 
end

