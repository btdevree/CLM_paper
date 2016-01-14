function [ pdf_map, center_coords, STORM_mask] = dots_pdf_map(number_of_dots, dot_radius, dot_to_background_ratio, parameter_struct)
%DOTS_PDF_MAP Makes a pdf map of circular regions with uniform density 
%   inside the main cell density of a simulated movie dataset.
%
%   Each dot is randomly placed within the boundries of the simulated cell. 
%
% Inputs:
%   number_of_dots: number of regions to put in the pdf map.
%   dot_radius: radius of each circular region.
%   dot_to_background_ratio: ratio of the density of events inside the dot
%       region to the general cellular background
%   parameter_struct: parameter structure made with
%       test_movie_parameters_dv.
% Outputs:
%   pdf_map: array of floating-point doubles, normalized sampling of the
%       analytical pdf.
%   center_coords: coordinates of the centers of the circular regions
%   STORM_mask: If requested, output a mask image that covers the simulated 
%       cell region in the STORM image. Not necessarily the same resolution 
%       as the mask. Default = false;

% Find the cell center, radius, points in the circle - copied from create_test_data_dv
%   Not ideal, but not sure how I want to functionalize/abstract this part yet. 

% Rename parameter structure for convenience
params = parameter_struct;

% Determine coordinate bounds
min_x_bound = params.bounds(1);
min_y_bound = params.bounds(2);
max_x_bound = params.bounds(3);
max_y_bound = params.bounds(4);
x_length = max_x_bound - min_x_bound;
y_length = max_y_bound - min_y_bound;

% Determine cell center
if strcmp(params.cell_center, 'centered')
    cell_center_x = (x_length)/2;
    cell_center_y = (y_length)/2;
else
    cell_center_x = params.cell_center(1);
    cell_center_y = params.cell_center(2);
end
cell_center = [cell_center_x, cell_center_y];

% Determine cell radius
if isnumeric(params.cell_radius)
    cell_radius = params.cell_radius;
else
    tokens = regexp(params.cell_radius, '(\d+)([%])', 'tokens'); % Look for a percentage value
    if ~isempty(tokens)
        max_radius = min([max_x_bound - cell_center_x, cell_center_x - min_x_bound,...
            max_y_bound - cell_center_y, cell_center_y - min_y_bound]);
        cell_radius = (str2num(tokens{1}{1})/100) * max_radius;
    end
end

% Keep all the pdf density inside the cell
cell_radius_for_dots = cell_radius - dot_radius;
if cell_radius_for_dots <= 0
    error('Cell is not large enough to hold circular regions of requested size');
end

% Initialize variables
center_coords = zeros(number_of_dots, 2);
event_counter = 0;

% Repeat loop until enough random events inside the circle are generated
while event_counter < number_of_dots 
    event_coord = (rand(1,2) * cell_radius_for_dots * 2) - cell_radius_for_dots;
    if sqrt(sum(event_coord.^2)) <= cell_radius_for_dots
       event_counter = event_counter + 1;
       center_coords(event_counter, :) = event_coord + cell_center;
    end
end

% Get the cell map
[pdf_map, xmesh, ymesh] = create_cell_map(params.ch1_distribution_params{2}, x_length, y_length, cell_center_x, cell_center_y, cell_radius);

% Multiply all pixels inside spot radius by deisred dot to background ratio
for dot_index = 1:size(center_coords, 1)
    center_x = center_coords(dot_index, 1);
    center_y = center_coords(dot_index, 2);
    
    % Calc delta distances from center
    delta_x = xmesh - center_x;
    delta_y = ymesh - center_y;
    delta_dist = sqrt(delta_x.^2 + delta_y.^2);
    
    % Add 1 to the pdf map for the pixels that are within the radius
    factors = ones(size(xmesh));
    factors(delta_dist <= dot_radius) = dot_to_background_ratio;
    pdf_map = pdf_map .* factors;
end

% Normalize the pdf map
map_sum = sum(pdf_map(:));
pdf_map = pdf_map / map_sum;

% Copy to the image mask if requested - Move this to another funciton in
% the next refactor, it doesn't really fit that well here.
if nargout == 3
    STORM_mask = create_cell_map(params.STORM_pixel_size, x_length, y_length, cell_center_x, cell_center_y, cell_radius);
end

end

function [cell_map, xmesh, ymesh] = create_cell_map(map_resolution, x_length, y_length, cell_center_x, cell_center_y, cell_radius)
% Local function that returns a double matrix with a circular region drawn 
%   in the cell region and the pixel coordinate meshgrid matrices.

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
cell_map = double(dist_bools); 
end
