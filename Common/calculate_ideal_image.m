function [ideal_image, exact_image] = calculate_ideal_image(parameters_structure, doughnut_width)
%CALCULATE_IDEAL_IMAGE Calculate the ideal STORM image for the given test
%   movie parameters. Makes a circular region of constant density or a ring
%   of specified width
%
%   Note: This is not the best way to create ideal images; better to
%   supersample exact values, blur, and then take the appropreate pixels to
%   get the desired resolution. Works for the current applications, though.
% 
% Inputs:
%   parameter_structure: parameter structure generated with
%       test_movie_parameters_dv
%   doughnut_width: optional, the width of the ring in nanometers. The
%       middle of the ring is put at the cell_radius value.
%   
% Outputs:
%   ideal_image: the ideal, no noise image that would be made given the
%       options in the given parameter_structure

% Set defaults
if nargin < 2; doughnut_width = []; end;

% Rename parameters_structure for convenience
params = parameters_structure; 

% Determine coordinate bounds
min_x_bound = params.bounds(1);
min_y_bound = params.bounds(2);
max_x_bound = params.bounds(3);
max_y_bound = params.bounds(4);

% Determine cell center
if strcmp(params.cell_center, 'centered')
    cell_center_x = (max_x_bound - min_x_bound)/2;
    cell_center_y = (max_y_bound - min_y_bound)/2;
else
    cell_center_x = params.cell_center(1);
    cell_center_y = params.cell_center(2);
end

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

% Parameters needed for STORM image creating function
STORM_resolution = params.STORM_pixel_size; % data is given in nm, so resolution is just the nm per pixel
STORM_sigma = params.STORM_precision; % data is given in nm, so precision is in nm
STORM_dims = [max_y_bound - min_y_bound, max_x_bound - min_x_bound]; % data is given in nm, so bounds are in nm

% Calc image parameters
total_number_pixels_x = ceil(STORM_dims(2) / STORM_resolution);
total_number_pixels_y = ceil(STORM_dims(1) / STORM_resolution);

% Calc mesh x and y values in original pixel units
x_vector = STORM_resolution .* [0.5:1:total_number_pixels_x-0.5];
y_vector = STORM_resolution .* [total_number_pixels_y-0.5:-1:0.5];
[x_mesh, y_mesh] = meshgrid(x_vector, y_vector);
distance = sqrt((x_mesh - cell_center_x).^2 + (y_mesh - cell_center_y).^2);

% Initialize the image
ideal_image = zeros(total_number_pixels_y, total_number_pixels_x);

% Making a straighforward circle region
if isempty(doughnut_width)

    % Set to 1.0 if distance is within the circle
    ideal_image(distance <= cell_radius) = 1.0;

% Making a ring region
else
    
    % Set to a cosine intensity profile if distance is within the ring
    selection = (distance <= cell_radius + doughnut_width / 2) & (distance >= cell_radius - doughnut_width / 2);
    cosine_image = cos(((distance - cell_radius) / (doughnut_width / 2)) * (pi / 2));
    ideal_image(selection) = cosine_image(selection);
end

if nargout > 1
    exact_image = ideal_image;
end

% Gaussian filter with specified radius
STORM_sigma_in_pixels = STORM_sigma / STORM_resolution;
filter_kernel = fspecial('gaussian', ceil(7 * STORM_sigma_in_pixels), STORM_sigma_in_pixels);
ideal_image = imfilter(ideal_image, filter_kernel);

end

