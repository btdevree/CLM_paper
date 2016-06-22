function [pdf_map, border_coords] = border_pdf_map(parameter_struct, displacement_factor, roughness_parameter, number_of_cycles, light_to_dark_ratio)
%BORDER_PDF_MAP Makes a pdf map of a complex light/dark border down the
% center of the image.
%
% Draws a light/dark border with a midpoint displacement algorithm using a
% Gaussian distribution of displacement lengths.
%
% Inputs:
%   parameter_struct: parameter structure made with
%       test_movie_parameters_dv.
%   displacement_factor: standard deviation of the midpoint displacement 
%       randomization values, given as a fraction of the full border
%       length.
%   roughness_parameter: exponent for the scaling factor (1/2)^x
%   number_of_cycles: number of midpoint division cycles to run
%   light_to_dark_ratio: ratio of the density of events on the light side 
%       of the border to the dark side. Use Inf for a image with no dark
%       background.
% Outputs:
%   pdf_map: array of floating-point doubles, normalized sampling of the
%       analytical pdf.
%   border_coords: n by 2 matrix of (x, y) coordinates of the points that 
%       define the border.

% Rename parameter structure for convenience
params = parameter_struct;

% Determine coordinate bounds
min_x_bound = params.bounds(1);
min_y_bound = params.bounds(2);
max_x_bound = params.bounds(3);
max_y_bound = params.bounds(4);
x_length = max_x_bound - min_x_bound;
y_length = max_y_bound - min_y_bound;

% Get a value for the standard deviation of the border displacements
displacement_stdev = y_length * displacement_factor;

% Initalize list of border coordinates
border_coords = [x_length / 2, y_length; x_length / 2, 0];

% Split and randomize border for specified number of cycles
for cycle_index = 1:number_of_cycles
    
    % Get midpoints of the curve
    midpoint_coords = (border_coords(1:end - 1, :) + border_coords(2:end, :)) / 2;
    
    % Displace midpoints
    displacements = 2.^(-roughness_paramater * cycle_index - 1) * normrnd(0, displacement_stdev);
    midpoint_coords(:, 2) = midpoint_coords(:, 2) + displacements;
    
    % Interweave original and midpoint values
    new_coords = zeros(4, size(border_coords, 1));
    new_coords(1:2, :) = border_coords.';
    new_coords(3:4, 1:end-1) = midpoint_coords.';
    border_coords = reshape(new_coords(1:end-2), size(border_coords, 1) + size(midpoint_coords, 1), 2);
end

% Calc needed number of pixels
map_resolution = params.ch1_distribution_params{2};
num_pixels_x = ceil(x_length/map_resolution);
num_pixels_y = ceil(y_length/map_resolution);

% Calc grid vectors and meshgrids
x_vec = [0.5: 1: num_pixels_x - 0.5] .* map_resolution;
y_vec = [num_pixels_y - 0.5: -1: 0.5] .* map_resolution;
[xmesh, ~] = meshgrid(x_vec, y_vec);

% Interpolate border onto the y_vec values
border_crossing_x = interp1(border_coords(:, 2), border_coords(:, 1), y_vec);

% Create pdf_map
pdf_map = zeros(length(y_vec), length(x_vec));
pdf_map = pdf_map + double(xmesh - repmat(border_crossing_x', 1, size(xmesh, 2)) >= 0);

% Normalize the pdf map
if ~isinf(region_to_background_ratio) % Any non-perfect image
    pdf_map = pdf_map * region_to_background_ratio; % Multiply by ratio
    pdf_map = pdf_map + ones(size(pdf_map)); % Add one for background
end
map_sum = sum(pdf_map(:));
pdf_map = pdf_map / map_sum;
end


