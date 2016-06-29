function [pdf_map, border_coords] = border_pdf_map(parameter_struct, displacement_factor, roughness_parameter,...
    light_to_background_ratio, number_of_cycles)
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
%   roughness_parameter: exponent for the scaling factor H^n for n cycles, 
%       typical range from .3 to .7. Choose higher values for a rougher 
%       curve and lower values for a smoother curve. 
%   light_to_background_ratio: ratio of the density of events on the light 
%       side of the border to the general background density. Use Inf for 
%       an image with no dark background.
%   number_of_cycles: number of midpoint division cycles to run. Optional,
%       default = [], which is enough cycles to randomize down to a 
%       single-pixel level.
% Outputs:
%   pdf_map: array of floating-point doubles, normalized sampling of the
%       analytical pdf.
%   border_coords: n by 2 matrix of (x, y) coordinates of the points that 
%       define the border.

% Set defaults
if nargin < 5; number_of_cycles = []; end;

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

% Calc number of cycles to run, if needed
if isempty(number_of_cycles)
    number_of_cycles = ceil(log(num_pixels_y - 1) / log(2)) + 2; % # points = 2^n + 1 for n cycles 
end
    
% Get a value for the standard deviation of the border displacements
displacement_stdev = y_length * displacement_factor;

% Initalize list of border coordinates
border_coords = [x_length / 2, y_length; x_length / 2, 0];

% Split and randomize border for specified number of cycles
for cycle_index = 1:number_of_cycles
    
    % Get midpoints of the curve
    midpoint_coords = (border_coords(1:end - 1, :) + border_coords(2:end, :)) / 2;
    
    % Displace midpoints
    displacements = roughness_parameter.^(cycle_index - 1) * normrnd(0, displacement_stdev, size(midpoint_coords, 1), 1);
    midpoint_coords(:, 1) = midpoint_coords(:, 1) + displacements;
    
    % Interweave original and midpoint values
    new_coords = zeros(4, size(border_coords, 1));
    new_coords(1:2, :) = border_coords.';
    new_coords(3:4, 1:end-1) = midpoint_coords.';
    new_coords = reshape(new_coords(:), 2, 2 * size(border_coords, 1));
    border_coords = new_coords(:, 1:end-1)';
end

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
if ~isinf(light_to_background_ratio) % Any non-perfect image
    pdf_map = pdf_map * light_to_background_ratio; % Multiply by ratio
    pdf_map = pdf_map + ones(size(pdf_map)); % Add one for background
end
map_sum = sum(pdf_map(:));
pdf_map = pdf_map / map_sum;
end


