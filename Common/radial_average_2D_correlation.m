function [distance_vector, mean_vector, stdev_vector, sem_vector] = radial_average_2D_correlation(image, radial_sampling_distance, arc_sampling_distance, max_distance, center_coords, radial_values)
%RADIAL_AVERAGE_2D_CORRELATION Calculates the radial average of a cross- or
% auto-correlation image.
%
% Calculates a radial average using linear interpolation of a cross- or
% auto-correlation image. Uses a Cartesian coordinate system with the
% origin at the lower left corner of the lower left pixel. Units in
% pixels
%
% Inputs:
%   image: A correlation image given as a matrix of floating-point doubles.
%       Can also be a 3D array containing several correlation images to be 
%       averaged.
%   radial_sampling_distance: Distance along the radius inbetween sampling 
%       points. Optional, default = 1 pixels. 
%   arc_sampling_distance: Distance along the arc inbetween sampling points. 
%       Optional, default = 0.3 pixels. 
%   max_distance: Maximum radial distance for the average calculation, 
%       given in pixels. Optional, default = [], computes radius of largest 
%       full circle that can be drawn.
%   center_coords: Position of the center of the correlation image, given
%       as a column vector [x_coord; y_coord]. Optional, default = exact
%       middle of image.
%   radial_values: Calculate the radius at only the specified radial
%       values, given as a column vector. Optional, default = [], computes 
%       the entire vector from 0 to max radius.
% Outputs:
%   distance_vector: Column vector of radial distances, in pixels.
%   mean_vector: Column vector of the radial averages at each radius.

% Set defaults
if nargin < 2; radial_sampling_distance = 1; end;
if nargin < 3; arc_sampling_distance = 0.3; end;
if nargin < 4; max_distance = []; end;
if nargin < 5 
    x_coord = size(image, 2)/2;
    y_coord = size(image, 1)/2; 
    center_coords = [x_coord; y_coord];
end
if nargin < 6; radial_values = []; end; 
if isempty(max_distance)
    right_max = size(image, 2) - center_coords(1);
    left_max = center_coords(1);
    up_max = size(image, 1) - center_coords(2); 
    down_max = center_coords(2);
    max_distance = floor(min([right_max; left_max; up_max; down_max]));
end

% Generate meshgrid of image coordinates
x_vector = [0.5:1:size(image, 2)-0.5];
y_vector = [size(image, 1)-0.5:-1:0.5];
[image_x, image_y] = meshgrid(x_vector, y_vector);

% Generate list of radial distances
if isempty(radial_values);
    distance_vector = [0:radial_sampling_distance:max_distance].';
else
    distance_vector = radial_values;
end

% Initialize result vectors
mean_vector = zeros(size(distance_vector));

% Loop through each radius
for rad_ind = 1:length(distance_vector) 
    radius = distance_vector(rad_ind);
    
    % Special case for radius == 0
    if radius == 0
        x_circle = 0;
        y_circle = 0;
    else
    
        % Spread out sampling arc evenly along circumference from 0 to 2*pi.
        total_circ_length = 2 .* pi .* radius; % circumfrence in pixels
        num_points = ceil(total_circ_length ./ arc_sampling_distance);
        actual_radian_interval = 2 .* pi ./ num_points;
        radians = actual_radian_interval * [0:num_points-1].';
        [x_circle, y_circle] = pol2cart(radians, radius);
    end
    
    % Convert to x,y coordinates
    x_coords = x_circle + center_coords(1);
    y_coords = y_circle + center_coords(2);
    
    % Loop through each image in the image stack
    values = zeros(size(x_coords, 1), size(image, 3));
    for image_index = 1:size(image, 3)
        
        % Calculate the value at each coordinate with linear interpolation
        values(:, image_index) = interp2(image_x, image_y, image(:, :, image_index), x_coords, y_coords);
    end
    
    % Calculate the average of the values
    mean_vector(rad_ind) = mean(values(:));
end
end

