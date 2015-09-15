function [distance_vector, mean_vector, stdev_vector, sem_vector] = radial_average_2D_correlation_binning(image, center_coords, max_distance)
%RADIAL_AVERAGE_2D_CORRELATION Calculates the radial average of a cross- or
%   auto-correlation image using a binning of pixel method
%
%   Calculates a radial average using assignment of pixels to the nearest 1
%   pixel wide ring. Uses a Cartesian coordinate system with the
%   origin at the lower left corner of the lower left pixel.
%   Inputs:
%   image: A correlation image given as a matrix of floating-point doubles.
%       Can also be a 3D array containing several correlation images to be 
%       averaged.
%   center_coords: Position of the center of the correlation image, given
%       as a column vector [x_coord; y_coord]. Optional, default = exact
%       middle of image.
%   max_distance: Maximum radial distance for the average calculation, 
%       given in pixels. Optional, default = radius of largest full circle 
%       that can be drawn.
%   Outputs:
%   distance_vector: Column vector of radial distances.
%   mean_vector: Column vector of the radial averages at each radius.
%   stdev_vector: Column vector of the standard deviation of each average 
%       value.

% Set defaults
if nargin < 2 
    x_coord = size(image, 2)/2;
    y_coord = size(image, 1)/2; 
    center_coords = [x_coord; y_coord];
end
if nargin < 3,
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

% Convert to polar coordiantes
centered_x = image_x - center_coords(1);
centered_y = image_y - center_coords(2);

% Generate list of radial distances
distance_vector = [0:max_distance].'; % Bin centers
bin_edges = [0, 0.5:1:max_distance + 0.5].'; % Bin edges 

% Initialize result vectors
mean_vector = zeros(size(distance_vector));
stdev_vector = zeros(size(distance_vector));
sem_vector = zeros(size(distance_vector));

% Loop through each image in the image stack
values = zeros(size(image, 1) * size(image, 2), size(image, 3));
for image_index = 1:size(image, 3)
    
    % Calculate the radius and get the value for each pixel on the image
    [~, radius, image_values] = cart2pol(centered_x, centered_y, image(:, :, image_index));
    values(:, image_index) = image_values(:); 
end
radius = radius(:);

% Loop through each bin
for bin_index = 1:length(distance_vector) 

    % Get binary indexes of pixels whose centers are inside the bin
    logical_index_vector = radius >= bin_edges(bin_index) & ...% inclusive of left bin edge
        radius < bin_edges(bin_index+1); % exclusive of right bin edge; mimics histc function
    
    % Get the subset of values that are going to be averaged
    value_subset = zeros(sum(logical_index_vector), size(image, 3));
    for image_index = 1:size(image, 3)
        value_subset(:, image_index) = values(logical_index_vector);
    end
    
    % Calculate the average and the standard deviation of the values
    mean_vector(bin_index) = mean(value_subset(:));
    stdev_vector(bin_index) = std(value_subset(:));
    cell_means = mean(value_subset, 2);
    sem_vector(bin_index) = std(cell_means)/sqrt(size(image, 3));
end
end

