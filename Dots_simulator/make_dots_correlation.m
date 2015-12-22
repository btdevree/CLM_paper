function [correlation_stack] = dots_correlation_individual(image, dot_center_coords, dot_radius, correlation_max, image_resolution, image_mask, image_origin)
%DOTS_CORRELATION_INDIVIDUAL Compute a cross-correlation between the center
%   of each individual dot and the given image. 
%
%   Assumes the coordinates are given in the Cartesian coordinate system.
%
% Inputs:
%   image: floating-point double STORM image. May be a sparse array.
%       Assumed to be zero padded and have the mask applied. 
%   dot_center_coords: n-by-2 array of (x, y) floating point double
%       coordinates.
%   dot_radius: sigma value for creating a Gaussian peak at the center of 
%       the dot region.
%   correlation_max: maximum correlation distance to return, in nanometers.
%   image_resolution: resolution of the image, in nanometers (or whatever 
%       unit the coordinates are given in). Optional, default = 1.
%   image_origin: 1-by-2 floating-point double array. The position of the
%       coodinate system origin relative to the bottom, left corner of the
%       the image's bottom, left pixel. Optional, default = (0, 0).
%   image_mask: floating-point double matrix with intensity range from 1 to
%       0. Assumed to be zero padded. Optional, defalult = all ones.
%   
% Output:
%   correlation_stack: 3D array of the n correlations for each individual 
%       dot region stacked along the 3rd dimension. Floating-point doubles.

% Set defaults
if nargin < 5; image_resolution = 1; end;
if nargin < 6; image_origin = [0, 0]; end;
if nargin < 7; image_mask = ones(size(image)); end;

% Convert sparse image, if needed
if issparse(image);
    image = full(image);
end

% Calculate the radius of the correlation in pixels 
max_radius = ceil(correlation_max / image_resolution);

% Initialize the results matrix
correlation_stack = zeros(2 * max_radius + 1, 2 * max_radius + 1, size(dot_center_coords, 1));

% Loop through each dot
for dot_index = 1:size(dot_center_coords, 1)
    dot_center = dot_center_coords(dot_index);
    
    % Create dot image
    


end

