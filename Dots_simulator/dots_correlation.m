function [correlation] = dots_correlation(STORM_image, dot_center_coords, dot_center_sigma, correlation_max, image_resolution, image_mask, image_origin)
%DOTS_CORRELATION Compute a cross-correlation between the center
%   of each individual dot and the given image. 
%
%   Assumes the coordinates are given in the Cartesian coordinate system.
%
% Inputs:
%   image: floating-point double STORM image. May be a sparse array. Does
%       not have mask applied.
%   dot_center_coords: n-by-2 array of (x, y) floating point double
%       coordinates.
%   dot_center_sigma: sigma value for creating a Gaussian peak at the 
%       center of the dot region.
%   correlation_max: maximum correlation distance to return, in nanometers.
%   image_resolution: resolution of the image, in nanometers (or whatever 
%       unit the coordinates are given in). Optional, default = 1.
%   image_mask: floating-point double matrix with intensity range from 1 to
%       0. Optional, defalult = all ones.
%   image_origin: 1-by-2 floating-point double array. The position of the
%       coodinate system origin relative to the bottom, left corner of the
%       the image's bottom, left pixel. Optional, default = (0, 0).
%   
% Output:
%   correlation: 2D array of the correlation between the dot regions and 
%       the labels in the STORM image. Floating-point doubles.

% Set defaults
if nargin < 5; image_resolution = 1; end;
if nargin < 6; image_mask = ones(size(STORM_image)); end;
if nargin < 7; image_origin = [0, 0]; end;

% Convert sparse image, if needed
if issparse(STORM_image);
    STORM_image = full(STORM_image);
end

% Calculate needed constants
max_radius = ceil(correlation_max / image_resolution);
pixel_dims = [size(STORM_image, 2), size(STORM_image, 1)];
dims =  pixel_dims * image_resolution;

% Create dot image 
data = struct();
data.x = dot_center_coords(:, 1) + image_origin(1, 1);
data.y = dot_center_coords(:, 2) + image_origin(1, 2);
dot_image = create_STORM_image(data, image_resolution, dot_center_sigma, dims);

% Calculate the correlation
correlation= calc_crosscorrelation(STORM_image, dot_image, max_radius, image_mask);
end

