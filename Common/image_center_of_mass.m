function [ CM_coords ] = image_center_of_mass(image, pixel_size, zero_baseline)
%IMAGE_CENTER_OF_MASS Finds center of mass of an image

% Assumes a Cartesian coordinate system, with the origin at the lower-left
%   corner of the lower-left pixel.
% Inputs:
%   image - matrix of image data, possibly sparse.
%   pixel_size - size of pixel, coordinates are returned in the same units.
%       Optional, default = 1;
%   zero_baseline - boolean flag, if true the minimum value in the image is
%       set to value zero. Optional, default = false.
% Output:
%   CM_coords = 2x1 matrix of floating-point doubles; [x_coord; y_coord].

% Set defaults
if nargin < 2; pixel_size = 1; end;
if nargin < 3; zero_baseline = false; end;

% Convert sparse image, if needed
if issparse(image)
    image = full(image);
end

% Zero baseline, if requested
if zero_baseline
    image = image - min(image(:));
end

% Get vectors of the intensity and the pixel indexes
intensity = image(:);
[row_index, column_index] = find(true(size(image)));

% Convert to cartesian coordinates
row_index = size(image, 1) - row_index + 0.5;
column_index = column_index - 0.5;

% Concatenate into a matrix 
index_coords = [column_index, row_index];

% Calculate first moment (mean)
intensity_sum = sum(intensity);
CM_coords = pixel_size * ((index_coords' * intensity)/intensity_sum);
end

