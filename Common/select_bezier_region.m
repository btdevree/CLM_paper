function [linear_indicies, x_coords, y_coords] = select_bezier_region(bezier_curve_points, xmesh, ymesh, target_distance)
%SELECT_BEZIER_REGION Select coordinates for pixels that lie next to a
% bezier curve.
%
% Returns indicies and coordinate values for the points that lie next to
% the curve defined by the given set of points belonging to the provided
% bezier curve. Uses morphology operations, so not very exact. Returns a 
% region that contains all distances up to the target distance. Assumes
% square pixels.
%
% Inputs:
%   bezier_curve_points: floating-point double n by 2 matrix of x and y 
%       coordinates sampled along bezier curve.
%   xmesh/ymesh: coordinates of the image pixels, generated with meshgrid.
%   target_distance: selected region will include all pixel coordinates at
%       or less than this distance to the curve. Given in the same units as
%       the xmesh/ymesh coordinates.
% Output:
%   row_indices: column vector with all the row index values (1st 
%       dimenstion) for the returned coordinates.
%   column_indices: column vector with all the column index values (2nd 
%       dimenstion) for the returned coordinates.
%   x_coords: column vector with all the x coordinate values.
%   y_coords: column vector with all the y coordinate values.

% Calculate the pixel size
pixel_size = xmesh(1, 2) - xmesh(1, 1);

% Convert x and y coordinates into index values
image_height = size(ymesh, 1);
curve_row_indices = ceil(image_height - (bezier_curve_points(:, 2) / pixel_size));
curve_column_indices = ceil(bezier_curve_points(:, 1) / pixel_size);
curve_linear_indices = curve_row_indices + image_height * (curve_column_indices - 1);

% Initalize an empty binary matrix and mark pixels that contain a curve point
map = false(size(xmesh));
map(curve_linear_indices) = true;

% Dialate the marked pixels
number_cycles = ceil(1.08 * (target_distance / pixel_size)) + 1; % Standard sequence dialates 2.8 pixels in the diagnal directions for every 3 cycles. Hit target, round up, and add one to gaurentee hitting target distance. 
map = standard_sequence_dilate(map, number_cycles);

% Get the index values and coordinate values of each marked pixel
linear_indicies = find(map);
x_coords = xmesh(map);
y_coords = ymesh(map);
end

