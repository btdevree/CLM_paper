function [ image ] = Gauss_STORM_image(xy_data, resolution, covar_inv, covar_det, calc_cutoff_pixels_x, calc_cutoff_pixels_y, x_vector, y_vector)
% Gauss_STORM_image Creates a STORM image from the list of (x,y) data; 
%   mirrors the C MEX function Gauss_STORM_image_MEX.cpp. 
% 
% Creates a STORM image as a matrix of floating-point doubles. Assumes a
%  Cartesian coordinate system, but the image does not need to contain the 
%  origin.
%
% Note: Type checking is done by the calling MATLAB function!!! The  
%  function can handle xy data points outside of the image boundries, but
%  most other errors are not checked. Don't give it bad values.
% 
% Inputs
%   xy_data: n by 2 array of floating-point doubles, center of each
%       gaussian pdf
%   resolution: floating-point double, number of nanometers per final image
%       pixel
%   covar_inv:  2 by 2 array of floating-point doubles, the inverse of the
%       2D gaussian covariance matrix
%   covar_inv: floating-point double, the determinant of the 2D gaussian
%       covariance matrix
%   calc_cutoff_pixels_x\y: 32 bit integer, extent in number of pixels on 
%       either side of the center pixel to calculate the pdf distribution
%       out to
%   total_number_pixels_x\y: 32 bit integer, number of pixels in the x and
%       y dimension of the final image.
%   x_vector\y_vector: Column vector arrays of floating-point doubles with 
%      the coordinate of the center of each pixel in the x or y direction.
%      Values are in increasing order (i.e., starting from the bottom-left 
%      corner of the Cartesian plane of the image). 
%   
% Outputs
%   image: array of floating-point doubles.

% Get the image dimensions
image_m = size(y_vector, 1);
image_n = size(x_vector, 1);

% Get coordinate values for the edges of the image
half_resolution = resolution / 2;
left_edge = min(x_vector) - half_resolution;
right_edge = max(x_vector) + half_resolution;
top_edge = max(y_vector) + half_resolution;
bottom_edge = min(y_vector) - half_resolution;

% Precalculate the 2D Gaussian scale factor
gauss_scale_factor = resolution^2 / (2 .* pi .* sqrt(covar_det));

% Initalize image array and loop through each point
image = zeros(image_m, image_n);
num_datapoints = size(xy_data, 1);
for data_ind = 1:num_datapoints
    
    % Get the x and y coordinate
    x = xy_data(data_ind, 1);
    y = xy_data(data_ind, 2);
    
    % Ignore points outside of the image boundries
    if x < left_edge || x > right_edge || y < bottom_edge || y > top_edge
        continue
    end
    
    % Calc index range for the point
    center_row = image_m - floor((y - bottom_edge)/resolution);
    center_column = ceil((x - left_edge)/resolution);
    min_row = center_row - calc_cutoff_pixels_y;
    max_row = center_row + calc_cutoff_pixels_y;
    min_column = center_column - calc_cutoff_pixels_x;
    max_column = center_column + calc_cutoff_pixels_x;

    % Make sure that the index range does not go off the edges of the image.
    if min_row < 1; min_row = 1; end
    if max_row > image_m; max_row = image_m; end
    if min_column < 1; min_column = 1; end
    if max_column > image_n; max_column = image_n; end

    % Calc x_mesh and y_mesh values in original pixel units
    [x_mesh, y_mesh] = meshgrid(x_vector(min_column:max_column, :), y_vector(min_row:max_row, :));
    X = [x_mesh(:).'; y_mesh(:).'];

    % Calc pdf
    X_shifted = X - repmat([x; y], 1, size(X, 2)); 
    pdf = gauss_scale_factor * exp(-0.5 .* sum(X_shifted.' * covar_inv .* X_shifted.', 2));

    % Add to total
    image(min_row:max_row, min_column:max_column) =...
        image(min_row:max_row, min_column:max_column) + reshape(pdf, size(x_mesh));
end
end

