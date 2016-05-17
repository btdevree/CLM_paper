function [ center, center_error, radius, radius_error ] = parametric_image_central_circle_finder(image)
%PARAMETRIC_IMAGE_CENTRAL_CIRCLE_FINDER Tries to find a centered circular region in the given image
%
% Radially averages the image intensity and finds the radius by searching 
% for the minimum of the first derivative. The value of the center is
% optimized by finding the minimum of the minimums, when the radial
% averages cross the border of steepest decent at all parts of the circle
% as the same radial value. Uses units of pixels in a Cartesian grid, with
% the orgin at the bottom left corner of the bottom left pixel.
%
% Input:
%   image - Image to look for the circular region in. May be sparse.
% Output:
%   center - optimized center value, given as a column matrix, in pixels.
%   center_error - optimized standard deviation of the center estimate, 
%       given as a column matrix, in pixels.
%   radius - optimized radius value, in pixels.
%   radius_error - optimized standard deviation of the radius estimate, in
%       pixels.

% Convert image to full, if needed
if issparse(image)
    image = full(image);
end

% Get initial estimates
initial_center = size(image)'./2;
initial_radius = find_radius(image, initial_center);

% Set up and run optimization
options = optimoptions('fminunc', 'Algorithm', 'quasi-newton', 'Display', 'iter-detailed', 'TolX', 1e-12);
initial_params = [initial_center(1, 1); initial_center(2, 1); initial_radius]; 
[optimized_params, ~, ~, ~, ~, params_hessian] = fminunc(@(params)find_nlog_nderivative(params, image), initial_params, options);

params_hessian

% Extract paramaters
center = optimized_params(1:2, :);
radius = optimized_params(3, :);

% Calculate error estimates from Hessian
errors = sqrt(diag(inv(params_hessian)));
center_error = errors(1:2, :);
radius_error = errors(3, :);
end

function [radius] = find_radius(image, center)
% Local function to get radius by finding the minimum of the first
% derivative of the radial average

% Take rough radial average with given center
[distance_vector, mean_vector] = radial_average_2D_correlation(image, 1, 1, [], center);

% Calculate derivative vector
d_dist = distance_vector(2:end) - distance_vector(1:end-1);
d_mean = mean_vector(2:end) - mean_vector(1:end-1);
derivative_vector = d_mean ./ d_dist;

% Minimum of derivative is at the radius value
[~, deriv_min_index] = min(derivative_vector);
radius = (distance_vector(deriv_min_index) + distance_vector(deriv_min_index + 1))/2;

% Take fine radial average with given center
[distance_vector, mean_vector] = radial_average_2D_correlation(image, .3, .3, [], center, [radius*.9:0.3:radius*1.1]);

% Calculate derivative vector
d_dist = distance_vector(2:end) - distance_vector(1:end-1);
d_mean = mean_vector(2:end) - mean_vector(1:end-1);
derivative_vector = d_mean ./ d_dist;

% Minimum of derivative is at the radius value
[~, deriv_min_index] = min(derivative_vector);
radius = (distance_vector(deriv_min_index) + distance_vector(deriv_min_index + 1))/2;
end

function [nlog_nderiv] = find_nlog_nderivative(params, image)
% Local function to get the negative log of the negative of the radial 
% average derivative value at the specified center and radial position.

% Extract parameters
center = params(1:2, :);
radius = params(3, :);

% Take radial average with given center and radial values
[distance_vector, mean_vector] = radial_average_2D_correlation(image, .3, .3, [], center, [radius - 0.15; radius + 0.15;]);

% Calculate derivative vector
d_dist = distance_vector(end) - distance_vector(1);
d_mean = mean_vector(2) - mean_vector(1);
derivative_value = d_mean ./ d_dist;
nlog_nderiv = -log(-derivative_value);
end