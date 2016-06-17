function [distances] = distance_to_bezier(bezier_curve_points, x_coords, y_coords, use_MEX_flag)
%DISTANCE_TO_BEZIER Approximates the distance between a sampled line segment
% and the given (x,y) coordinates.

% Calculate the minimium distance between a set of x and y coordinates and
% any of the supplied bezier curve points. Function is a wrapper that calls 
% min_dist_to_curve or min_dist_to_curve_MEX, depending on the 
% use_MEX_flag.

% Inputs:
%   bezier_curve_points: floating-point double n by 2 matrix of x and y 
%       coordinates sampled along bezier curve.
%   x_coords/y_coords: matrix of x and y coordinate values, can be any 
%       shape (i.e., meshgrid format or just vectors), but x and y need to 
%       be the same length and shape.
%   use_MEX_flag: Binary flag to perform the distance calculations with a 
%       MEX routine. Optional, default = true. 
% Outputs:
%   distances: Minimum distance between the coordinates and any testpoint 
%       on the length of the curve. Returned as the same shape as the 
%       x_coords argument.

% Set defaults
if nargin < 4; use_MEX_flag = true; end;

% Give list of points to the core distance measurement function
if ~use_MEX_flag % Don't use C++ MEX routine
    
    % Run measurement function
    distances = reshape(min_dist_to_curve([x_coords(:), y_coords(:)], bezier_curve_points), size(x_coords));
    
else % Use the C++ MEX routine
    
    % Check inputs to avoid segfaults
    
    % Run measurement MEX function
    distances = reshape(min_dist_to_curve_MEX([x_coords(:), y_coords(:)], bezier_curve_points), size(x_coords));
    
end

