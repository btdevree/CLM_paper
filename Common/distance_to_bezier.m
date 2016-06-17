function [distances] = distance_to_bezier(bezier_control_points, x_coords, y_coords, number_line_testpoints, use_MEX_flag)
%DISTANCE_TO_BEZIER Approximates the distance between a bezier line segment
% and the given (x,y) coordinates.

% Calculate the minimium distance between a set of x and y coordinates and
% a bezier curve with the specified control points. Use two control points 
% for a simple line segment, can hande up to 4 control points. Approximates 
% the measure by testing points all along the line with the specified 
% coordinates, so the true minimum is slightly lower in value.

% Inputs:
%   bezier_control_points: n by 2 matrix of (x, y) control points, n must
%       be 4, 3, or 2. 
%   x_coords/y_coords: matrix of x and y coordinate values, can be any 
%       shape (i.e., meshgrid format or just vectors), but x and y need to 
%       be the same length and shape.
%   number_line_testpoints: the number of points to test along the length 
%       of the curve.
%   use_MEX_flag: Binary flag to perform the distance calculations with a 
%       MEX routine. Optional, default = true. 
% Outputs:
%   distances: Minimum distance between the coordinates and any testpoint 
%       on the length of the curve. Returned as the same shape as the 
%       x_ccords argument.

% Set defaults
if nargin < 5; use_MEX_flag = true; end;

% Get lists of testpoints
testpoints = calc_bezier_line(bezier_control_points, number_line_testpoints);

% Give list of points to the core distance measurement function
if ~use_MEX_flag % Don't use C++ MEX routine
    
    % Run measurement function
    distances = reshape(min_dist_to_curve([x_coords(:), y_coords(:)], testpoints), size(x_coords));
    
else % Use the C++ MEX routine
    
    % Check inputs to avoid segfaults
    
    % Run measurement MEX function
    distances = reshape(min_dist_to_curve_MEX([x_coords(:), y_coords(:)], testpoints), size(x_coords));
    
end

