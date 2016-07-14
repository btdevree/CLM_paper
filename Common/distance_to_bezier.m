function [distances, closest_bezier_indices] = distance_to_bezier(bezier_curve_points, x_coords, y_coords, use_MEX_flag)
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
%   closest_bezier_indices: the row index value of the bezier point that is
%       closest to the specified x and y coordinates. Optional.

% Set defaults
if nargin < 4; use_MEX_flag = true; end;

% Give list of points to the core distance measurement function
if ~use_MEX_flag % Don't use C++ MEX routine
    
    % Run measurement function
    if nargout > 1
        [distances_linear, closest_bezier_indices_linear] = min_dist_to_curve([x_coords(:), y_coords(:)], bezier_curve_points, true); % Set flag to return the indices as well
        distances = reshape(distances_linear, size(x_coords));
        closest_bezier_indices = reshape(closest_bezier_indices_linear, size(x_coords));
    else
        distances = reshape(min_dist_to_curve([x_coords(:), y_coords(:)], bezier_curve_points), size(x_coords));
    end
    
else % Use the C++ MEX routine
    
    % MEX requires doubles
    if ~isa(bezier_curve_points, 'double')
        bezier_curve_points = double(bezier_curve_points);
    end
    if ~isa(x_coords, 'double')
        x_coords = double(x_coords);
    end
    if ~isa(y_coords, 'double')
        y_coords = double(y_coords);
    end
    
    % Prep coords for measurement function
    coords_matrix = [x_coords(:), y_coords(:)];
    
    % Type checking to make sure we don't end up giving the MEX file bad inputs and cause a seg-fault
    assert(ismatrix(coords_matrix) && size(coords_matrix, 2) == 2 && (isa(coords_matrix, 'double')));
    assert(ismatrix(bezier_curve_points) && size(bezier_curve_points, 2) == 2 && (isa(bezier_curve_points, 'double')));
    
    if nargout > 1
         % Run measurement MEX function with index tracking
        [distances_linear, closest_bezier_indices_linear] = min_dist_to_curve_with_indices_MEX(coords_matrix, bezier_curve_points);
        distances = reshape(distances_linear, size(x_coords));
        closest_bezier_indices = reshape(closest_bezier_indices_linear, size(x_coords));
    else
        % Run measurement MEX function
        distances = reshape(min_dist_to_curve_MEX(coords_matrix, bezier_curve_points), size(x_coords));
    end
end

