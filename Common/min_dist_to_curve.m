function [distances, curve_indices] = min_dist_to_curve(coords, curve, keep_index_flag)
%MIN_DIST_TO_CURVE Calculates the minimum distance between each of the x
% and y coordinates to any point on the sampled curves.
%
% Inputs:
%   coords: floating-point double n by 2 matrix of x and y coordinates to
%       measure the distances from.
%   curve: floating-point double m by 2 matrix of x and y coordinates that
%       define a curve to measure the coords against.
%   keep_index_flag: if set to true, the index of the curve point that is
%       closest to each coord value is returned as well. Optional, default
%       = false;
% Output:
%   distances: floating-point double n by 2 matrix of Eucledian distances 
%       from the points given in coords to any point contained in curve.
%   curve_indices; n by 2 matrix of the row index value for the curve point
%       that is closest to each coord value.

% Set default
if nargin < 3; keep_index_flag = false; end;

% Get lengths of matrices
coord_length = size(coords, 1);
curve_length = size(curve, 1);

% Calculate the distances from all coord points to the curve points
for curve_index = 1:curve_length
    
    % Find distances to the current curve point
    curve_point = curve(curve_index, :);
    new_distances = sqrt(sum((repmat(curve_point, coord_length, 1) - coords).^2, 2));
    
    % Keep the distance measurement if it's smaller than the last
    if curve_index == 1
        distances = new_distances; % Copy all on the first iteration
        if keep_index_flag
            curve_indices = repmat(curve_index, size(new_distances));
        end
    else
        if keep_index_flag
            [distances, min_indices] = min([distances, new_distances], [], 2);
            curve_indices(min_indices == 2) = curve_index;
        else
            distances = min(distances, new_distances);
        end
    end
end    
end

