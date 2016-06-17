function [distances] = min_dist_to_curve(coords, curve)
%MIN_DIST_TO_CURVE Calculates the minimum distance between each of the x
% and y coordinates to any point on the sampled curves.
%
% Inputs:
%   coords: floating-point double n by 2 matrix of x and y coordinates to
%       measure the distances from.
%   curve: floating-point double m by 2 matrix of x and y coordinates that
%       define a curve to measure the coords against.
% Output:
%   distances: floating-point double n by 2 matrix of Eucledian distances 
%       from the points given in coords to any point contained in curve.\

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
    else
        distances = min(distances, new_distances);
    end
end    
end

