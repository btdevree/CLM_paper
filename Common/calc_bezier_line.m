function [line_points] = calc_bezier_line(bezier_control_points, number_points)
%CALC_BEZIER_LINE Calculates points along a bezier line of order = 3 or
%   less (4 or less total control points).

% Inputs:
%   bezier_control_points: n by 2 matrix of (x, y) control points, n must
%       be 4, 3, or 2. 
%   number_points: number of points to draw along the bezier line. Minimum
%       = 2 points.
% Output:
%   number_points by 2 matrix of (x, y) control points.

% Check that the line is at least two points long
if number_points < 2; number_points = 2; end;

% Generate parameter vector
t = linspace(0, 1, number_points);

% Evalulate cubic bezier
if size(bezier_control_points, 1) == 4
    
    % Get points as column vectors
    p1 = bezier_control_points(1, :)';
    p2 = bezier_control_points(2, :)';
    p3 = bezier_control_points(3, :)';
    p4 = bezier_control_points(4, :)';
    
    % Compute as kronecker product
    line_points = (kron((1-t).^3, p1) + kron(3*(1-t).^2.*t, p2) + kron(3*(1-t).*t.^2, p3) + kron(t.^3, p4))';

% Evalulate quadratic bezier
elseif size(bezier_control_points, 1) == 3
    
    % Get points as column vectors
    p1 = bezier_control_points(1, :)';
    p2 = bezier_control_points(2, :)';
    p3 = bezier_control_points(3, :)';
    
    % Compute as kronecker product
    line_points = (kron((1-t).^2, p1) + kron(2*(1-t).*t, p2) + kron(t.^2, p3))';

% Evalulate line segment
elseif size(bezier_control_points, 1) == 2
    
    % Get points as column vectors
    p1 = bezier_control_points(1, :)';
    p2 = bezier_control_points(2, :)';
    
    % Compute as kronecker product
    line_points = (kron((1-t), p1) + kron(t, p2))';
end
end

