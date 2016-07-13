function [area] = triangle_area(point_1, point_2, point_3)
%TRIANGLE_AREA Area of triangle given its three vertices.

% Inputs: 
%   point_1/2/3: 1 x 2 floating point matrix with the (x, y) values of the
%       vertex coordinates
% Outputs:
%   area: area, given in the same units as the coordinates squared

% A = |(x1y2 + x2y3 + x3y1 – x1y3 – x2y1 – x3y2) / 2|
area = abs((point_1(1) * point_2(2) + point_2(1) * point_3(2) + point_3(1) * point_1(2) -...
        point_1(1) * point_3(2) - point_2(1) * point_1(2) - point_3(1) * point_2(2)) / 2);
end

