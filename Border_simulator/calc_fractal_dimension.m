function [fractal_dim, fractal_dim_stdev] = calc_fractal_dimension(line_points, end_behavior, number_repeats)
%CALC_FRACTAL_DIMENSION Estimates the fractal dimension of a linear series
% of (x, y) points.
%
% Estimates the fractal dimension with a divider type algorithm.
%
% Input:
%   line_points: n by 2 matrix of (x, y) points.
%   end_behavior: string representing the method used to deal with the
%       incomplete measurement segments at the ends of the line. Choices
%       are: 
%       'finite': Assume the sequence simply ends and we draw segments 
%           directly to the first and last points. 
%       'circular': Sequence is circular and can be repeated to get 
%           more points at the beginning and end as necessary. Assumes the
%           first and last points overlap as the same point, and returned
%           measurements are over one period of points only.
%       Optional, default = 'finite'
%   number_repeats: number of times to repeat the calculation on different
%       starting points. Optional, default = 30.
% Output:
%   fractal_dim: estimate of the fractal dimension.
%   fractal_dim_stdev: estimate of the standard deviation of the estimated 
%       fractal dimension, uses Fisher information from fitting
%       Richardson's plot.

% Set defaults
if nargin < 2; end_behavior = 'finite'; end;
if nargin < 3; number_repeats = 30; end;

% Set parameters
number_divider_lengths = 10;

% Determine divider sizes
total_number_points = size(line_points, 1);
maximum_length = sum(sqrt(sum((line_points(2:end, :) - line_points(1:end-1, :)).^2, 2)), 1);
min_divider = 0.5 * maximum_length / total_number_points; % Half the mean step distance
max_divider = 0.1 * maximum_length; % 1/10th of max length
divider_lengths = logspace(log10(min_divider), log10(max_divider), number_divider_lengths);

% Initalize distance results matrix
curve_distance = zeros(number_divider_lengthes, number_repeats);

% Repeat the measurement process from multiple starting points
for repeat_index = 1:number_repeats
    
    % Choose a point at random to start at
    starting_point_index = randi(total_number_points);
        
    % Repeat measurement for all divider lengths
    for divider_index = 1:number_divider_lengths
       curve_distance(divider_index, repeat_index) = calc_distance(line_points, divider_length, starting_point_index, end_behavior);
    end
end

end

function [distance, divider_points] = calc_distance(line_points, divider_length, start_index, end_behavior)
% Local function to calculate the distance of a curve for a particular
% divider length. Optional output gives all the divider points.

% See if we'll need to log points
if nargout > 1
    log_points_flag = true;
else
    log_points_flag = false;
end

% Get starting coordinates
start_coords = line_points(start_index, :);
if log_points_flag
    divider_points = start_coords;
end

% ---- Count number of divider lengths backwards ----

% Initalize flags, counters, and coordinates
back_count = 0;
reached_beginning_flag = false;
current_index = start_index;
current_divider_index = start_index;
current_coords = start_coords;
current_divider_coords = start_coords;
first_divider_coords = [];

% Loop through divider lengths until the line points are used up
while ~reached_beginning_flag
    
    % Reset current step distance
    current_distance = 0; 
    
    % Go through line points until the step length is larger than the divider
    while current_distance < divider_length
        
        % Try to get values at the previous index
        try % Usually, get the coordinates of the previous point on the line
            current_coords = line_points(current_index - 1, :);
            current_index = current_index - 1; % Will not evaluate if error is thrown
        catch ME % If the point doesn't exist, we've reached the beginning of the line
            if strcmp(ME.identifier, 'MATLAB:badsubscript');
                reached_beginning_flag = true;
                first_divider_coords = current_divider_coords;
                first_divider_index = current_divider_index;
            else
                throw(ME); % Rethrow any other error
            end
        end % End try 
                        
        % Calculate the distance between the divider point and the new coordinates
        current_distance = sqrt(sum((current_coords - current_divider_coords).^2, 2)); 
    
    end % End distance loop, the divider coordinate lies between the points at the current index and the previous index
    
    % 

end
end
        
function [divider_point] = divide_line_segment(segment_point_1, segment_point_2, third_point, divider_line_length)
% Local function to solve the law of sines problem in order to divide a
% line segment so that the divider length is equal to the exact value
% requested.
%
% Inputs:
%   segment_point_1/2: the two points that define the line segment along
%       the line_points curve which are known to contain the exact divider
%       endpoint.
%   third_point: the known endpoint of divider line, i.e. the
%       current_divider_coords.
%   divider_line_length: the length of the divider line.
% Outputs:
%   divider_point: the point along the line segment between points 1 and 2
%       that creates the a divider line of desired length.

% Find the angle across from the divider line
vec_12 = segment_point_2 - segment_point_1;
vec_13 = third_point - segment_point_1;
length_12 = sqrt(sum(vec_12.^2, 2));
length_13 = sqrt(sum(vec_13.^2, 2));
divider_angle = acos(dot(vec_12, vec_13) / (length_12 * length_13));

% Find the angle between the divider line and the line segment
law_of_sines_ratio = sin(divider_angle) / divider_line_length;
angle_13 = asin(length_13 * law_of_sines_ratio);

% 
    
end

