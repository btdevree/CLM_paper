function [fractal_dim] = calc_fractal_dimension(line_points, number_repeats)
%CALC_FRACTAL_DIMENSION Estimates the fractal dimension of a linear series
% of (x, y) points.
%
% Estimates the fractal dimension with a hand and divider algorithm.
%
% Input:
%   line_points: n by 2 matrix of (x, y) points.
%   number_repeats: number of times to repeat the calculation on different
%       starting points. Optional, default = 30.
% Output:
%   fractal_dim: estimate of the fractal dimension.

% Set defaults
if nargin < 2; number_repeats = 30; end;

% Set parameters
number_divider_lengths = 10;

% Determine divider sizes
total_number_points = size(line_points, 1);
maximum_length = sum(sqrt(sum((line_points(2:end, :) - line_points(1:end-1, :)).^2, 2)), 1);
min_divider = 0.5 * maximum_length / total_number_points; % Half the mean step distance
max_divider = 0.1 * maximum_length; % 1/10th of max length
divider_lengths = logspace(log10(min_divider), log10(max_divider), number_divider_lengths)';

% Initalize distance results matrix
curve_distance = zeros(number_divider_lengthes, number_repeats);

% Repeat the measurement process from multiple starting points
for repeat_index = 1:number_repeats
    
    % Choose a point at random to start at
    starting_point_index = randi(total_number_points);
        
    % Repeat measurement for all divider lengths
    for divider_index = 1:number_divider_lengths
        divider_length = divider_lengths(divider_index);
        curve_distance(divider_index, repeat_index) = calc_distance(line_points, divider_length, starting_point_index);
    end
end

%  Transform to Richardson's plot
log_divider_matrix = [ones(size(curve_distance(:))), repmat(log(divider_lengths), number_repeats)];
log_distance_matrix = curve_distance(:);

% Fit a line to the points
fit_parameters = log_divider_matrix \ log_distance_matrix;
fractal_dim = 1 - fit_parameters(2);
end
        
function [distance, divider_points] = calc_distance(line_points, divider_length, start_index)
% Local function to calculate the distance of a curve for a particular
% divider length. Optional output gives all the divider points.

% Enable/disable graphical feedback, false = off
draw_plots = false;

% See if we'll need to log points
if nargout > 1
    log_points_flag = true;
else
    log_points_flag = false;
end

% Create new figure and turn on point logging if feedback is on
if draw_plots
    figure('Units', 'pixels', 'Position', [50, 50, 1000, 1000], 'PaperPositionMode', 'auto', 'InvertHardcopy', 'off');
    log_points_flag = true;
end

% Get starting coordinates
start_coords = line_points(start_index, :);
if log_points_flag
    divider_points = start_coords;
end

% Walk the divider lengths forwards and backwards from the start point
for step_direction = [-1, 1]

    %Initalize flags, counters, and coordinates
    step_count = 0;
    reached_end_flag = false;
    current_index = start_index;
    current_coords = start_coords;
    current_divider_coords = start_coords;
    current_distance = 0;
       
    % Loop through divider lengths until the line points are used up
    while true

        % Go through line points until the step length is larger than the divider
        while current_distance < divider_length

            % Try to get values at the previous index
            try % Usually, get the coordinates of the previous point on the line
                current_coords = line_points(current_index + step_direction, :);
                current_index = current_index + step_direction; % Will not evaluate if error is thrown
            catch ME % If the point doesn't exist, we've reached the beginning or end of the line
                if strcmp(ME.identifier, 'MATLAB:badsubscript');
                    reached_end_flag = true;
                    break % exit line point walking loop
                else
                    throw(ME); % Rethrow any other error
                end
            end % End try 
            
            % Calculate the distance between the divider point and the new coordinates
            current_distance = sqrt(sum((current_coords - current_divider_coords).^2, 2)); 
            if draw_plots
                draw_plot(line_points, divider_points, current_coords, current_divider_coords);
            end
            
        end % End distance loop, the divider coordinate lies between the points at the current index and the previous index

        % Quit looping when we hit the end of the line points
        if reached_end_flag
            break % Exit divider finding loop
        end

        % Get the dividing point
        new_divider_coords = divide_line_segment(line_points(current_index - step_direction, :), line_points(current_index, :), current_divider_coords, divider_length);

        % Add the new point to the list, if requested
        if log_points_flag
            if step_direction == -1 % stepping backwards
                divider_points = [new_divider_coords; divider_points];
            elseif step_direction == 1 % stepping forwards
                divider_points = [divider_points; new_divider_coords];
            end
        end

        % Set the new divider coords and increment counter 
        current_divider_coords = new_divider_coords;
        step_count = step_count + 1;

        % Recalculate current step distance
        current_distance = sqrt(sum((current_coords - current_divider_coords).^2, 2));
        if draw_plots
            draw_plot(line_points, divider_points, current_coords, current_divider_coords);
        end

        % Check if the next divider falls on the current line segment; if so,
        %   add divider points until less than one divider length remains.
        if current_distance > divider_length

            % Get an offest vector for the additional lengths
            segment_vector = current_coords - current_divider_coords;        
            offset_vector = segment_vector * (divider_length / current_distance);

            % Figure out how many divider lengths will fit on the remander of the line segment 
            number_additional_lengths = floor(current_distance / divider_length);

            % Add any additional divider lengths on the current line segment
            for addlen_index = 1:number_additional_lengths
                new_divider_coords = current_divider_coords + offset_vector;

                % Add the new point to the list, if requested
                if log_points_flag
                    if step_direction == -1 % stepping backwards
                        divider_points = [new_divider_coords; divider_points];
                    elseif step_direction == 1 % stepping forwards
                        divider_points = [divider_points; new_divider_coords];
                    end
                end

                % Set the new divider coords and increment counter 
                current_divider_coords = new_divider_coords;
                step_count = step_count + 1;
            end

            % Recalculate current step distance
            current_distance = sqrt(sum((current_coords - current_divider_coords).^2, 2));
            if draw_plots
                draw_plot(line_points, divider_points, current_coords, current_divider_coords);
            end
            
        end % end colinear segment adding  
    end % end divider function finding

    % Get the distance between the last divider coordinate and the last point of the line
    extra_distance = current_distance;
    
    % Add the new point to the list, if requested
    if log_points_flag
        if step_direction == -1 % stepping backwards
            divider_points = [line_points(current_index, :); divider_points];
        elseif step_direction == 1 % stepping forwards
            divider_points = [divider_points; line_points(current_index, :)];
        end
    end
    
    % Add up distance
    direction_distance = step_count * divider_length + extra_distance;
    
    % Record the distances
    if step_direction == -1 % stepping backwards
        backwards_distance = direction_distance;
    elseif step_direction == 1 % stepping forwards
        forwards_distance = direction_distance;
    end
    
end % End forward or backward for loop

% Add backwards and forwards distances together
distance = backwards_distance + forwards_distance;
end

function [divider_point] = divide_line_segment(segment_point_1, segment_point_2, third_point, divider_line_length)
% Local function to solve the law of sines problem in order to divide a
% line segment so that the divider length is equal to the exact value
% requested.
%
% Inputs:
%   segment_point_1/2: the two points that define the line segment along
%       the line_points curve which are known to contain the exact divider
%       endpoint. Point 1 should be the one closer to the third point. 
%       Given as a 1 by 2 (x, y) matrix.
%   third_point: the known endpoint of divider line, i.e. the
%       current_divider_coords. Given as a 1 by 2 (x, y) matrix.
%   divider_line_length: the length of the divider line.
% Outputs:
%   divider_point: the point along the line segment between points 1 and 2
%       that creates the a divider line of desired length.

% Find the angle across from the divider line
vec_12 = segment_point_2 - segment_point_1;
vec_13 = third_point - segment_point_1;
vec_23 = third_point - segment_point_2;
length_12 = sqrt(sum(vec_12.^2, 2));
length_13 = sqrt(sum(vec_13.^2, 2));
length_23 = sqrt(sum(vec_23.^2, 2));
angle_213 = acos(dot(vec_12, vec_13) / (length_12 * length_13));
angle_123 = acos(dot(-vec_12, vec_23) / (length_12 * length_23));
angle_132 = acos(dot(-vec_13, -vec_23) / (length_13 * length_23));

% Find the two solutions to the angle between the divider line and the line segment.
law_of_sines_ratio = sin(angle_213) / divider_line_length;
angle_1d3_s1 = asin(length_13 * law_of_sines_ratio); % asin always returns angles between -pi/2 and pi/2
angle_1d3_s2 = pi - angle_1d3_s1;

% Get the equations of the two intersecting lines
slope_12 = vec_12(2) / vec_12(1);
intercept_12 = segment_point_1(2) - slope_12 * segment_point_1(1);
slope_divider_s1 = tan(angle_1d3_s1 + atan(slope_12) - pi);
slope_divider_s2 = tan(angle_1d3_s2 + atan(slope_12) - pi);
intercept_divider_s1 = third_point(2) - slope_divider_s1 * third_point(1);
intercept_divider_s2 = third_point(2) - slope_divider_s2 * third_point(1);

% Solve the equations
divider_x_s1 = (intercept_divider_s1 - intercept_12) / (slope_12 - slope_divider_s1);
divider_x_s2 = (intercept_divider_s2 - intercept_12) / (slope_12 - slope_divider_s2);
divider_y_s1 = slope_12 * divider_x_s1 + intercept_12;
divider_y_s2 = slope_12 * divider_x_s2 + intercept_12;

% Test if the solution is for a point that lies between points 1 and 2 
min_x = min(segment_point_1(1), segment_point_2(1));
max_x = max(segment_point_1(1), segment_point_2(1));
min_y = min(segment_point_1(2), segment_point_2(2));
max_y = max(segment_point_1(2), segment_point_2(2));
if min_x <= divider_x_s1 <= max_x && min_y <= divider_y_s1 <= max_y
    s1_OK_flag = true;
else
    s1_OK_flag = false;
end
if min_x <= divider_x_s2 <= max_x && min_y <= divider_y_s2 <= max_y
    s2_OK_flag = true;
else
    s2_OK_flag = false;
end

% Usually on one solution will be valid. It is possible that both are OK if the 1d3 angle
% is close to 90 degrees, in which case we use solution that is closer to the closer to 
% point 1 and thus the first point that the divider line crosses as we trace along line_points.
if s1_OK_flag && ~s2_OK_flag
    divider_point = [divider_x_s1, divider_y_s1];
elseif ~s1_OK_flag && s2_OK_flag
    divider_point = [divider_x_s2, divider_y_s2];
elseif s1_OK_flag && s2_OK_flag;
    dist_s1 = sqrt(sum((segment_point_1 - [divider_x_s1, divider_y_s1]).^2, 2));
    dist_s2 = sqrt(sum((segment_point_1 - [divider_x_s2, divider_y_s2]).^2, 2));
    if dist_s1 <= dist_s2
        divider_point = [divider_x_s1, divider_y_s1];
    else
        divider_point = [divider_x_s2, divider_y_s2];
    end
end
end

function draw_plot(line_points, divider_points, current_point, current_divider_point)
plot(gca, line_points(1:100, 1), line_points(1:100, 2), '-or');
hold on
plot(divider_points(:, 1), divider_points(:, 2), '-xb');
scatter(current_point(1), current_point(2), 5, 'g');
scatter(current_divider_point(1), current_divider_point(2), 5, 'k')
hold off
waitforbuttonpress
end

