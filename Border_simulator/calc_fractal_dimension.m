function [fractal_dim] = calc_fractal_dimension(line_points, number_repeats, divider_length_range, show_Richardson_plot)
%CALC_FRACTAL_DIMENSION Estimates the fractal dimension of a linear series
% of (x, y) points.
%
% Estimates the fractal dimension with a hand and divider algorithm.
%
% Input:
%   line_points: n by 2 matrix of (x, y) points.
%   number_repeats: number of times to repeat the calculation on different
%       starting points. Optional, default = 30.
%   divider_length_range: 2 by 1 matrix [min_divider, max_divider].
%       Optional, default = [], automatic range choice.
%   show_Richardson_plot: create a figure that shows the log(divider) vs. 
%       log(distance) plot and the fitted line. Optional, default = false.
% Output:
%   fractal_dim: estimate of the fractal dimension.

% NOTE: Not the greatest algorithm; it was suprisingly complex to code. I
% suspect either vectors or intersections of line segments and circles can
% be used simplify it, if the code is being refactored. A plotting function
% in the calc_distance local function can be turned on to aid in
% troubleshooting.

% Set defaults
if nargin < 2; number_repeats = 30; end;
if nargin < 3; show_Richardson_plot = false; end;

% Set parameters
number_divider_lengths = 10;

% Determine divider sizes if not given
if isempty(divider_length_range)
   
    % Auto estimate sizes
    total_number_points = size(line_points, 1);
    maximum_length = sum(sqrt(sum((line_points(2:end, :) - line_points(1:end-1, :)).^2, 2)), 1);
    min_divider = .05 * maximum_length / total_number_points; % half the mean step distance
    max_divider = .2 * maximum_length; % 1/5th of max length
    divider_lengths = logspace(log10(min_divider), log10(max_divider), number_divider_lengths)';
else
     min_divider = divider_length_range(1);
     max_divider = divider_length_range(2);
end    

% Initalize distance results matrix
curve_distance = zeros(number_divider_lengths, number_repeats);

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
log_divider_matrix = [ones(size(curve_distance(:))), repmat(log(divider_lengths), number_repeats, 1)];
log_distance_matrix = log(curve_distance(:));

% Fit a line to the points
fit_parameters = log_divider_matrix \ log_distance_matrix;
fractal_dim = 1 - fit_parameters(2);

% Show the plot if requested
if show_Richardson_plot
    
    % Calc line to display
    display_line_x = linspace(log(min_divider), log(max_divider));
    display_line_y = fit_parameters(2) .* display_line_x + fit_parameters(1);
    scatter(log_divider_matrix(:, 2), log_distance_matrix)
    hold on
    plot(display_line_x, display_line_y, 'r');
    xlabel('ln(divider length)')
    ylabel('ln(total distance)')
end
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
                draw_plot(line_points, divider_points, current_coords, current_divider_coords, start_index);
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
            draw_plot(line_points, divider_points, current_coords, current_divider_coords, start_index);
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
                draw_plot(line_points, divider_points, current_coords, current_divider_coords, start_index);
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
length_12 = sqrt(sum(vec_12.^2, 2));
length_13 = sqrt(sum(vec_13.^2, 2));
angle_213 = acos(dot(vec_12, vec_13) / (length_12 * length_13));

% Special case if points 1 and 3 are the same
if length_13 == 0;
    % Add new point on the line segment between 1 and 2
    divider_point = segment_point_1 + vec_12 .* (divider_line_length / length_12);
else

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
    if isinf(slope_12) % Special solutions for vertical lines
        divider_x_s1 = segment_point_1(1);
        divider_x_s2 = segment_point_1(1);
        divider_y_s1 = slope_divider_s1 * divider_x_s1 + intercept_divider_s1;
        divider_y_s2 = slope_divider_s2 * divider_x_s2 + intercept_divider_s2;
    elseif isinf(slope_divider_s1) && ~isinf(slope_divider_s2)
        divider_x_s1 = third_point(1);
        divider_x_s2 = (intercept_divider_s2 - intercept_12) / (slope_12 - slope_divider_s2);
        divider_y_s1 = slope_12 * divider_x_s1 + intercept_12;
        divider_y_s2 = slope_12 * divider_x_s2 + intercept_12;
    elseif ~isinf(slope_divider_s1) && isinf(slope_divider_s2)
        divider_x_s1 = (intercept_divider_s1 - intercept_12) / (slope_12 - slope_divider_s1);
        divider_x_s2 = third_point(1);
        divider_y_s1 = slope_12 * divider_x_s1 + intercept_12;
        divider_y_s2 = slope_12 * divider_x_s2 + intercept_12;
    elseif isinf(slope_divider_s1) && isinf(slope_divider_s2) % Pretty unlikley, but just in case...
        divider_x_s1 = third_point(1);
        divider_x_s2 = third_point(1);
        divider_y_s1 = slope_12 * divider_x_s1 + intercept_12;
        divider_y_s2 = slope_12 * divider_x_s2 + intercept_12;
    else
        divider_x_s1 = (intercept_divider_s1 - intercept_12) / (slope_12 - slope_divider_s1);
        divider_x_s2 = (intercept_divider_s2 - intercept_12) / (slope_12 - slope_divider_s2);
        divider_y_s1 = slope_12 * divider_x_s1 + intercept_12;
        divider_y_s2 = slope_12 * divider_x_s2 + intercept_12;
    end

    % Test if the solution is for a point that lies between points 1 and 2 
    min_x = min(segment_point_1(1), segment_point_2(1));
    max_x = max(segment_point_1(1), segment_point_2(1));
    min_y = min(segment_point_1(2), segment_point_2(2));
    max_y = max(segment_point_1(2), segment_point_2(2));
    if min_x <= divider_x_s1 && divider_x_s1 <= max_x && min_y <= divider_y_s1 && divider_y_s1 <= max_y
        s1_OK_flag = true;
    else
        s1_OK_flag = false;
    end
    if min_x <= divider_x_s2 && divider_x_s2 <= max_x && min_y <= divider_y_s2 && divider_y_s2 <= max_y
        s2_OK_flag = true;
    else
        s2_OK_flag = false;
    end

    % Usually only one solution will be valid. It is possible that both are OK if the 1d3 angle
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
end % end else for special case with points 1 and 3
end

function draw_plot(line_points, divider_points, current_point, current_divider_point, start_index)
% function for aiding in visualization of progress during troubleshooting
plot(gca, line_points(max(1, start_index - 100):start_index, 1), line_points(max(1, start_index - 100):start_index, 2), '-or');
hold on
plot(divider_points(:, 1), divider_points(:, 2), '-xb');
scatter(current_point(1), current_point(2), 5, 'g');
scatter(current_divider_point(1), current_divider_point(2), 5, 'k')
hold off
waitforbuttonpress
end
