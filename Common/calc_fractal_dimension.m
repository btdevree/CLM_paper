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

% Set defaults
if nargin < 2; number_repeats = 3; end;
if nargin < 3; divider_length_range = []; end;
if nargin < 4; show_Richardson_plot = false; end;

% Set parameters and get constants
number_divider_lengths = 7;
total_number_points = size(line_points, 1);
displacement = sqrt(sum((line_points(1, :) - line_points(end-1, :)).^2, 2));
distances = sqrt(sum((line_points(2:end, :) - line_points(1:end-1, :)).^2, 2));
maximum_length = sum(distances, 1);
stdev_length = std(distances, 0, 1);
skew_length = skewness(distances, 1, 1);
mean_length = maximum_length / total_number_points;
path_to_displacement_ratio = maximum_length / displacement;

% Determine divider sizes
if isempty(divider_length_range)
   
    % Auto estimate sizes
    min_divider = .5 * mean_length; % half the mean step distance
    max_divider = .001 * maximum_length; % 1/100th of max length
else
    
    % Explicit divider sizes
    min_divider = divider_length_range(1);
    max_divider = divider_length_range(2);
end    
divider_lengths = logspace(log10(min_divider), log10(max_divider), number_divider_lengths)';

% Initalize results matrix
curve_distance = zeros(number_divider_lengths, number_repeats);

% Prepare arguments for parallel execution
args = struct('divider_length', [], 'starting_point_index', [], 'points_required_per_divider', []);
for repeat_index = 1:number_repeats
    
    % Choose a point at random to start at
    starting_point_index = randi(total_number_points);
        
    % Repeat measurement for all divider lengths
    for divider_index = 1:number_divider_lengths
        
        % Calculate the number of points that we should calculate the distance from for each cycle of the divider algorithm
        num_points = calc_number_poins_per_divider(divider_lengths(divider_index), mean_length, stdev_length, skew_length, path_to_displacement_ratio);
        
        % Record arguments
        args(divider_index, repeat_index).divider_length = divider_lengths(divider_index);
        args(divider_index, repeat_index).start_index = starting_point_index;
        args(divider_index, repeat_index).points_required_per_divider = num_points;
    end
end

% % Nonparallel execution - Run instead of parallel feval for debuging/graphing
% for eval_index = 1:size(args(:), 1)
%     curve_distance(eval_index) = calc_distance(line_points, args(eval_index), true);
% end

% Start up a pool of future calc_distance evaluations
for eval_index = 1:size(args(:), 1)
    future_results(eval_index) = parfeval(@calc_distance, 1, line_points, args(eval_index));
end

% Create an onCleanup to ensure we do not leave any futures running when we exit.
cancelFutures = onCleanup(@() cancel(future_results));

% Collect the future_results and put them in the distance result matrix
for eval_index = 1:size(args(:), 1)
   [completed_index, new_result] = fetchNext(future_results);
   curve_distance(completed_index) = new_result;
end

%  Transform to Richardson's plot
log_divider_matrix = [ones(size(curve_distance(:))), repmat(log(divider_lengths), number_repeats, 1)];
log_distance_matrix = log(curve_distance(:));

% Fit a line to the points
fit_parameters = log_divider_matrix \ log_distance_matrix;
fractal_dim = 1 - fit_parameters(2);

% Show the plot if requested
if show_Richardson_plot
    
    %colors_mat = repmat(linspace(0, 1, number_repeats), size(divider_lengths, 1), 1);
    colors_mat = vertcat(args.start_index);
    %colors = colors_mat(:);
    colors = colors_mat / max(colors_mat);
    
    % Calc line to display
    display_line_x = linspace(log(min_divider), log(max_divider));
    display_line_y = fit_parameters(2) .* display_line_x + fit_parameters(1);
    scatter(log_divider_matrix(:, 2), log_distance_matrix, 20, colors)
    hold on
    plot(display_line_x, display_line_y, 'r');
    xlabel('ln(divider length)')
    ylabel('ln(total distance)')
end
end
        
function [distance, divider_points] = calc_distance(line_points, args, draw_plots)
% Local function to calculate the distance of a curve for a particular
% divider length. Optional output gives all the divider points.

% Default is no plotting
if nargin < 3; draw_plots = false; end;

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

% Unpack args
divider_length = args.divider_length;
start_index = args.start_index;
points_required = args.points_required_per_divider; 

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
    current_divider_coords = start_coords;
    points_required_multiplier = 1;
       
    % Loop through adding dividers until we've reached the end and break the loop
    while true

        % Try to get a set of point values that is expected to contain the next divider point
        try % Usually, get the coordinates of the next/previous point and the following points in the set
             current_coord_set = line_points(current_index + step_direction : step_direction : current_index + step_direction * points_required * points_required_multiplier, :);

        % If the points don't exist, we've reached the beginning or end of the line
        catch ME 
            if strcmp(ME.identifier, 'MATLAB:badsubscript');

                % Take the set from the current_index to the end of the curve instead, current_index should always be > 1 or < end
                if step_direction == -1
                    current_coord_set = line_points(current_index + step_direction : -1 : 1, :);
                elseif step_direction == 1
                    current_coord_set = line_points(current_index + step_direction : 1 : end, :);
                end
                reached_end_flag = true;

            % Rethrow any other error
            else
                throw(ME);
            end
        end

        % Calculate the distances between the divider point and the new coordinates
        current_distances = sqrt(sum((current_coord_set - repmat(current_divider_coords, size(current_coord_set, 1), 1)).^2, 2)); 
               
        % Check if we have points that include the divider length
        crosspoint_index = find(current_distances >= divider_length, 1); % returns empty 1 x 0 matrix if there's no point that's far enough
        
        % If the first point in the coord set is already far enough away, add as many colinear dividers as possible.
        if crosspoint_index == 1
 
            % Get an offset vector for the additional lengths
            segment_vector = current_coord_set(1, :) - current_divider_coords;        
            offset_vector = segment_vector * (divider_length / current_distances(1));
 
            % Figure out how many divider lengths will fit on the remander of the line segment 
            number_additional_lengths = floor(current_distances(1) / divider_length);

            % Add any additional divider lengths on the current line segment
            for addlen_index = 1:number_additional_lengths
                current_divider_coords = current_divider_coords + offset_vector;
                step_count = step_count + 1;
                
                % Add the new point to the list, if requested
                if log_points_flag
                    if step_direction == -1 % stepping backwards
                        divider_points = [current_divider_coords; divider_points];
                    elseif step_direction == 1 % stepping forwards
                        divider_points = [divider_points; current_divider_coords];
                    end
                end
            end
            
            % Recalculate the crosspoint from the new divider coordinates
            current_distances = sqrt(sum((current_coord_set - repmat(current_divider_coords, size(current_coord_set, 1), 1)).^2, 2)); 
            crosspoint_index = find(current_distances >= divider_length, 1);
        end                
        
        % If we found a valid index, split the line segment and set new current points and indices
        if ~isempty(crosspoint_index)

            % Get the dividing point
            point_1_index = current_index + crosspoint_index * step_direction - step_direction;
            point_2_index = current_index + crosspoint_index * step_direction;
            new_divider_coords = divide_line_segment(line_points(point_1_index, :), line_points(point_2_index, :), current_divider_coords, divider_length);
            
            if isempty(new_divider_coords)
                disp('What is happening here?');
            end
            
            % Count, set indices and flags 
            current_divider_coords = new_divider_coords;
            current_index = point_1_index;
            step_count = step_count + 1;
            reached_end_flag = false; % If we found a point, we should check for another between the current index and the end
            points_required_multiplier = 1;

            % Add the new point to the list, if requested
            if log_points_flag
                if step_direction == -1 % stepping backwards
                    divider_points = [current_divider_coords; divider_points];
                elseif step_direction == 1 % stepping forwards
                    divider_points = [divider_points; current_divider_coords];
                end
            end

        % If no point in the set was long enough but we didn't get to the end of the curve yet, try again with more points
        elseif isempty(crosspoint_index) && ~reached_end_flag

            % Add another multiple of points to the set that's taken
            points_required_multiplier = points_required_multiplier + 1;

        % If no point in the set was long enough and we're at the end of the curve, add the endpoint distances and quit looking for points 
        elseif isempty(crosspoint_index) && reached_end_flag

            % Record the extra distances
            extra_distance = sqrt(sum((current_coord_set(end, 1) - current_divider_coords).^2, 2));
            if step_direction == -1 % stepping backwards
                backwards_extra = extra_distance;
            elseif step_direction == 1 % stepping forwards
                forwards_extra = extra_distance;
            end

            % Quit looking for points in this direction
            break
        end
        
        % Plot changes if requested
        if draw_plots
            draw_plot(line_points, divider_points, current_coord_set, current_divider_coords, start_index);
        end

    end
end

% Add up the distances
distance = step_count * divider_length + backwards_extra + forwards_extra;
end

function [divider_point] = divide_line_segment(segment_point_1, segment_point_2, third_point, divider_line_length)
% Local function to solve the intersection of a line and circle in order to
% divide a line segment so that the divider length is equal to the exact 
% value requested.
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
%       that creates the a divider line of desired length. Returns empty
%       matrix if there is no solution

% Solves the intersection of
% y = mx + b and (x - h)^2 + (y - k)^2 = r^2

% Find required constants - enforce conversion to double because we have a lot of adding to do
segment_point_2 = double(segment_point_2);
segment_point_1 = double(segment_point_1);
vec_12 = segment_point_2 - segment_point_1;
m = vec_12(2) / vec_12(1); % m = slope of line
b = segment_point_1(2) - m * segment_point_1(1); % b = intercept of line
h = double(third_point(1));
k = double(third_point(2));
r = double(divider_line_length);

% Normal solving, but a special case is required if the line is vertical
if ~isinf(m)

    % Subsitute y = mx + b into circle equation, expanding all products and collecting terms gives:
    % x^2(1 + m^2) + x(-2h + 2mb - 2mk) + (h^2 + b^2 - 2bk + k^2 - r^2) = 0
    % Problem is quadratic, solve for two solutions of x, then find y with line equation
    A = 1 + m^2;
    B = 2 * (-h + m * b - m * k);
    C = h^2 + b^2 - 2 * b * k + k^2 - r^2;
    discriminant = B^2 - 4 * A * C;

    % Results based on discriminant
    % Two intersection points
    if discriminant > 0
        x_1 =  (-B - sqrt(discriminant)) / (2 * A);
        x_2 =  (-B + sqrt(discriminant)) / (2 * A);
        y_1 = m * x_1 + b;
        y_2 = m * x_2 + b;
        
        % Choose the correct point
        [divider_point] = choose_point([x_1, y_1], [x_2, y_2], segment_point_1, segment_point_2);

    % Line is tangent
    elseif discriminant == 0 
        x = -B / (2 * A);
        y = m * x + b;
        divider_point = [x, y];

    % No intersection
    elseif discriminant < 0
        divider_point = [];
    end

% Special case when line is vertical, x is a constant, so solve circle equation for y as a quadratic equation: y^2 - 2ky + k^2 + (x - h)^2 -r^2 = 0
else
    x = double(segment_point_1(1));
    A = 1;
    B = -2 * k;
    C = k^2 + (x - h)^2 - r^2;
    discriminant = B^2 - 4 * A * C;

    % Results based on discriminant
    % Two intersection points
    if discriminant > 0
        y_1 =  (-B - sqrt(discriminant)) / (2 * A);
        y_2 =  (-B + sqrt(discriminant)) / (2 * A);
        
        % Choose the correct point
        [divider_point] = choose_point([x, y_1], [x, y_2], segment_point_1, segment_point_2);
        
    % Line is tangent
    elseif discriminant == 0 
        y = -B / (2 * A);
        divider_point = [x, y];

    % No intersection
    elseif discriminant < 0
        divider_point = [];
    end
end

% Return divider point as a single if point 1 is also single
if isa(segment_point_1, 'single')
    divider_point = single(divider_point);
end
end

function [point] = choose_point(p1, p2, segment_point_1, segment_point_2)
% Local function to choose the correct point when given two solutions

% Test which solution is for a point that lies between points 1 and 2 
min_x = min(segment_point_1(1), segment_point_2(1));
max_x = max(segment_point_1(1), segment_point_2(1));
min_y = min(segment_point_1(2), segment_point_2(2));
max_y = max(segment_point_1(2), segment_point_2(2));
tol_x = 2*eps(single(max_x)); % Allow for precision issues with single-point coordinates 
tol_y = 2*eps(single(max_y));
if min_x - p1(1) <= tol_x && p1(1) - max_x <= tol_x && min_y - p1(2) <= tol_y && p1(2) - max_y <= tol_y
    p1_OK_flag = true;
else
    p1_OK_flag = false;
end
if min_x - p2(1) <= tol_x && p2(1) - max_x <= tol_x && min_y - p2(2) <= tol_y && p2(2) - max_y <= tol_y
    p2_OK_flag = true;
else
    p2_OK_flag = false;
end

% Usually only one solution will be valid, but it is possible that both are OK. In this case we use the solution that is closer to 
% point 1 as this is the first point that the divider line crosses as we trace along line_points.
if p1_OK_flag && ~p2_OK_flag
    point = p1;
elseif ~p1_OK_flag && p2_OK_flag
    point = p2;
elseif p1_OK_flag && p2_OK_flag;
    dist_p1 = sqrt(sum((segment_point_1 - p1).^2, 2));
    dist_p2 = sqrt(sum((segment_point_1 - p2).^2, 2));
    if dist_p1 <= dist_p2
        point = p1;
    else
        point = p2;
    end
end
end

function [number_points] = calc_number_poins_per_divider(divider, mean, stdev, skew, path_length_to_displacement_ratio)
% Local function to estimate the number of points we'll need to assay in
% order to include the next divider length 99.5% of the time. Accounts for
% skewed step size distributions and small numbers of steps. 

% Loop cap
max_iterations = 1e4;

% Increase the number of points until we find enough to satisfy the problem
number_points = 1;
while true 
    
    % Determine the statistics of the path length sum
    path_mean = number_points * mean; % Wasserman, "All of Statistics" page 28, iii.
    path_stdev = sqrt(number_points) * stdev; % Wasserman, "All of Statistics" page 28, iii.
    path_skew = (1/sqrt(number_points)) * skew; % From skew definition, given in Eriksson, "A simulation method for skewness correction" (Master's thesis, Uppasla U.), Appendex 1, Corollary 4. 
    
    % Approximate as normal with a corrected mean due to the skew
    adjusted_path_mean = path_mean - (path_skew / (6 * path_stdev.^2 * number_points)); % Norman J. Johnson, "Modified t Tests and Confidence Intervals for Asymmetrical Populations" JASA, v73:p536-44, eqn 2.7, but I swear the plus sign should be a minus

    % Get value of the path length that is smaller than 99.5% of the expected paths
    cutoff_path_length = adjusted_path_mean - 2.78 * path_stdev; % p = 0.005 for single tailed test
    
%     % Assume path is 2D random walk: rmsd = <r> * sqrt(n) --> (path length / n) * sqrt(n) --> path length / sqrt(n)
%     % NOTE: behavior as n --> infinity: rmsd = sqrt(n)* mu - 2.78 * sigma, we shouldn't get stuck in an infinite loop, but this is a foolish
%     %   algorithm for large dividers and small steps. In such a case, use the limiting behavior equation instead. Asks for far too many
%     %   points for this application, we have generally directed motion.
%     adjusted_cutoff_estimate = cutoff_path_length / sqrt(number_points);

    % Adjust path length by a simple ratio of the overall path length to displacement.
    adjusted_cutoff_estimate = cutoff_path_length / path_length_to_displacement_ratio;
    
    % Stop loop if we've gone enough steps or reached the cap
    if adjusted_cutoff_estimate > divider || number_points >= max_iterations;
        break
    else
        % Add another point
        number_points = number_points + 1;
    end
end

end

function draw_plot(line_points, divider_points, current_point_set, current_divider_point, start_index)
% function for aiding in visualization of progress during troubleshooting
plot(gca, line_points(max(1, start_index - 100):start_index, 1), line_points(max(1, start_index - 100):start_index, 2), '-or');
hold on
plot(divider_points(:, 1), divider_points(:, 2), '-xb');
scatter(current_point_set(:, 1), current_point_set(:, 2), 5, 'g');
scatter(current_divider_point(1), current_divider_point(2), 5, 'k')
axis equal
hold off
waitforbuttonpress
end
