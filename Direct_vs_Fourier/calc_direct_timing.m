function [timing_vector, whos_struct_cell_vector] = calc_direct_timing(number_points_vector, method, verbose_flag)
%CALC_DIRECT_TIMING Calculate the run time and memory use for a direct 
%   calculation method of autocorrelations 
%
%   Calculates a radial pairwise autocorrelation curve from points
%   uniformaly distributed in a unit square and records the time of
%   execution and the memory use for the calculation.
%   Input:
%   number_points_vector: Column vector of doubles; the number of points to 
%       use in the autocorrelation
%   method: string denoting the method used for calculation. Choices are
%       'pdist2', 'custom_memory', and 'custom_func_calls'.
%   verbose_flag: logical value, default = false. Set to true so that the
%       function report progress to the console.
%   Output:
%   timing_vector: Time of execution for calculating the correlation and
%       the radial average. Best of 10 repetitions for number_of_points <= 
%       1e4, best of 3 repetitions for number_of_points <= 1e6, and a 
%       single repetition for all else. Column vector of doubles.
%   whos_struct_cell_vector: output of the whos command at the end of each
%       requested number of points. Column vector of cells containing 
%       strucures.

% Set default values
if nargin < 3; verbose_flag = false; end; 

% Starting message
if verbose_flag
    fprintf(['Start direct with method ', method, ': \n']);
end

% Define repetition limits
repeat_limits = [1e4, 5e5];

% Define radial resolution and max length
radial_resolution = 2e-4; % 1/5000
number_histogram_bins = 143; % 1um if each bin = 7nm

% Loop through each number_of_points value
for number_points_index = 1:size(number_points_vector, 1);
    number_points = number_points_vector(number_points_index);

    % Create list of points
    points = rand(number_points, 2);
    
    % Determine the number of repeats
    if number_points <= repeat_limits(1)
        number_repeats = 5;
    elseif number_points <= repeat_limits(2)
        number_repeats = 4; 
    else
        number_repeats = 3;
    end
    
    % Repeat the specified number of times
    for repeat_index = 1:number_repeats
        
        % Start timer
        tic
        
        % Calculate distances with pdist2
        if strcmp(method, 'pdist2')
            
            % Calculate the distance
            dist_matrix = pdist2(points, points);
            
            % Take distance as a vector and exclude the distance of a point to itself
            dist_vector = dist_matrix(~diag(true(1, length(dist_matrix))));
        
            
        % Calculate distances with custom algorthim that minimzes memory use    
        elseif strcmp(method, 'custom_memory')
            
            % Initialize vector of distances
            dist_vector = zeros(nchoosek(number_points, 2), 1);
            
            % Build distance on subset of measurements at a time
            last_filled_dist_index = 0;
            for subset_index = 1:number_points - 1
                
                % Subtract the subset that start fromthe first point from the subset that starts at the last point
                xy_dist = points(1:end - subset_index, :) - points(1 + subset_index:end, :);
                
                % Square the x and y distances
                xy_squared = xy_dist.^2;
                
                % Add columns and take the square root
                dist_subset = sqrt(xy_squared(:, 1) + xy_squared(:, 2));
                
                % Add to distance vector
                start_index = last_filled_dist_index + 1;
                end_index = last_filled_dist_index + size(dist_subset, 1);
                dist_vector(start_index:end_index, :) = dist_subset;
                
                % Update last_filled_dist_index
                last_filled_dist_index = end_index;
            end
           
        % Calculates distances with custom method that minimizes the number of function calls
        elseif strcmp(method, 'custom_func_calls')
            
            % Create arrays of the x and y dimensions of the points
            point_array_x = repmat(points(:, 1), [1, number_points]);
            point_array_y = repmat(points(:, 2), [1, number_points]);
            
            % Subtract array and its transposes
            x_dist_array = point_array_x - point_array_x.';
            y_dist_array = point_array_y - point_array_y.';
            
            % Calculate the distance
            dist_matrix = hypot(x_dist_array, y_dist_array);
            
            % Take distance as a vector and exclude the distance of a point to itself
            dist_vector = dist_matrix(~diag(true(1, length(dist_matrix))));
        end
        
        % Calculate the radial autocorrelation
        edges = [0:number_histogram_bins + 1] .* radial_resolution;
        counts = histc(dist_vector, edges);
        
        % Double the counts if we're using the memory effecient method
        if strcmp(method, 'custom_memory')
            counts = 2 * counts;
        end
        
        % We're done, get the timing
        timing_repeats(repeat_index) = toc;
        
        % Report the repeat completion
        if verbose_flag
            fprintf('.');
        end
    end
    
    % Take the best timing
    timing_vector(number_points_index) = min(timing_repeats);
    
    % Report the completion
    if verbose_flag
        fprintf([num2str(number_points), ' points: best time = ', num2str(min(timing_repeats)), ' seconds \n']);
    end
    
    % Record the whos vector
    whos_struct_cell_vector{number_points_index} = whos;
    
    % Clear the repeat vector
    clear timing_repeats
end

% Convert outputs to column vectors
timing_vector = timing_vector.';
whos_struct_cell_vector = whos_struct_cell_vector.';
end
   

