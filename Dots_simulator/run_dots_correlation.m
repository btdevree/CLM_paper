% Script to collect correlations from 1D "dots" 

% Set directory information
figure_path_parts = strsplit(pwd, 'CLM_paper');
figure_path = [figure_path_parts{1}, 'CLM_figures_and_data/'];
data_name = 'stdev_vs_num_dots_no_over.mat';
figure_name = 'stdev_vs_num_dots.png';

% Load a parameter structure from current dirctory with deisred settings
load('dots_params.mat') % loads as 'params'

% Set parameters
number_replicates = 5;
number_bootstraps = 5;
number_dots_vector = [10];%, 20, 30, 50, 70, 100, 150, 200, 300, 400, 500];
max_correlation_radius = 1000; %nm
dots_per_cell_mean = 25;
dot_radius = 100; %nm
dot_correlation_value = 5;
label_density_mean = 30; %events/um^2
label_density_stdev = 5; %events/um^2
label_SN_ratio = 'no_noise';
event_overcounting = 0;
dot_center_precision = 25; %nm
event_precision = 25; %nm
STORM_pixel_resolution = 7; %nm

mean_vector_stack_cells = cell(size(number_dots_vector));
% Repeat for each number of dots
for number_dot_index = 1:length(number_dots_vector)
    number_dots = number_dots_vector(number_dot_index);
    
    % Calculate the size of the correlation stack
    max_radius_px = ceil(max_correlation_radius / STORM_pixel_resolution);
    correlation_width = 2*max_radius_px + 1;

    % Repeat dot collection and averageing 
    for index = 1:number_replicates

        % Evaluate the function asynchronously
        future_results(index) = parfeval(@eval_true_mean, 2, params, number_dots, max_correlation_radius, dots_per_cell_mean,...
            dot_radius, dot_correlation_value, label_density_mean, label_density_stdev, label_SN_ratio, event_overcounting,...
            dot_center_precision, event_precision, STORM_pixel_resolution, 'pdf'); 
    end
    
    % Collect results
    mean_vector_stack = zeros(max_radius_px + 1, number_replicates);
    for index = 1:number_replicates

        % Get next available result, fetchNext blocks until next results are available.
        [completed_index, distance_vector, mean_vector] = fetchNext(future_results);
      
        % Save in appropreate slice
        mean_vector_stack(:, completed_index) = mean_vector;
    end
    
    % Store the mean vector stack in a cell
    mean_vector_stack_cells{number_dot_index} = mean_vector_stack;
    fprintf('finished true mean number_dots = %d\n', number_dots);
end

mean_bootstrap_vector_stack_cells = cell(number_bootstraps, size(number_dots_vector, 2));
% Repeat for each number of dots
for number_dot_index = 1:length(number_dots_vector)
    number_dots = number_dots_vector(number_dot_index);
    
    % Calculate the size of the correlation stack
    max_radius_px = ceil(max_correlation_radius / STORM_pixel_resolution);
    correlation_width = 2*max_radius_px + 1;
    
    % Repeat entire bootstrapping process in parallel
    for index = 1:number_bootstraps

        % Evaluate the function asynchronously
        future_results(index) = parfeval(@eval_bootstrap_mean, 2, number_replicates, params, number_dots, max_correlation_radius, dots_per_cell_mean,...
            dot_radius, dot_correlation_value, label_density_mean, label_density_stdev, label_SN_ratio, event_overcounting,...
            dot_center_precision, event_precision, STORM_pixel_resolution, 'pdf'); 
    end
    
    % Collect results
    mean_vector_stack = zeros(max_radius_px + 1, number_replicates);
    for index = 1:number_replicates

        % Get next available result, fetchNext blocks until next results are available.
        [completed_index, distance_vector, mean_vector] = fetchNext(future_results);
      
        % Save in appropreate slice
        mean_vector_stack(:, completed_index) = mean_vector;
    end
    % Repeat entire bootstrapping process
    for bootstrap_index = 1:number_bootstraps
        
        % Get a results stack
        [correlation_stack, number_cells] = collect_dot_correlations(params, number_dots, max_correlation_radius, dots_per_cell_mean,...
            dot_radius, dot_correlation_value, label_density_mean, label_density_stdev, label_SN_ratio, event_overcounting,...
            dot_center_precision, event_precision, STORM_pixel_resolution, 'pdf');

        % Repeat getting the mean with resampling
        mean_vector_stack = zeros(max_radius_px + 1, number_replicates);
        for replicate_index = 1:number_replicates

            % Permute the stack
            permuted_indices = randi(size(correlation_stack, 3));
            permuted_stack = cat(3, correlation_stack(:, :, permuted_indices));
            permuted_image = mean(permuted_stack, 3);

            % Get a radial average
            [distance_vector, mean_vector, ~, ~] = radial_average_2D_correlation(permuted_image);

            % Add to the results stack
            mean_vector_stack(:, replicate_index) = mean_vector;
        end
        
        % Store the mean vector stack in a cell
        mean_bootstrap_vector_stack_cells{bootstrap_index, number_dot_index} = mean_vector_stack;
    end
    fprintf('finished bootstrap number_dots = %d\n', number_dots);
end


% Save data
save([figure_path, data_name], 'distance_vector', 'mean_vector_stack_cells', 'number_dots_vector', 'max_correlation_radius', 'dots_per_cell_mean', 'dot_radius',...
    'dot_correlation_value', 'label_density_mean', 'label_density_stdev', 'label_SN_ratio', 'event_overcounting', 'dot_center_precision',...
    'event_precision', 'STORM_pixel_resolution', 'mean_bootstrap_vector_stack_cells');


