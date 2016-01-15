% Script to collect correlations from 1D "dots" 

% Set directory information
figure_path_parts = strsplit(pwd, 'CLM_paper');
figure_path = [figure_path_parts{1}, 'CLM_figures_and_data/'];
data_name = 'stdev_vs_num_dots_no_over.mat';

% Load a parameter structure from current dirctory with deisred settings
load('dots_params.mat') % loads as 'params'

% Set parameters
number_replicates = 30;
number_dots_vector = [10, 20, 30, 50, 70, 100, 150, 200, 300, 400, 500];
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
    mean_vector_stack = zeros(max_radius_px + 1, number_replicates);
    parfor replicate_index = 1:number_replicates

        % Get a results stack
        [correlation_stack, number_cells] = collect_dot_correlations(params, number_dots, max_correlation_radius, dots_per_cell_mean,...
            dot_radius, dot_correlation_value, label_density_mean, label_density_stdev, label_SN_ratio, event_overcounting,...
            dot_center_precision, event_precision, STORM_pixel_resolution, 'pdf');
    
        % Average the stack
        correlation_image = mean(correlation_stack, 3);

        % Get a radial average
        [distance_vector, mean_vector, ~, ~] = radial_average_2D_correlation(correlation_image);

        % Add to the results stack
        mean_vector_stack(:, replicate_index) = mean_vector;
    end
    
    % Store the mean vector stack in a cell
    mean_vector_stack_cells{number_dot_index} = mean_vector_stack;
    fprintf('finished number_dots = %d\n', number_dots);
end

% Save data
save([figure_path, data_name], 'distance_vector', 'mean_vector_stack_cells', 'number_dots_vector', 'max_correlation_radius', 'dots_per_cell_mean', 'dot_radius',...
    'dot_correlation_value', 'label_density_mean', 'label_density_stdev', 'label_SN_ratio', 'event_overcounting', 'dot_center_precision',...
    'event_precision', 'STORM_pixel_resolution');


