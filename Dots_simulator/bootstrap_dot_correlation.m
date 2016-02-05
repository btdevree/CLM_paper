function [correlation_stack, number_cells_screened] = bootstrap_dot_correlation(parameter_struct, number_resamplings, number_dots, max_corr_length,...
    dots_per_cell_mean, dot_radius, dot_correlation_value, label_density_mean, label_density_stdev, label_SN_ratio, event_overcounting,...
    dot_center_precision, event_precision, STORM_pixel_resolution, STORM_method)
%BOOTSTRAP_DOT_CORRELATIONS Returns 2D correlations from simulated cells 
%   with "dot" type well-defined centers of interest
%
%   This function handles the logic needed to run a dot simulation with the
%   given parameters. Strategy is a bit hacky, need to think about a
%   refactor once it's working. We create the specified number of labeled
%   spots as channel 1, and then sample them one (or more) times in channel
%   2. We only use channel 2 in the end. 
%   
% Inputs:
%   param_struct: structure that holds the parameters needed for simulating
%       a cell with the create_test_data/create_test_STORM_image
%       functions. Parameters in this structure may be superceded due to 
%       requirements of the other inputs.
%   number_resamplings: number of bootstrapping resamplings to perform
%   number_dots: the number of dots to simulate
%   max_corr_length: the maximum radius out to which the correlation is 
%       calculated, given in nanometers. 
%   dots_per_cell_mean: mean number of dots to simulate per each simulated
%       cell. Assumes a Poission distribuiton of the number of dots per
%       cell.
%   dot_radius: radius of each dot region, given in nanometers.
%   dot_correlation_value: factor for how more/less frequently events are
%       found in the dot regions
%   label_density_mean: mean density of true events in the main body of the 
%      cell, given in events/micrometer^2
%   label_density_stdev: standard deviation of the density of true events 
%       in the main body of the cell, given in events/micrometer^2
%   label_SN_ratio: signal to noise ratio of labels, given as # true
%       events / # spurious events. Can also be the string 'no_noise'.
%   event_overcounting: average number of times each event is overcounted. 
%   dot_center_precision: estimate of the standard deviation of the dot's 
%       true center coordinate. 
%   event_precision: estimate of the localization precision for each 
%       detected event.     
%   STORM_pixel_resolution: size of pixels in the STORM images that are 
%       created, given in nanometers.
%   STORM_method: string indicating the method used to create the STORM
%       images. Options are: 'pdf' or 'binning'
% Outputs:
%   correlation_stack: 3D stack of 2D correlation from each bootstrap.
%   number_cells_screened: the number of cells nedded in to collect all the
%       required dots

% Rename parameter structure for convenience
params = parameter_struct;

% Determine coordinate bounds
min_x_bound = params.bounds(1);
min_y_bound = params.bounds(2);
max_x_bound = params.bounds(3);
max_y_bound = params.bounds(4);
x_length = max_x_bound - min_x_bound;
y_length = max_y_bound - min_y_bound;

% Determine cell center
if strcmp(params.cell_center, 'centered')
    cell_center_x = (x_length)/2;
    cell_center_y = (y_length)/2;
else
    cell_center_x = params.cell_center(1);
    cell_center_y = params.cell_center(2);
end

% Determine cell radius
if isnumeric(params.cell_radius)
    cell_radius = params.cell_radius;
else
    tokens = regexp(params.cell_radius, '(\d+)([%])', 'tokens'); % Look for a percentage value
    if ~isempty(tokens)
        max_radius = min([max_x_bound - cell_center_x, cell_center_x - min_x_bound,...
            max_y_bound - cell_center_y, cell_center_y - min_y_bound]);
        cell_radius = (str2num(tokens{1}{1})/100) * max_radius;
    end
end

% Calculate required constant values
cell_area = pi * (cell_radius / 1000)^2; % in micrometers^2
logn_mu = log((label_density_mean^2) / sqrt(label_density_stdev^2 + label_density_mean^2));
logn_sigma = sqrt(log(label_density_stdev^2 / (label_density_mean^2) + 1));
correlation_max_radius_px = ceil(max_corr_length / STORM_pixel_resolution); 
correlation_width = 2 * correlation_max_radius_px + 1;

% Edit param values for all cells
params.STORM_pixel_size = STORM_pixel_resolution;
params.ch2_crosscor_params = [0, event_precision];

% Initialize counter variables
number_cells_screened = 0;
number_cells_used = 0;
number_dots_complete = 0;

% Create and correlate cells until we have collected enough cells
while number_dots_complete < number_dots
    
    % Choose a number of dots and label density
    num_dots_cell = poissrnd(dots_per_cell_mean);
    if num_dots_cell == 0
        num_cells_screened = num_cells_screened + 1;
        continue % Try again if there are no dots in cell
    end
    label_density_cell = lognrnd(logn_mu, logn_sigma);
    
    % Determine how many dots are still needed
    num_needed = number_dots - number_dots_complete;
    
    % Count cell
    number_cells_screened = number_cells_screened + 1;
    number_cells_used = number_cells_used + 1;
    
    % Calculate required values per cell
    num_base_events = round(cell_area * label_density_cell);   
    num_overcount_events = round(num_base_events * event_overcounting);
    num_true_events = num_base_events + num_overcount_events;
    if isnumeric(label_SN_ratio)
        num_noise_events = num_true_events / label_SN_ratio;
    elseif strcmp(label_SN_ratio, 'no_noise');
        num_noise_events = 0;
    end
    
    % Edit param values for each individual cell
    params.number_events_ch1 = num_base_events;
    params.number_events_ch2 = num_true_events;
    params.number_background_events_ch2 = round(num_noise_events);
    
    % Create a 2D pdf for finding the cells 
    [pdf_map, dot_center_coords, cell_mask] = dots_pdf_map(num_dots_cell, dot_radius, dot_correlation_value, params);
    
    % Get event data for the image
    [~, event_data, ~, STORM_vars ] = create_test_data_dv(params);

    % Create the event STORM image
    data = struct();
    data.x = event_data(:, 1);
    data.y = event_data(:, 2);
    STORM_image = create_STORM_image(data, params.STORM_pixel_size, event_precision, [x_length, y_length], false, false, true);
    
    % Randomly remove extra dot coordinates on last cell
    if num_dots_cell > num_needed
        needed_indices = randperm(num_dots_cell, num_needed);
        dot_center_coords = dot_center_coords(needed_indices, :);
    end
    
    % Save images and coordinates
    STORM_image_cells{number_cells_used} = STORM_image;
    dot_center_cells{number_cells_used} = dot_center_coords;
    
     % Record the number of dots finished in the current cell
    if num_dots_cell <= num_needed
        number_dots_complete = number_dots_complete + num_dots_cell;
    else
        number_dots_complete = number_dots_complete + num_needed;
    end
end

% Repeat the bootstrapping
correlation_stack = zeros(correlation_width, correlation_width, number_resamplings);
for resample_index = 1:number_resamplings

    % Create a new dataset by resampling the existing one
    resampled_center_coords_cells = resample_coordinates(dot_center_cells);
    
    % Correlate the dots in each cell and take a weighted average
    weighted_corr_stack = zeros(correlation_width, correlation_width, length(resampled_center_coords_cells(:)));
    for cell_index = 1:length(resampled_center_coords_cells(:))
        
        % Correlate the dots in the cell and save the correlation
        [cell_correlation] = dots_correlation(STORM_image_cells{cell_index}, resampled_center_coords_cells{cell_index},...
            dot_center_precision, max_corr_length, params.STORM_pixel_size, cell_mask, [0, 0]);
        weighted_corr_stack(:, :, cell_index) = cell_correlation * size(resampled_center_coords_cells{cell_index}, 1);
    end
    
    % Average the weighted cell stack and add to resampling stack
    correlation_stack(:, :, resample_index) = sum(weighted_corr_stack, 3) / number_dots;
    
    % Report to console
    fprintf('.');
end

% Report to console
fprintf('\n');
end

