function [correlation_stack, number_cells] = collect_dot_correlations(parameter_struct, number_dots, max_corr_length,...
    dots_per_cell_mean, dot_radius, dot_correlation_value, label_density_mean, label_density_stdev, label_SN_ratio, event_overcounting,...
    dot_center_precision, event_precision, STORM_pixel_resolution, STORM_method)
%COLLECT_DOT_CORRELATIONS Returns a stack of 2D correlations from simulated
%   cells with "dot" type well-defined centers of interest
%
%   This function handles the logic needed to run a dot simulation with the
%   given parameters. 
%   
% Inputs:
%   param_struct: structure that holds the parameters needed for simulating
%       a cell with the create_test_data/create_test_STORM_image
%       functions. Parameters in this structure may be superceded due to 
%       requirements of the other inputs.
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
%       events / # spurious events.
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
%   correlation_stack: 3D stack of the 2D correlations for each dot region.
%   number_cells: the number of cells required to collect the requested
%       number of correlations.

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
max_corr_px = ceil(max_corr_length / STORM_pixel_resolution);
cell_area = pi * cell_radius^2;

% Initialize the results matrix
correlation_stack = zeros(2 * max_corr_px + 1, 2 * max_corr_px + 1, number_dots);

% Edit param values for all cells
params.STORM_pixel_size = STORM_pixel_resolution;
params.ch1_distribution_params = {'pdf_map', 5, [0,0]}; % Constant 5nm pdf map resolution
params.ch2_distribution = 'once_then_random';
params.ch2_distribution_params = 1;
params.ch2_crosscor = 'Gaussian';
params.ch2_crosscor_params = [0, event_precision];
params.number_background_events_ch1 = 0;

% Initialize counter variables
num_cells_screened = 0;
number_dots_complete = 0;

% Create and correlate cells until we have collected enough correlations
while number_dots_complete < number_dots
    
    % Choose a number of dots and label density
    num_dots_cell = poissrnd(dots_per_cell_mean);
    if num_dots_cell == 0
        num_cells_screened = num_cells_screened + 1;
        continue % Try again if there are no dots in cell
    end
    label_density_cell = lognrnd(label_density_mean, label_density_stdev);
    
    % Calculate required values per cell
    num_base_events = round(cell_area * label_density_cell);   
    num_overcount_events = num_base_events * event_overcounting;
    num_true_events = num_base_events + num_overcount_events;
    num_noise_events = num_true_events / label_SN_ratio;
    
    % Edit param values for each individual cell
    params.number_events_ch1 = num_base_events;
    params.number_events_ch2 = num_overcount_events;
    params.number_background_events_ch2 = num_noise_events;
    
    % Create a 2D pdf for finding the cells 
    [pdf_map, dot_center_coords, cell_mask] = dots_pdf_map(num_dots_cell, dot_radius, dot_correlation_value, params);
    
    % Get event data for the image
    [~, event_data, ~, STORM_vars ] = create_test_data_dv(params);
    
    
    
end


end

