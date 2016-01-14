% Script to collect correlations from 1D "dots" 

% Set directory information
figure_path_parts = strsplit(pwd, 'CLM_paper');
figure_path = [figure_path_parts{1}, 'CLM_figures_and_data/'];

% Load a parameter structure from current dirctory with deisred settings
load('dots_params.mat') % loads as 'params'

% Get a results stack
%parameter_struct, number_dots, max_corr_length, dots_per_cell_mean, dot_radius, dot_correlation_value, 
%   label_density_mean, label_density_stdev, label_SN_ratio, event_overcounting, dot_center_precision, event_precision,
%   STORM_pixel_resolution, STORM_method
[correlation_stack, number_cells] = collect_dot_correlations(params, 50, 1000, 20, 125, 5, 30, 5, 10, 1, 25, 25, 7, 'pdf');


