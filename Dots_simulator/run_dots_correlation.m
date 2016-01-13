% Script to collect correlations from 1D "dots" 

% Set directory information
figure_path_parts = strsplit(pwd, 'CLM_paper');
figure_path = [figure_path_parts{1}, 'CLM_figures_and_data/'];

% Load a parameter structure from current dirctory with deisred settings
load('dots_params.mat') % loads as 'params'

% Get a results stack
[correlation_stack, number_cells] = collect_dot_correlations(params, 100, 1000, 20, 125, 5, 30, 5, 10, 1, 25, 25, 7, 'pdf');


