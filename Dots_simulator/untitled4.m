% Script to collect correlations from 1D "dots" 

% Set directory information
figure_path_parts = strsplit(pwd, 'CLM_paper');
figure_path = [figure_path_parts{1}, 'CLM_figures_and_data/'];

% Load a parameter structure from current dirctory with deisred settings
load('dots_params.mat')

% 














print([figure_path, 'NPIF_SN', num2str(SN_ratio), '_', method_name, '.png'], '-dpng');