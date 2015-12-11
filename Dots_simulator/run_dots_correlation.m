% Script to collect correlations from 1D "dots" 

% Set directory information
figure_path_parts = strsplit(pwd, 'CLM_paper');
figure_path = [figure_path_parts{1}, 'CLM_figures_and_data/'];

% Load a parameter structure from current dirctory with deisred settings
load('dots_params.mat') % loads as 'params'



% Create a pdf map for the image
[testimage_ch1_distribution_map, ch1_centers_list] = dots_pdf_map(25, 100, 5, params);

% 














print([figure_path, 'NPIF_SN', num2str(SN_ratio), '_', method_name, '.png'], '-dpng');