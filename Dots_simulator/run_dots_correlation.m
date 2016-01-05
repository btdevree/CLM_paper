% Script to collect correlations from 1D "dots" 

% Set directory information
figure_path_parts = strsplit(pwd, 'CLM_paper');
figure_path = [figure_path_parts{1}, 'CLM_figures_and_data/'];

% Load a parameter structure from current dirctory with deisred settings
load('dots_params.mat') % loads as 'params'

% Create a pdf map for the image
[pdf_map, center_list, mask] = dots_pdf_map(25, 150, 5, params, true);

% Create a dataset
map_resolution = params.ch1_distribution_params{2};
[ xy_data ] = sample_2D_pdf( 1e5, pdf_map, map_resolution, true);

% Convert xy data and create a STORM image
data = struct();
data.x = xy_data(:, 1);
data.y = xy_data(:, 2);
STORM_dims = params.bounds(3:4);
STORM_image = create_STORM_image(data, 7, 25, STORM_dims, false, true, true);

% Get a stack of correlations for the dots
[correlation_stack] = dots_correlation_individual(STORM_image, center_list, 25, 1000, 7, mask);












print([figure_path, 'NPIF_SN', num2str(SN_ratio), '_', method_name, '.png'], '-dpng');