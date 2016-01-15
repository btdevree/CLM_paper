% Set directory information
figure_path_parts = strsplit(pwd, 'CLM_paper');
figure_path = [figure_path_parts{1}, 'CLM_figures_and_data/'];
data_name_no = 'stdev_vs_num_dots_no_over.mat';
data_name_low = 'stdev_vs_num_dots_low_over.mat';
data_name_high = 'stdev_vs_num_dots_high_over.mat';
figure_prefix = 'stdev_vs_num_dots.mat';

% Load data from no overcounting
load([figure_path, data_name_no])

% 

