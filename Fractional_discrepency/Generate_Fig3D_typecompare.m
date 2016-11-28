% Script to generate the picures used in figure 3 
% Run in Fractional_discrepency folder

% Make simulated images of different types
% Use tophat circle image, microtubule images, dots, and a cosine doughnut

% Get path for output
figure_path_parts = strsplit(pwd, 'CLM_paper');
figure_path = [figure_path_parts{1}, 'CLM_figures_and_data/Fig3_v3/'];

% Define Signal/Noise ratio for spurious localizations
SNratio = 4;

% Define various parameters
doughnut_width = 100; % nanometers
number_dots = 12;
dot_radius = 100; % nanometers
number_lines = 6;
line_width = 26; % nanometers, microtubule width

% Load in parameters structure
load('parameters_Fig3.mat'); % loads 'params' into local namespace

% Prep circle and doughnut images

% Randomize the center and radius
center_offset_x = 200 * (.5-rand); % +- 200 nm offset
center_offset_y = 200 * (.5-rand); % +- 200 nm offset
radius_offset = 200 * (.5-rand); % +- 100 nm offset

% Edit parameters
params.cell_center = [750 + center_offset_x; 750 + center_offset_y];
params.cell_radius = 500 + radius_offset;
pdf_map_size = 3; % nanometers
pixel_size = 21; % 7 x supersampeled map
params.ch1_distribution = 'mapped';
params.ch1_distribution_params = {'pdf_map', pdf_map_size, [0, 0]};
params.bounds = [0, 0, 1512, 1512]; % Set image size to be evenly divisible number of pixels

% Get the supersampeled ideal image and the exact image for the circle and doughnut 
params.STORM_pixel_size = pdf_map_size; 
[ideal_image_circle_ss, exact_image_circle] = calculate_ideal_image(params);
[ideal_image_doughnut_ss, exact_image_doughnut] = calculate_ideal_image(params, doughnut_width);

% Get the supersampeled exact image for lines and dots 
params.STORM_pixel_size = pixel_size; 
exact_image_dots = dots_pdf_map(params, number_dots, dot_radius, SNratio);
exact_image_lines = lines_pdf_map(params, number_lines, line_width, SNratio, 'cubic');

% Create supersampeled ideal images for lines and dots
Gauss_radius_in_pixels = params.STORM_precision / pdf_map_size;
filter_kernel = fspecial('gaussian', ceil(7 * Gauss_radius_in_pixels), Gauss_radius_in_pixels);
ideal_image_dots_ss = imfilter(exact_image_dots, filter_kernel);
ideal_image_lines_ss = imfilter(exact_image_lines, filter_kernel);

% Create ideal images
ideal_image_circle = bin_matrix(ideal_image_circle_ss, 7); % Set for 3nm maps and 21 nm STORM images
ideal_image_doughnut = bin_matrix(ideal_image_doughnut_ss, 7);
ideal_image_dots = bin_matrix(ideal_image_dots_ss, 7);
ideal_image_lines = bin_matrix(ideal_image_lines_ss, 7);

% ---- Figure TCI and ECI curves ----

% Define fraction of events, number of events, and replicates
fraction_vector = [0; .005; .01; .02; .04; .06; .08; .1; .15; .2; .25; .3; .35; .4; .5; .6; .7; .8; .9; 1];
number_events_vector = round(logspace(1, 5, 15));
number_replicates = 10; % Number of new datsets
number_pseudoreplicates = 3; % Number of times to split up each dataset

% Get the IIC and ideal discrepency results
pdf_map = exact_image_circle;
[IIC_results_circle, ~, ideal_discrepency_results_circle] = generate_IIC_curve(params, fraction_vector, number_events_vector,...
    Inf, number_replicates, number_pseudoreplicates, 'sum_of_squares', ideal_image_circle, true, true);
pdf_map = exact_image_doughnut;
[IIC_results_doughnut, ~, ideal_discrepency_results_doughnut] = generate_IIC_curve(params, fraction_vector, number_events_vector,...
    Inf, number_replicates, number_pseudoreplicates, 'sum_of_squares', ideal_image_doughnut, true, true);
pdf_map = exact_image_dots;
[IIC_results_dots, ~, ideal_discrepency_results_dots] = generate_IIC_curve(params, fraction_vector, number_events_vector,...
    Inf, number_replicates, number_pseudoreplicates, 'sum_of_squares', ideal_image_dots, true, true);
pdf_map = exact_image_lines;
[IIC_results_lines, ~, ideal_discrepency_results_lines] = generate_IIC_curve(params, fraction_vector, number_events_vector,...
    Inf, number_replicates, number_pseudoreplicates, 'sum_of_squares', ideal_image_lines, true, true);

% Calculate the ECI and TCI values
[ECI_mean_circle, ECI_stdev_circle, TCI_mean_circle, TCI_stdev_circle] =...
    calc_ECI_and_TCI(IIC_results_circle, ideal_discrepency_results_circle, fraction_vector);
[ECI_mean_doughnut, ECI_stdev_doughnut, TCI_mean_doughnut, TCI_stdev_doughnut] =...
    calc_ECI_and_TCI(IIC_results_doughnut, ideal_discrepency_results_doughnut, fraction_vector);
[ECI_mean_dots, ECI_stdev_dots, TCI_mean_dots, TCI_stdev_dots] =...
    calc_ECI_and_TCI(IIC_results_dots, ideal_discrepency_results_dots, fraction_vector);
[ECI_mean_lines, ECI_stdev_lines, TCI_mean_lines, TCI_stdev_lines] =...
    calc_ECI_and_TCI(IIC_results_lines, ideal_discrepency_results_lines, fraction_vector);

% Save data
filename = [figure_path, 'Fig3D_typecomare_ECI_TCI_data.mat'];
save(filename, 'fraction_vector', 'number_events_vector', 'number_replicates', 'number_pseudoreplicates',...
    'ECI_mean_circle', 'ECI_stdev_circle', 'TCI_mean_circle', 'TCI_stdev_circle',...
    'ECI_mean_doughnut', 'ECI_stdev_doughnut', 'TCI_mean_doughnut', 'TCI_stdev_doughnut',...
    'ECI_mean_dots', 'ECI_stdev_dots', 'TCI_mean_dots', 'TCI_stdev_dots',...
    'ECI_mean_lines', 'ECI_stdev_lines', 'TCI_mean_lines', 'TCI_stdev_lines',...
    'ideal_image_circle', 'ideal_image_doughnut', 'ideal_image_dots', 'ideal_image_lines')    

% Combine for graphing
ECI_mean = {ECI_mean_circle; ECI_mean_doughnut; ECI_mean_dots; ECI_mean_lines};
ECI_stdev = {ECI_stdev_circle; ECI_stdev_doughnut; ECI_stdev_dots; ECI_stdev_lines};
TCI_mean = {TCI_mean_circle; TCI_mean_doughnut; TCI_mean_dots; TCI_mean_lines};
TCI_stdev = {TCI_stdev_circle; TCI_stdev_doughnut; TCI_stdev_dots; TCI_stdev_lines};

% Make figure
filename = [figure_path, 'Fig3D_typecompare_TCI_vs_ECI_graph.png'];
make_ECI_vs_TCI_plot(filename, ECI_mean, TCI_mean, {'Circle'; 'Ring'; 'Dots'; 'Lines'}, {'Circle'; 'Ring'; 'Dots'; 'Lines'}, ECI_stdev, TCI_stdev);

% Make images
filename = [figure_path, 'Fig3D_typecompare_circle.png'];
make_nonparameteric_plot(filename, ideal_image_circle, pixel_size, true)
filename = [figure_path, 'Fig3D_typecompare_ring.png'];
make_nonparameteric_plot(filename, ideal_image_doughnut, pixel_size, true)
filename = [figure_path, 'Fig3D_typecompare_dots.png'];
make_nonparameteric_plot(filename, ideal_image_dots, pixel_size, true)
filename = [figure_path, 'Fig3D_typecompare_lines.png'];
make_nonparameteric_plot(filename, ideal_image_lines, pixel_size, true)
