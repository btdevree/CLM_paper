% Script to generate the picures used in figure 2 
% Run in Parametric_fitting folder

% Make simulated circular structures with uniform label distribution.
% Circle center = middle of image +- .1 um in both directions, circle 
% radius = .5 um +- .1 um.  

% Get path for output
figure_path_parts = strsplit(pwd, 'CLM_paper');
figure_path = [figure_path_parts{1}, 'CLM_figures_and_data/Fig1/'];

% Load in parameters structure
load('parameters_Fig2.mat'); % loads 'params' into local namespace

% Define Signal/Noise ratio for spurious localizations
SNratio = 4;

% Initialize results structure
fit_data = struct();
reference_data = struct();

% Initialize RNG
rng('shuffle')
reference_data.starting_RNG = rng;

% ---- Part 1 - Exact image fitting -----
% Just simulate the fitting, the true exact process is not really achevable
% using numerical, not analytical, techniques

% Randomize the center and radius
center_offset_x = 200 * (.5-rand); % +- 200 nm offset
center_offset_y = 200 * (.5-rand); % +- 200 nm offset
radius_offset = 200 * (.5-rand); % +- 100 nm offset

% Edit parameters
params.cell_center = [750 + center_offset_x; 750 + center_offset_y];
params.cell_radius = 500 + radius_offset;
params.STORM_pixel_size = 2; % Very fine grid

% Calculate image
[~, image] = calculate_ideal_image(params);

% Pass image and parameter fit info to graphing function
filename = [figure_path, 'Fig2A_exact.png'];
make_nonparameteric_plot(filename, image, pixel_size);

% Record info (in nanometers)
fit_data.exact_center_x = params.cell_center(1);
fit_data.exact_center_x_error = 0;
fit_data.exact_center_y = params.cell_center(2);
fit_data.exact_center_y_error = 0;
fit_data.exact_radius = params.cell_radius;
fit_data.exact_radius_error = 0;
reference_data.exact_center_x = params.cell_center(1);
reference_data.exact_center_y = params.cell_center(2);
reference_data.exact_radius = params.cell_radius;

% ---- Part 2 - Ideal image fitting ----
% Use same image as part 1

% Calculate image
[image] = calculate_ideal_image(params);

% Run parametric fitting
[center, center_error, radius, radius_error] = parametric_image_central_circle_finder(image);

% Pass image and parameter fit info to graphing function
filename = [figure_path, 'Fig1A_ideal_fit.png'];
pixel_size = params.STORM_pixel_size;
make_parameteric_circle_plot(filename, image, center, radius, pixel_size);

% Record info (in nanometers)
fit_data.ideal_center_x = pixel_size * center(1);
fit_data.ideal_center_x_error = pixel_size * center_error(1);
fit_data.ideal_center_y = pixel_size * center(2);
fit_data.ideal_center_y_error = pixel_size * center_error(2);
fit_data.ideal_radius = pixel_size * radius;
fit_data.ideal_radius_error = pixel_size * radius_error;
reference_data.ideal_center_x = params.cell_center(1);
reference_data.ideal_center_y = params.cell_center(2);
reference_data.ideal_radius = params.cell_radius;

% ---- Part 3 - Low data image fitting ----
% Define number of events
event_number = 3e2;

% Randomize the center and radius
center_offset_x = 200 * (.5-rand); % +- 200 nm offset
center_offset_y = 200 * (.5-rand); % +- 200 nm offset
radius_offset = 200 * (.5-rand); % +- 100 nm offset

% Edit parameters
params.cell_center = [750 + center_offset_x; 750 + center_offset_y];
params.cell_radius = 500 + radius_offset;
params.number_events_ch1 = event_number;
params.number_background_events_ch1 = event_number/SNratio;

% Get a RNG seed
seed = randi(1000);

% Create some data
[data_ch1, data_ch2, passed_vars] = create_test_data_dv(params, seed);

% Create image 
[image] = create_test_STORM_images_dv(params, data_ch1, data_ch2, passed_vars, false, true, true);

% Run parametric fitting
[center, center_error, radius, radius_error] = parametric_image_central_circle_finder(image);

% Pass image and parameter fit info to graphing function
filename = [figure_path, 'Fig1B_low_fit.png'];
pixel_size = params.STORM_pixel_size;
make_parameteric_circle_plot(filename, image, center, radius, pixel_size);

% Record info (in nanometers)
fit_data.low_number_events = event_number;
fit_data.low_center_x = pixel_size * center(1);
fit_data.low_center_x_error = pixel_size * center_error(1);
fit_data.low_center_y = pixel_size * center(2);
fit_data.low_center_y_error = pixel_size * center_error(2);
fit_data.low_radius = pixel_size * radius;
fit_data.low_radius_error = pixel_size * radius_error;
reference_data.low_center_x = params.cell_center(1);
reference_data.low_center_y = params.cell_center(2);
reference_data.low_radius = params.cell_radius;

% ---- Part 4 - Medium data image fitting ----
% Define number of events
event_number = 3e3;

% Randomize the center and radius
center_offset_x = 200 * (.5-rand); % +- 200 nm offset
center_offset_y = 200 * (.5-rand); % +- 200 nm offset
radius_offset = 200 * (.5-rand); % +- 100 nm offset

% Edit parameters
params.cell_center = [750 + center_offset_x; 750 + center_offset_y];
params.cell_radius = 500 + radius_offset;
params.number_events_ch1 = event_number;
params.number_background_events_ch1 = event_number/SNratio;

% Get a RNG seed
seed = randi(1000);

% Create some data
[data_ch1, data_ch2, passed_vars] = create_test_data_dv(params, seed);

% Create image 
[image] = create_test_STORM_images_dv(params, data_ch1, data_ch2, passed_vars, false, true, true);

% Run parametric fitting
[center, center_error, radius, radius_error] = parametric_image_central_circle_finder(image);

% Pass image and parameter fit info to graphing function
filename = [figure_path, 'Fig1B_med_fit.png'];
pixel_size = params.STORM_pixel_size;
make_parameteric_circle_plot(filename, image, center, radius, pixel_size);

% Record info (in nanometers)
fit_data.med_event_number = event_number;
fit_data.med_center_x = pixel_size * center(1);
fit_data.med_center_x_error = pixel_size * center_error(1);
fit_data.med_center_y = pixel_size * center(2);
fit_data.med_center_y_error = pixel_size * center_error(2);
fit_data.med_radius = pixel_size * radius;
fit_data.med_radius_error = pixel_size * radius_error;
reference_data.med_center_x = params.cell_center(1);
reference_data.med_center_y = params.cell_center(2);
reference_data.med_radius = params.cell_radius;

% ---- Part 5 - High data image fitting ----
% Define number of events
event_number = 3e4;

% Randomize the center and radius
center_offset_x = 200 * (.5-rand); % +- 200 nm offset
center_offset_y = 200 * (.5-rand); % +- 200 nm offset
radius_offset = 200 * (.5-rand); % +- 100 nm offset

% Edit parameters
params.cell_center = [750 + center_offset_x; 750 + center_offset_y];
params.cell_radius = 500 + radius_offset;
params.number_events_ch1 = event_number;
params.number_background_events_ch1 = event_number/SNratio;

% Get a RNG seed
seed = randi(1000);

% Create some data
[data_ch1, data_ch2, passed_vars] = create_test_data_dv(params, seed);

% Create image 
[image] = create_test_STORM_images_dv(params, data_ch1, data_ch2, passed_vars, false, true, true);

% Run parametric fitting
[center, center_error, radius, radius_error] = parametric_image_central_circle_finder(image);

% Pass image and parameter fit info to graphing function
filename = [figure_path, 'Fig1B_high_fit.png'];
pixel_size = params.STORM_pixel_size;
make_parameteric_circle_plot(filename, image, center, radius, pixel_size);

% Record info (in nanometers)
fit_data.high_event_number = event_number;
fit_data.high_center_x = pixel_size * center(1);
fit_data.high_center_x_error = pixel_size * center_error(1);
fit_data.high_center_y = pixel_size * center(2);
fit_data.high_center_y_error = pixel_size * center_error(2);
fit_data.high_radius = pixel_size * radius;
fit_data.high_radius_error = pixel_size * radius_error;
reference_data.high_center_x = params.cell_center(1);
reference_data.high_center_x = params.cell_center(2);
reference_data.high_radius = params.cell_radius;

% ---- Finish script ----
% Save data
save([figure_path, 'Fig1_data.mat'], 'fit_data', 'reference_data');

% Save data as .csv file for easy reading
temp_table = struct2table(fit_data);
writetable(temp_table, [figure_path, 'Fig1_fit_data.csv'])
temp_table = struct2table(reference_data);
writetable(temp_table, [figure_path, 'Fig1_reference_data.csv'])