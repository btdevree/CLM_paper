% Script to generate the picures used in figure 2 
% Run in Nonparametric_examples folder

% Make simulated circular structures with uniform label distribution.
% Circle center = middle of image +- .1 um in both directions, circle 
% radius = .5 um +- .1 um.  

% Get path for output
figure_path_parts = strsplit(pwd, 'CLM_paper');
figure_path = [figure_path_parts{1}, 'CLM_figures_and_data/Fig2/'];

% Load in parameters structure
load('parameters_Fig2.mat'); % loads 'params' into local namespace

% Define Signal/Noise ratio for spurious localizations
SNratio = 4;

% Initialize results structure
SSQ_data = struct();

% Initialize RNG
%rng('shuffle')
rng(138290719);
SSQ_data.starting_RNG = rng;

% ---- Part 1 - Exact image -----
% Just simulate the fitting, the true exact process is not really achevable
% using numerical, not analytical, techniques

% Randomize the center and radius
center_offset_x = 200 * (.5-rand); % +- 200 nm offset
center_offset_y = 200 * (.5-rand); % +- 200 nm offset
radius_offset = 200 * (.5-rand); % +- 100 nm offset

% Edit parameters
params.cell_center = [750 + center_offset_x; 750 + center_offset_y];
params.cell_radius = 500 + radius_offset;
pixel_size = 2; % Very fine grid
params.STORM_pixel_size = pixel_size; 

% Calculate image
[~, image] = calculate_ideal_image(params);

% Pass image and parameter fit info to graphing function
filename = [figure_path, 'Fig2A_exact_image.png'];
make_nonparameteric_plot(filename, image, pixel_size);
filename = [figure_path, 'Fig2A_exact_profile.png'];
make_nonparameteric_profile_smooth(filename, image, pixel_size);

% Record image generation info
SSQ_data.example_center_x = params.cell_center(1);
SSQ_data.example_center_y = params.cell_center(2);
SSQ_data.example_radius = params.cell_radius;

% ---- Part 2 - Ideal image ----
% Use same image as part 1

% Calculate image
[image] = calculate_ideal_image(params);

% Pass image and parameter fit info to graphing function
filename = [figure_path, 'Fig2A_ideal_image.png'];
make_nonparameteric_plot(filename, image, pixel_size);
filename = [figure_path, 'Fig2A_ideal_profile.png'];
make_nonparameteric_profile_smooth(filename, image, pixel_size);

% ---- Part 3 - Pixelated ideal image ----

% Edit parameters
pixel_size = 20; % Normal grid
params.STORM_pixel_size = pixel_size; 

% Calculate image
[image] = calculate_ideal_image(params);

% Pass image and parameter fit info to graphing function
filename = [figure_path, 'Fig2A_pixel_image.png'];
make_nonparameteric_plot(filename, image, pixel_size);
filename = [figure_path, 'Fig2A_pixel_profile.png'];
make_nonparameteric_profile_pixels(filename, image, pixel_size);

% ---- Part 4 - No data image ----

% Define number of events
event_number = 0;

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
[STORM_image] = create_test_STORM_images_dv(params, data_ch1, data_ch2, passed_vars, false, true, true);

% Calculate image
[ideal_image] = calculate_ideal_image(params);

% Calculate SSQ distance
[SSQ_value, image_scale_factor] = calculate_discrepency(STORM_image, ideal_image, 'sum_of_squares');

% Pass image to graphing function
pixel_size = params.STORM_pixel_size;
filename = [figure_path, 'Fig2B_no_data_image.png'];
make_nonparameteric_plot(filename, STORM_image, pixel_size);

% Pass images and scaling factor to profile graphing function
filename = [figure_path, 'Fig2B_no_data_profile.png'];
make_nonparameteric_SSQ_profile(filename, STORM_image, ideal_image, image_scale_factor, pixel_size);

% Record info (in nanometers)
SSQ_data.no_number_events = event_number;
SSQ_data.no_SSQ = SSQ_value;
SSQ_data.no_scale = image_scale_factor;
SSQ_data.no_center_x = params.cell_center(1);
SSQ_data.no_center_y = params.cell_center(2);
SSQ_data.no_radius = params.cell_radius;

% ---- Part 5 - Low data image ----

% Define number of events
event_number = 3e2;

% Edit parameters
params.number_events_ch1 = event_number;
params.number_background_events_ch1 = event_number/SNratio;

% Get a RNG seed
seed = randi(1000);

% Create some data
[data_ch1, data_ch2, passed_vars] = create_test_data_dv(params, seed);

% Create image 
[STORM_image] = create_test_STORM_images_dv(params, data_ch1, data_ch2, passed_vars, false, true, true);

% Calculate image
[ideal_image] = calculate_ideal_image(params);

% Calculate SSQ distance
[SSQ_value, image_scale_factor] = calculate_discrepency(STORM_image, ideal_image, 'sum_of_squares');

% Pass image to graphing function
pixel_size = params.STORM_pixel_size;
filename = [figure_path, 'Fig2B_low_data_image.png'];
make_nonparameteric_plot(filename, STORM_image, pixel_size);

% Pass images and scaling factor to profile graphing function
filename = [figure_path, 'Fig2B_low_data_profile.png'];
make_nonparameteric_SSQ_profile(filename, STORM_image, ideal_image, image_scale_factor, pixel_size);

% Record info (in nanometers)
SSQ_data.low_number_events = event_number;
SSQ_data.low_SSQ = SSQ_value;
SSQ_data.low_scale = image_scale_factor;
SSQ_data.low_center_x = params.cell_center(1);
SSQ_data.low_center_y = params.cell_center(2);
SSQ_data.low_radius = params.cell_radius;

% ---- Part 6 - Medium data image ----

% Define number of events
event_number = 3e3;

% Edit parameters
params.number_events_ch1 = event_number;
params.number_background_events_ch1 = event_number/SNratio;

% Get a RNG seed
seed = randi(1000);

% Create some data
[data_ch1, data_ch2, passed_vars] = create_test_data_dv(params, seed);

% Create image 
[STORM_image] = create_test_STORM_images_dv(params, data_ch1, data_ch2, passed_vars, false, true, true);

% Calculate image
[ideal_image] = calculate_ideal_image(params);

% Calculate SSQ distance
[SSQ_value, image_scale_factor] = calculate_discrepency(STORM_image, ideal_image, 'sum_of_squares');

% Pass image to graphing function
pixel_size = params.STORM_pixel_size;
filename = [figure_path, 'Fig2B_med_data_image.png'];
make_nonparameteric_plot(filename, STORM_image, pixel_size);

% Pass images and scaling factor to profile graphing function
filename = [figure_path, 'Fig2B_med_data_profile.png'];
make_nonparameteric_SSQ_profile(filename, STORM_image, ideal_image, image_scale_factor, pixel_size);

% Record info (in nanometers)
SSQ_data.med_number_events = event_number;
SSQ_data.med_SSQ = SSQ_value;
SSQ_data.med_scale = image_scale_factor;
SSQ_data.med_center_x = params.cell_center(1);
SSQ_data.med_center_y = params.cell_center(2);
SSQ_data.med_radius = params.cell_radius;

% ---- Part 7 - High data image ----

% Define number of events
event_number = 3e4;

% Edit parameters
params.number_events_ch1 = event_number;
params.number_background_events_ch1 = event_number/SNratio;

% Get a RNG seed
seed = randi(1000);

% Create some data
[data_ch1, data_ch2, passed_vars] = create_test_data_dv(params, seed);

% Create image 
[STORM_image] = create_test_STORM_images_dv(params, data_ch1, data_ch2, passed_vars, false, true, true);

% Calculate image
[ideal_image] = calculate_ideal_image(params);

% Calculate SSQ distance
[SSQ_value, image_scale_factor] = calculate_discrepency(STORM_image, ideal_image, 'sum_of_squares');

% Pass image to graphing function
pixel_size = params.STORM_pixel_size;
filename = [figure_path, 'Fig2B_high_data_image.png'];
make_nonparameteric_plot(filename, STORM_image, pixel_size);

% Pass images and scaling factor to profile graphing function
filename = [figure_path, 'Fig2B_high_data_profile.png'];
make_nonparameteric_SSQ_profile(filename, STORM_image, ideal_image, image_scale_factor, pixel_size);

% Record info (in nanometers)
SSQ_data.high_number_events = event_number;
SSQ_data.high_SSQ = SSQ_value;
SSQ_data.high_scale = image_scale_factor;
SSQ_data.high_center_x = params.cell_center(1);
SSQ_data.high_center_y = params.cell_center(2);
SSQ_data.high_radius = params.cell_radius;

% ---- Finish script ----
% Save data
save([figure_path, 'Fig2_SSQ_data.mat'], 'SSQ_data');

% Save data as .csv file for easy reading
temp_table = struct2table(SSQ_data);
writetable(temp_table, [figure_path, 'Fig2_SSQ_data.csv'])