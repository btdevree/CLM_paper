% Script to generate the picures used in figure 3 
% Run in Nonarametric_examples folder

% Make simulated circular structures with uniform label distribution.
% Circle center = middle of image +- .1 um in both directions, circle 
% radius = .5 um +- .1 um.  

% Get path for output
figure_path_parts = strsplit(pwd, 'CLM_paper');
figure_path = [figure_path_parts{1}, 'CLM_figures_and_data/Fig3/'];

% Load in parameters structure
load('parameters_Fig3.mat'); % loads 'params' into local namespace

% Define Signal/Noise ratio for spurious localizations
SNratio = 4;

% Initialize results structure
SSQ_data = struct();

% Initialize RNG
rng('shuffle')
%rng(138290719);
SSQ_data.starting_RNG = rng;

% ---- Part 1 - Zero image -----

% Randomize the center and radius
center_offset_x = 200 * (.5-rand); % +- 200 nm offset
center_offset_y = 200 * (.5-rand); % +- 200 nm offset
radius_offset = 200 * (.5-rand); % +- 100 nm offset

% Define number of events
total_event_number = 3e3;

% Edit parameters
params.cell_center = [750 + center_offset_x; 750 + center_offset_y];
params.cell_radius = 500 + radius_offset;
pixel_size = 20; % Normal grid
params.STORM_pixel_size = pixel_size; 
params.number_events_ch1 = event_number;
params.number_background_events_ch1 = event_number/SNratio;

% Get a RNG seed
seed = randi(1000);

% Create some data
[data_ch1, data_ch2, passed_vars] = create_test_data_dv(params, seed);

% Create image 
[image_100pct] = create_test_STORM_images_dv(params, data_ch1, data_ch2, passed_vars, false, true, true);
image_0pct = zeros(size(image_100pct));

% Calculate SSQ distance
[SSQ_value, image_scale_factor] = calculate_discrepency(image_0pct, image_100pct, 'sum_of_squares');

% Pass image to graphing function
pixel_size = params.STORM_pixel_size;
filename = [figure_path, 'Fig3A_no_data_image.png'];
make_nonparameteric_plot(filename, image_0pct, pixel_size);

% Pass images and scaling factor to profile graphing function
filename = [figure_path, 'Fig3A_no_data_profile.png'];
make_nonparameteric_SSQ_profile(filename, image_0pct, image_100pct, image_scale_factor, pixel_size);

% Record info (in nanometers)
SSQ_data.no_number_events = event_number;
SSQ_data.no_SSQ = SSQ_value;
SSQ_data.no_scale = image_scale_factor;
SSQ_data.no_center_x = params.cell_center(1);
SSQ_data.no_center_y = params.cell_center(2);
SSQ_data.no_radius = params.cell_radius;

% ---- Part 2 - Low data image ----

% Define split fraction
fraction_datapoints = 0.05;

% Generate a logical index array for splitting the data
number_datapoints = size(data_ch1, 1);
number_datapoints_taken = round(number_datapoints * fraction_datapoints);
bool_selection = false(size(data_ch1));
index_integers = randsample(number_datapoints, number_datapoints_taken);                
bool_selection(index_integers) = true;
bool_selection(index_integers + number_datapoints) = true; % select all the y points for corrosponding x in linear representation

% Split the dataset into larger and smaller datasets
data_ch1_taken = reshape(data_ch1(bool_selection), [number_datapoints_taken, 2]);

% Create image with partial datasets
[image_10pct] = create_test_STORM_images_dv(params, data_ch1_taken, data_ch2, passed_vars, false, true, true);

% Calculate SSQ distance
[SSQ_value, image_scale_factor] = calculate_discrepency(image_10pct, image_100pct, 'sum_of_squares');

% Pass image to graphing function
pixel_size = params.STORM_pixel_size;
filename = [figure_path, 'Fig3A_low_data_image.png'];
make_nonparameteric_plot(filename, image_10pct, pixel_size);

% Pass images and scaling factor to profile graphing function
filename = [figure_path, 'Fig3A_low_data_profile.png'];
make_nonparameteric_SSQ_profile(filename, image_10pct, image_100pct, image_scale_factor, pixel_size);

% Record info (in nanometers)
SSQ_data.low_number_events = number_datapoints_taken;
SSQ_data.low_SSQ = SSQ_value;
SSQ_data.low_scale = image_scale_factor;
SSQ_data.low_center_x = params.cell_center(1);
SSQ_data.low_center_y = params.cell_center(2);
SSQ_data.low_radius = params.cell_radius;

% ---- Part 3 - Medium data image ----

% Define split fraction
fraction_datapoints = 0.3;

% Generate a logical index array for splitting the data
number_datapoints = size(data_ch1, 1);
number_datapoints_taken = round(number_datapoints * fraction_datapoints);
bool_selection = false(size(data_ch1));
index_integers = randsample(number_datapoints, number_datapoints_taken);                
bool_selection(index_integers) = true;
bool_selection(index_integers + number_datapoints) = true; % select all the y points for corrosponding x in linear representation

% Split the dataset into larger and smaller datasets
data_ch1_taken = reshape(data_ch1(bool_selection), [number_datapoints_taken, 2]);

% Create image with partial datasets
[image_50pct] = create_test_STORM_images_dv(params, data_ch1_taken, data_ch2, passed_vars, false, true, true);

% Calculate SSQ distance
[SSQ_value, image_scale_factor] = calculate_discrepency(image_50pct, image_100pct, 'sum_of_squares');

% Pass image to graphing function
pixel_size = params.STORM_pixel_size;
filename = [figure_path, 'Fig3A_med_data_image.png'];
make_nonparameteric_plot(filename, image_50pct, pixel_size);

% Pass images and scaling factor to profile graphing function
filename = [figure_path, 'Fig3A_med_data_profile.png'];
make_nonparameteric_SSQ_profile(filename, image_50pct, image_100pct, image_scale_factor, pixel_size);

% Record info (in nanometers)
SSQ_data.med_number_events = number_datapoints_taken;
SSQ_data.med_SSQ = SSQ_value;
SSQ_data.med_scale = image_scale_factor;
SSQ_data.med_center_x = params.cell_center(1);
SSQ_data.med_center_y = params.cell_center(2);
SSQ_data.med_radius = params.cell_radius;

% ---- Part 4 - High data image ----

% Define split fraction
fraction_datapoints = 0.9;

% Generate a logical index array for splitting the data
number_datapoints = size(data_ch1, 1);
number_datapoints_taken = round(number_datapoints * fraction_datapoints);
bool_selection = false(size(data_ch1));
index_integers = randsample(number_datapoints, number_datapoints_taken);                
bool_selection(index_integers) = true;
bool_selection(index_integers + number_datapoints) = true; % select all the y points for corrosponding x in linear representation

% Split the dataset into larger and smaller datasets
data_ch1_taken = reshape(data_ch1(bool_selection), [number_datapoints_taken, 2]);

% Create image with partial datasets
[image_90pct] = create_test_STORM_images_dv(params, data_ch1_taken, data_ch2, passed_vars, false, true, true);

% Calculate SSQ distance
[SSQ_value, image_scale_factor] = calculate_discrepency(image_90pct, image_100pct, 'sum_of_squares');

% Pass image to graphing function
pixel_size = params.STORM_pixel_size;
filename = [figure_path, 'Fig3A_high_data_image.png'];
make_nonparameteric_plot(filename, image_90pct, pixel_size);

% Pass images and scaling factor to profile graphing function
filename = [figure_path, 'Fig3A_high_data_profile.png'];
make_nonparameteric_SSQ_profile(filename, image_90pct, image_100pct, image_scale_factor, pixel_size);

% Record info (in nanometers)
SSQ_data.med_number_events = number_datapoints_taken;
SSQ_data.med_SSQ = SSQ_value;
SSQ_data.med_scale = image_scale_factor;
SSQ_data.med_center_x = params.cell_center(1);
SSQ_data.med_center_y = params.cell_center(2);
SSQ_data.med_radius = params.cell_radius;

% ---- Finish script ----
% Save data
save([figure_path, 'Fig3_SSQ_data.mat'], 'SSQ_data');

% Save data as .csv file for easy reading
temp_table = struct2table(SSQ_data);
writetable(temp_table, [figure_path, 'Fig3_SSQ_data.csv'])