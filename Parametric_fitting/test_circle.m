% Script to test circle fitting

% Load parameters
load('parameters_Fig1.mat'); % loads in 'params'
seed = randi(1000);

% Randomize center and radius by a pixel
pixel = params.STORM_pixel_size;
circle_center = [12800 - pixel/2 + pixel*rand; 12800 - pixel/2 + pixel*rand];
circle_center_pixels = circle_center./pixel
params.cell_center = circle_center;
circle_radius = 10240 - pixel/2 + pixel*rand;
circle_radius_pixel = circle_radius/pixel
params.cell_radius = circle_radius;

% Edit parameters
params.number_events_ch1 =1e5;
params.STORM_pixel_resolution = 20;
params.STORM_precision = 25;
params.number_background_events_ch1 = params.number_events_ch1/4;

% Create some data
[data_ch1, data_ch2, passed_vars] = create_test_data_dv(params, seed);

% Create image 
[image] = create_test_STORM_images_dv(params, data_ch1, data_ch2, passed_vars, false, true, true);

% Run parametric fitting
[ center, center_error, radius, radius_error ] = parametric_image_central_circle_finder(image)

% Ideal image
[ ideal_image ] = calculate_ideal_image( params );

% Run parametric fitting
[ center, center_error, radius, radius_error ] = parametric_image_central_circle_finder(ideal_image)