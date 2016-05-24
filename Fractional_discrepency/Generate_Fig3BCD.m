% Script to generate the picures used in figure 3 
% Run in Fractional_discrepency folder

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

% Randomize the center and radius
center_offset_x = 200 * (.5-rand); % +- 200 nm offset
center_offset_y = 200 * (.5-rand); % +- 200 nm offset
radius_offset = 200 * (.5-rand); % +- 100 nm offset

% Edit parameters
params.cell_center = [750 + center_offset_x; 750 + center_offset_y];
params.cell_radius = 500 + radius_offset;
pixel_size = 20; % Normal grid
params.STORM_pixel_size = pixel_size; 

% ---- IIC curves ----

% Define fraction of events, number of events, and replicates
fraction_vector = [0; .002; .005; .01; .02; .03; .04; .06; .08; .1; .15; .2; .25; .3; .35; .4; .5; .6; .7; .8; .9; 1];
number_events_vector = [1e3; 3e3; 1e4; 3e4; 1e5; 3e5; 1e6];
number_replicates = 2;

% Initialize IIC result matrix
IIC_results = zeros(length(fraction_vector), length(number_events_vector), number_replicates);

% Loop through each set of parameters to make an IIC curve
for num_events_index = 1:length(number_events_vector)
    number_events = number_events_vector(num_events_index);
    for replicate_index = 1:number_replicates
        
    % Edit parameters
    params.number_events_ch1 = number_events;
    params.number_background_events_ch1 = number_events/SNratio;
    
    % Generate new data
    seed = randi(1e6);
    [dataset] = create_test_data_dv(params, seed);
    
    % Get IIC curve
    IIC_results(:, num_events_index, replicate_index) = calculate_IIC(params, dataset, fraction_vector);
    end
end

% Calculate mean and standard deviation
IIC_mean = mean(IIC_results, 3);
IIC_stdev = std(IIC_results, 0, 3);

% Create figure


    



