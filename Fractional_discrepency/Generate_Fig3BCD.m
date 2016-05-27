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

% ---- Fig3B IIC curves ----

% Define fraction of events, number of events, and replicates
fraction_vector = [0; .005; .01; .02; .04; .06; .08; .1; .15; .2; .25; .3; .35; .4; .5; .6; .7; .8; .9; 1];
number_events_vector = [1e2; 3e2; 1e3; 3e3; 1e4; 3e4];
number_replicates = 30; % Number of new datsets
number_pseudoreplicates = 3; % Number of times to split up each dataset

% Get the IIC results
[IIC_results] = generate_IIC_curve(params, number_events_vector, SNratio, number_replicates, number_pseudoreplicates, 'sum_of_squares', true, true);

% Calculate mean and standard deviation
IIC_mean = mean(IIC_results, 3);
IIC_stdev = std(IIC_results, 0, 3);

% Save data
filename = [figure_path, 'Fig3B_SSQ_IIC_data.mat'];
save(filename, 'fraction_vector', 'number_events_vector', 'number_replicates', 'number_pseudoreplicates', 'IIC_mean', 'IIC_stdev');

% Make figure
filename = [figure_path, 'Fig3B_SSQ_IIC_graph.png'];
make_IIC_plot(filename, fraction_vector, number_events_vector, IIC_mean, IIC_stdev);
    
% ---- Figure 3C TCI vs ECI curves ----

% Define fraction of events, number of events, and replicates
fraction_vector = [0; .005; .01; .02; .04; .06; .08; .1; .15; .2; .25; .3; .35; .4; .5; .6; .7; .8; .9; 1];
number_events_vector = logspace(2, 5, 12);
number_replicates = 30; % Number of new datsets
number_pseudoreplicates = 3; % Number of times to split up each dataset

% Get the IIC and ideal discrepency results
[IIC_results, ~, ideal_discrepency_results] = generate_IIC_curve(params, number_events_vector, SNratio, number_replicates, number_pseudoreplicates, 'sum_of_squares', ideal_image, true, true);

% Calculate the AOC and ECI of the IIC curves
dFrac = fraction_vector(2:end) - event_fractions(1:end-1);
midpoint_II = (IIC_results(2:end, :, :) + IIC_results(1:end-1, :, :))/2;
AOC = sum(repmat(dFrac, 1, size(midpoint_II, 2), size(midpoint_II, 3)) .* midpoint_II, 1);
ECI_data = (2 * AOC - 1);
ECI_mean = squeeze(mean(ECI_data, 3))';
ECI_stdev = squeeze(std(ECI_data, 0, 3))';

% Calculate the corrosponding TCI 
TCI_data = 1 - (ideal_discrepency_results(end, :, :) ./ ideal_discrepency_results(1, :, :));
TCI_mean = squeeze(mean(TCI_data, 3))';
TCI_stdev = squeeze(std(TCI_data, 0, 3))';


