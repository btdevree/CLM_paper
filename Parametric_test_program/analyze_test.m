function analyze_test(test_directory, output_directory)
%ANALYZE_TEST Compare responses of a parametric STORM image test to the 
% ground truth and save the summary statistics.
%
% Requires functions from Common, Actin_simulator, Border_simulator,
% Dots_simulator, and Region_simulator folders; make sure these are on the
% MATLAB path. 
%   
% Input:
%   test_directory: string, filepath for the test folder where the test
%       archive and response are kept. Optional, default = current 
%       directory.
%   output_directory: string, filepath for saving the test results.
%       Optional, default = test_directory.
% Output:
%   Summary is saved to the output directory as "test_summary.mat"

% Set defaults
if nargin < 1; test_directory = ''; end;
if nargin < 2; output_directory = test_directory; end;    

% Start up a new parallel pool using default settings, if needed.
gcp

% Expand and load the test, load the answers
untar([test_directory, '/test_archive.tar.gz'], output_directory);
load([test_directory, '/test_archive.mat']); % Loads test_info
load([test_directory, '/response_info.mat']); % Load responses

% Get required constants and variables
number_images = size(test_info.image_info, 1);
params = test_info.params;

% Initalize variables
ts = struct(); % For test_summary, shorten for convenience
ts.TCI_mean = cell(number_images, 1);
ts.ECI_mean = cell(number_images, 1);
ts.TCI_stdev = cell(number_images, 1);
ts.ECI_stdev = cell(number_images, 1);

% ------Calculate the TCI and ECI of each image-------

% Define fraction of events, number of events, and replicates
fraction_vector = [0; .005; .01; .02; .04; .06; .08; .1; .15; .2; .25; .3; .35; .4; .5; .6; .7; .8; .9; 1];
number_pseudoreplicates = 30; % Number of times to split up each dataset

% Run IIC curve script for each image
for image_index = 1:number_images
    
    % Get the ideal image and dataset
    ideal_image = test_info.ideal_images{image_index};
    dataset = test_info.event_data{image_index};

    % Get the IIC and ideal discrepency results
    [IIC_results, ~, ideal_discrepency_results] = calculate_IIC(params, dataset, fraction_vector, number_pseudoreplicates, 'sum_of_squares', ideal_image, true);

    % Calculate the AOC and ECI of the sum of squares IIC curves
    dFrac = fraction_vector(2:end) - fraction_vector(1:end-1);
    midpoint_II = (IIC_results(2:end, :, :) + IIC_results(1:end-1, :, :))/2;
    AOC = sum(repmat(dFrac, 1, size(midpoint_II, 2), size(midpoint_II, 3)) .* midpoint_II, 1);
    ECI_data = (2 * AOC - 1);
    ts.ECI_mean{image_index} = squeeze(mean(ECI_data, 3))';
    ts.ECI_stdev{image_index} = squeeze(std(ECI_data, 0, 3))';

    % Calculate the corrosponding TCI 
    TCI_data = 1 - (ideal_discrepency_results(end, :, :) ./ ideal_discrepency_results(1, :, :));
    ts.TCI_mean{image_index} = squeeze(mean(TCI_data, 3))';
    ts.TCI_stdev{image_index} = squeeze(std(TCI_data, 0, 3))';
end

% ------Get measurements for Regions---------------

% Find region image indices
region_indices = [];
for image_index = 1:number_images
    info = test_info.ideal_images{image_index};
    if strcmp(info.image_type, 'region')
        region_indices = [region_indices; image_index];
    end
end

% Initalize



end
