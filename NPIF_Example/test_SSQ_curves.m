% Check out SSQ-100 behavior as we remove more and more data
%
% Want to determine if a slope at 1, 2, 3% removal can be used to make a
% approximated CI index

% Calculates the completeness index data for the part C figure. Assumes we
% can find the images and other data at ..CLM_figures_and_data under the 
% names 'NPIF_part_C_data.mat' and 'NPIF_part_C_images.h5'.

% We assume that we're in the CLM_paper repository, and we want to save the
% big binary figure and data files to CLM_figures_and_data folder, not on 
% the repository but on the same file level as CLM_paper.
%binary_path_parts = strsplit(pwd, 'CLM_paper');
%binary_path = [binary_path_parts{1}, 'CLM_figures_and_data/'];
binary_path = '/home/btdevree/large_file_temp/'; % Network drive is just too slow and causes process to get killed

% Use the most accurate images
SN_ratio = 10;

% Initialize hdf5 array to hold images directly on disc
image_filepath = [binary_path, 'NPIF_part_C_images_SN', num2str(SN_ratio), '.h5'];
data_filepath = [binary_path, 'NPIF_part_C_data_SN', num2str(SN_ratio), '.mat']; 

% Read in the data file
load(data_filepath);

% Pick the event number dataset to use
event_num_index = 7; % 1e5 events
dataset_replicate_index = 1; 

% Fraction of events to calculate ideal-SSQ and SSQ-100 measures at
event_fractions_larger = [.99; .98; .97; .96; .95; .90; .85; .8; .75; .7; .65; .6; .55]; % Special calculation for 1, 0, and .5
event_fractions_smaller = 1-event_fractions_larger;
event_fractions = [1; event_fractions_larger; .5; flipud(event_fractions_smaller); 0];

% Define repeats
number_repeats = 3;

% Get image size
image_info = h5info(image_filepath, '/images_array');
image_dims = image_info.Dataspace.Size; % Note that Matlab translates this to row-major standard
image_height = image_dims(1);
image_width = image_dims(2);

% Load in the image
image_100 = h5read(image_filepath, '/images_array', [1, 1, 1, event_num_index, dataset_replicate_index], [image_height, image_width, 1, 1, 1]);
zero_image = zeros(size(image_100));

% Load in the dataset and params
dataset = datasets{1, event_num_index, dataset_replicate_index};
params = param_array{1, event_num_index, dataset_replicate_index};

 % Get a copy of the STORMvars structure and channel 2 data created with the given parameter structures
testparams = param_array{1, 1, 1};
testparams.number_events_ch1 = 10;
testparams.number_background_events_ch1 = 10;
[ ~, data_ch2, ~, STORMvars] = create_test_data_dv(testparams, 10);

% Initialize result matrix
ideal_SSQ_results = zeros(size(event_fractions, 1), number_repeats);
approx_SSQ_results = zeros(size(event_fractions, 1), number_repeats);

% Repeat each measurement several times
for repeat_index = 1:number_repeats
    
    % Calculate the result with full image
    ideal_SSQ_results(1, repeat_index) = calculate_discrepency(image_100, ideal_image, 'sum_of_squares');
    approx_SSQ_results(1, repeat_index) = calculate_discrepency(image_100, image_100, 'sum_of_squares');
    
    % Calculate the result with zero image
    ideal_SSQ_results(end, repeat_index) = calculate_discrepency(zeros(size(image_100)), ideal_image, 'sum_of_squares');
    approx_SSQ_results(end, repeat_index) = calculate_discrepency(zeros(size(image_100)), image_100, 'sum_of_squares');
    
    % Split dataset and calculate indices with the larger and smaller halves
    fractions = [event_fractions_larger; 0.5]; % The half-half split will be the one from the "smaller" half image
    for split_index = 1:size(fractions, 1)
        
        % Generate a logical index array for splitting the data
        number_datapoints = size(dataset, 1);
        number_datapoints_larger = round(number_datapoints * fractions(split_index));
        number_datapoints_smaller = number_datapoints - number_datapoints_larger;
        bool_selection = false(size(dataset));
        integers_larger = randsample(number_datapoints, number_datapoints_larger);
        for int_index = 1:length(integers_larger);
            integer_value = integers_larger(int_index);
            bool_selection(integer_value, :) = true;
        end

        % Split the dataset into larger and smaller datasets
        dataset_larger = reshape(dataset(bool_selection), [number_datapoints_larger, 2]);
        dataset_smaller = reshape(dataset(~bool_selection), [number_datapoints_smaller, 2]);

        % Create images for the larger and smaller datasets
        [image_larger, ~, ~] = create_test_STORM_images_dv(params, dataset_larger, data_ch2, STORMvars, false, true, true);
        [image_smaller, ~, ~] = create_test_STORM_images_dv(params, dataset_smaller, data_ch2, STORMvars, false, true, true);

        % Get the discrepency for the larger and smaller dataset images
        ideal_SSQ_results(split_index + 1, repeat_index) = calculate_discrepency(image_larger, ideal_image, 'sum_of_squares');
        ideal_SSQ_results(end - split_index, repeat_index) = calculate_discrepency(image_smaller, ideal_image, 'sum_of_squares');
        approx_SSQ_results(split_index + 1, repeat_index) = calculate_discrepency(image_larger, image_100, 'sum_of_squares');
        approx_SSQ_results(end - split_index, repeat_index) = calculate_discrepency(image_smaller, image_100, 'sum_of_squares');
     end   
end

save('test_SSQ_data_v2.mat', 'ideal_SSQ_results', 'approx_SSQ_results', 'event_fractions',...
    'number_repeats', 'event_num_index', 'dataset_replicate_index')
