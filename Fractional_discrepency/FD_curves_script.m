% We assume that we're in the CLM_paper repository, and we want to save the
% big binary figure and data files to CLM_figures_and_data folder, not on 
% the repository but on the same file level as CLM_paper. Assumes we
% can find the images and other data at ..CLM_figures_and_data under the 
% names 'NPIF_part_C_data.mat' and 'NPIF_part_C_images.h5'.

%binary_path_parts = strsplit(pwd, 'CLM_paper');
%binary_path = [binary_path_parts{1}, 'CLM_figures_and_data/'];
binary_path = '/home/btdevree/large_file_temp/'; % Network drive is just too slow and causes process to get killed

% Use the most accurate images
SN_ratio = 10;
optimize_flag = true;

% Choose discrepency method
method = 'sum_of_squares';

% Initialize hdf5 array to hold images directly on disc
image_filepath = [binary_path, 'FD_images_SN', num2str(SN_ratio), '.h5'];
data_filepath = [binary_path, 'FD_data_SN', num2str(SN_ratio), '.mat']; 

% Read in the data file
load(data_filepath);

% Pick the event number datasets to use
event_num_indices = [1:13].';
dataset_replicate_indices = [1:3].';

% Define repeats
number_repeats = 1; 

% Fraction of events to use, must include full and zero fractions
event_fractions = [1; .9; .8; .7; .6; .5; .4; .35; .3; .25; .2; .15; .1; .08; .06; .04; .03; .02; .01; .005; 0];

% Get image size
image_info = h5info(image_filepath, '/images_array');
image_dims = image_info.Dataspace.Size; % Note that Matlab translates this to row-major standard
image_height = image_dims(1);
image_width = image_dims(2);
zero_image = zeros(image_height, image_width);

% Get a copy of the STORMvars structure and channel 2 data created with the given parameter structures
testparams = param_array{1, 1, 1}; % param_array from loading in 'FD_data_SNxx'
testparams.number_events_ch1 = 10;
testparams.number_background_events_ch1 = 10;
[ ~, data_ch2, STORMvars] = create_test_data_dv(testparams, 10);

% Initialize result matrix
ideal_SSQ_results = zeros(length(event_fractions), number_repeats, length(dataset_replicate_indices), length(event_num_indices));
approx_SSQ_results = zeros(length(event_fractions), number_repeats, length(dataset_replicate_indices), length(event_num_indices));

% Loop through each event number
for event_num_index_index = 1:length(event_num_indices)
    event_num_index = event_num_indices(event_num_index_index);
    
    % Loop through each image_replicate
    for dataset_replicate_index_index = 1:length(dataset_replicate_indices)
        dataset_replicate_index = dataset_replicate_indices(dataset_replicate_index_index);

        % Load in the image
        image_100 = h5read(image_filepath, '/images_array', [1, 1, 1, event_num_index, dataset_replicate_index], [image_height, image_width, 1, 1, 1]);

        % Load in the dataset and params
        dataset = datasets{1, event_num_index, dataset_replicate_index};
        params = param_array{1, event_num_index, dataset_replicate_index};
        
        % Give feedback
        tic
        fprintf('Working on event number = %.1g, replicate = %u, repeat = 1', true_event_numbers(event_num_index), dataset_replicate_index);
        
        % Repeat each measurement several times
        for repeat_index = 1:number_repeats
            
            % Give feedback
            fprintf('\b%u', repeat_index);
                      
            % Calculate the result with full image
            ideal_SSQ_results(1, repeat_index, dataset_replicate_index, event_num_index) = calculate_discrepency(image_100, ideal_image, method, optimize_flag);
            approx_SSQ_results(1, repeat_index, dataset_replicate_index, event_num_index) = calculate_discrepency(image_100, image_100, method, optimize_flag);

            % Calculate the result with zero image
            ideal_SSQ_results(end, repeat_index, dataset_replicate_index, event_num_index) = calculate_discrepency(zeros(size(image_100)), ideal_image, method, optimize_flag);
            approx_SSQ_results(end, repeat_index, dataset_replicate_index, event_num_index) = calculate_discrepency(zeros(size(image_100)), image_100, method, optimize_flag);

            % Split dataset and calculate discrepencys, skip first and last fraction 
            for fraction_index = 1:size(event_fractions, 1)-2

                % Generate a logical index array for splitting the data
                number_datapoints = size(dataset, 1);
                number_datapoints_taken = round(number_datapoints * event_fractions(fraction_index + 1));
                bool_selection = false(size(dataset));
                index_integers = randsample(number_datapoints, number_datapoints_taken);                
                bool_selection(index_integers) = true;
                bool_selection(index_integers + number_datapoints) = true; % select all the y points for corrosponding x in linear representation

                % Split the dataset into larger and smaller datasets
                dataset_taken = reshape(dataset(bool_selection), [number_datapoints_taken, 2]);

                % Create images for the larger and smaller datasets
                [partial_image] = create_test_STORM_images_dv(params, dataset_taken, data_ch2, STORMvars, false, true, true);

                % Get the discrepency for the larger and smaller dataset images
                ideal_SSQ_results(fraction_index + 1, repeat_index, dataset_replicate_index, event_num_index) = calculate_discrepency(partial_image, ideal_image, method, optimize_flag);
                approx_SSQ_results(fraction_index + 1, repeat_index, dataset_replicate_index, event_num_index) = calculate_discrepency(partial_image, image_100, method, optimize_flag);
            end     
        end

        % Give feedback
        sec_time = toc;
        fprintf('\bDone. Elapsed time = %.3f min\n', sec_time/60);
    end
end

save(['fractional_discrepency_curves_doublecheck_ssq_SN', num2str(SN_ratio), '.mat'], 'ideal_SSQ_results', 'approx_SSQ_results', 'event_fractions',...
    'number_repeats', 'event_num_indices', 'dataset_replicate_indices', 'true_event_numbers')
