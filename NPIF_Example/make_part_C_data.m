% Make part C figure data

% Calculates the completeness index data for the part C figure. Assumes we
% can find the images and other data at ..CLM_figures_and_data under the 
% names 'NPIF_part_C_data.mat' and 'NPIF_part_C_images.h5'.

% We assume that we're in the CLM_paper repository, and we want to save the
% big binary figure and data files to CLM_figures_and_data folder, not on 
% the repository but on the same file level as CLM_paper.
binary_path_parts = strsplit(pwd, 'CLM_paper');
binary_path = [binary_path_parts{1}, 'CLM_figures_and_data/'];
% binary_path = '/home/btdevree/large_file_temp/'; % Network drive is just too slow and causes process to get killed

% Read in the data file
load([binary_path, 'NPIF_part_C_data_SN10.mat']);

% Define number of pseudoreplicates for the approximation completeness index, get number of replicates
approx_pseudoreplicates = 5;
replicates = size(param_array, 5);

% Get a copy of the STORMvars structure and channel 2 data created with the given parameter structures
testparams = param_array{1, 1, 1};
testparams.number_events_ch1 = 10;
testparams.number_background_events_ch1 = 10;
[ ~, data_ch2, ~, STORMvars] = create_test_data_dv(testparams, 10);

% Generate seeds 
seeds = round(rand(size(param_array)).*1e6); % Can specify to reproduce the same images

% Initialize a matrix for the completeness index data
ideal_completeness_index = zeros(size(param_array));
approx_completeness_index = zeros([size(param_array), approx_pseudoreplicates]);

% Get image size
imagefile_pathname = [binary_path, 'NPIF_part_C_images_SN10.h5'];
image_info = h5info(imagefile_pathname, '/images_array');
image_dims = image_info.Dataspace.Size; % Note that Matlab translates this to row-major standard
image_height = image_dims(1);
image_width = image_dims(2);

% Loop through each parameter struture and corrosponding image
for param_index = 1:length(param_array(:));
    params = param_array{param_index};
    
    % Channge the random number generator seed
    rng(seeds(param_index));
    
    % Load in the image
    image = h5read(imagefile_pathname, '/images_array', [1, 1, params.output_indices], [image_height, image_width, 1, 1, 1]);
    
    % Calculate the ideal completeness index
    ideal_discrepency = calculate_discrepency(image, ideal_image, 'normalized_variation_of_information', 100, false);
    ideal_completeness_index(param_index) = 1-ideal_discrepency;
    
    % Get the full dataset
    data_100pct = datasets{param_index};
    
    % Calculate the approximated completeness index
    % Repeat multiple times for pseudo-replicates 
    for pseudorep_index = 1:approx_pseudoreplicates
            
        % Generate a logical index array for splitting the data
        number_datapoints_100pct = size(data_100pct, 1);
        number_datapoints_1pct = ceil(number_datapoints_100pct/100);
        number_datapoints_99pct = number_datapoints_100pct - number_datapoints_1pct;
        bool_selection = false(size(data_100pct));
        integers_1pct = randsample(number_datapoints_100pct, number_datapoints_1pct);
        for int_index = 1:length(integers_1pct);
            integer_value = integers_1pct(int_index);
            bool_selection(integer_value, :) = true;
        end
        
        % Split the dataset into 1% and 99% datasets
        data_1pct = reshape(data_100pct(bool_selection), [number_datapoints_1pct, 2]);
        data_99pct = reshape(data_100pct(~bool_selection), [number_datapoints_99pct, 2]);
        
        % Create images for the 1% and 99% datasets
        [image_1pct, ~, ~] = create_test_STORM_images_dv(params, data_1pct, data_ch2, STORMvars, false);
        [image_99pct, ~, ~] = create_test_STORM_images_dv(params, data_99pct, data_ch2, STORMvars, false);
        
        % Get the inherent information increase in the 1% image
        discrepency_1pct = calculate_discrepency(image_1pct, image_99pct, 'normalized_variation_of_information', 100, false);
        
        % Get the information increase by adding the 1% image to the 99% image
        discrepency_100pct = calculate_discrepency(image_1pct, image_99pct, 'normalized_variation_of_information', 100, false);
        
        % Calculate completeness approximation
        approx_completeness_index([params.output_indices, pseudorep_index]) = 1-(discrepency_100pct/discrepency_1pct);
    end
end

save('NPIF_part_C_completeness_data.mat', 'ideal_completeness_index', 'approx_completeness_index', 'param_array', 'SN_ratios', 'true_event_numbers');