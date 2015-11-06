% Make part C figure data

% We assume that we're in the CLM_paper repository, and we want to save the
% big binary figure and data files to CLM_figures_and_data folder, not on 
% the repository but on the same file level as CLM_paper.
binary_path_parts = strsplit(pwd, 'CLM_paper');
binary_path = [binary_path_parts{1}, 'CLM_figures_and_data/'];

% Define variable quantities for the S/N and event_number
SN_ratios = [1; 4; 10];
true_event_numbers = [1e3; 2e3; 5e3; 1e4; 2e4; 5e4; 1e5; 2e5; 5e5; 1e6; 2e6; 5e6; 1e7];
replicates = 3;

% Generate seeds 
num_ratios = size(SN_ratios, 1);
num_eventnums = size(true_event_numbers, 1);
seeds = round(rand(num_ratios, num_eventnums, replicates).*1e6); % Can specify to reproduce the same images

% Initialize cell arrays to hold the images and datasets, do each SN ratio
% seperate in order to keep from using up all the heap memory
images_array = cell(1, num_eventnums, replicates);
ch1_datasets = cell(num_ratios, num_eventnums, replicates);
ch2_datasets = cell(num_ratios, num_eventnums, replicates);
param_array = cell(num_ratios, num_eventnums, replicates);

% Load in the image creation parameters as variable params
load 'params_part_B'

% Copy the parameter structure into a cell array so that we can make the images in parallel
% Loop through all SN_ratio and event_number options
for SN_index = 1:size(SN_ratios, 1)
    SN_ratio = SN_ratios(SN_index);
    for event_num_index = 1:size(true_event_numbers, 1)
        event_number = true_event_numbers(event_num_index);
        
        % Edit the parameter structure and make a copy in the parameter array
        params.number_events_ch1 = event_number;
        params.number_background_events_ch1 = ceil(event_number/SN_ratio);
        params.SN_ratio = SN_ratio;
        for replicate_index = 1:replicates
            param_array{SN_index, event_num_index, replicate_index} = params;
        end
    end
end

% Calculate the images in parallel; run and save each SN ration seperate
for SN_index = 1:size(SN_ratios, 1)
    SN_ratio = SN_ratios(SN_index);
    
    % Take a single SN ratio group of param structs and seeds
    param_array_group = param_array(SN_index, :, :); % Same size as images_array
    seed_group = seeds(SN_index, :, :);
    ch1_datasets_group = ch1_datasets(SN_index, :, :);
    ch2_datasets_group = ch2_datasets(SN_index, :, :);
    
    % Generate images for the group
    parfor index = 1:length(param_array_group(:))

            % Get a new seed
            seed = seeds(index); % Uses a linear index
            param_struct = param_array_group{index};

            % Create data and save it
            [data_ch1, data_ch2, ~, STORMvars] = create_test_data_dv(param_struct, seed);
            ch1_datasets_group{index} = data_ch1;
            ch2_datasets_group{index} = data_ch2;

            % Create image and save it in cell array
            [image, ~, ~] = create_test_STORM_images_dv(param_struct, data_ch1, data_ch2, STORMvars);
            images_array{index} = image;
    end
    
    % Let user know about progress
    disp(['finished S/N ratio = ', num2str(SN_ratio)]);
    
    % Copy the datasets into the main cell arrays
    ch1_datasets(SN_index, :, :) = ch1_datasets_group;
    ch2_datasets(SN_index, :, :) = ch2_datasets_group;
        
    % Save the group of images
    save([binary_path, 'NPIF_part_C_images_SN', num2str(SN_ratio), '.mat'], 'images_array');
end

% Get the ideal image
ideal_image = calculate_ideal_image(param_array{1}); % Assume same relevent parameters for the ideal image of all images

% Save the datafiles because they can take a long time to generate
save([binary_path, 'NPIF_part_C_data.mat'], 'binary_path', 'param_array', 'SN_ratios', 'true_event_numbers', 'ideal_image', 'ch1_datasets', 'ch2_datasets');

