% Make figure part B, example images with S/N and completeness differences

% We assume that we're in the CLM_paper repository, and we want to save the
% big binary figure and data files to CLM_figures_and_data folder, not on 
% the repository but on the same file level as CLM_paper.
binary_path_parts = strsplit(pwd, 'CLM_paper');
binary_path = [binary_path_parts{1}, 'CLM_figures_and_data/'];

% Define variable quantities for the S/N and event_number
SN_ratios = [1; 4; 10];
true_event_numbers = [3e3; 1e5; 3e6];

% Generate seeds 
num_ratios = size(SN_ratios, 1);
num_eventnums = size(true_event_numbers, 1);
seeds = round(rand(num_ratios, num_eventnums).*1e6); % Can specify to reproduce the same images

% Initialize cell array to hold the images and titles
images_array = cell(num_ratios, num_eventnums);
subtitle_array = cell(num_ratios, num_eventnums);

% Load in the image creation parameters as variable params
load 'params_part_B'

% Copy the parameter structure into a cell array so that we can make the images in parallel
% Loop through all SN_ratio and event_number options
param_array = cell(num_ratios, num_eventnums);
for SN_index = 1:size(SN_ratios, 1)
    SN_ratio = SN_ratios(SN_index);
    for event_num_index = 1:size(true_event_numbers, 1)
        event_number = true_event_numbers(event_num_index);
        
        % Edit the parameter structure and make a copy in the parameter array
        params.number_events_ch1 = event_number;
        params.number_background_events_ch1 = ceil(event_number/SN_ratio);
        params.SN_ratio = SN_ratio;
        param_array{SN_index, event_num_index} = params;
    end
end

% Create image and subtitle for each parameter set
parfor index = 1:length(param_array(:))
           
        % Get a new seed
        seed = seeds(index); % Uses a linear index
        param_struct = param_array{index};
        
        % Create data
        [data_ch1, data_ch2, ~, STORMvars] = create_test_data_dv(param_struct, seed);
        
        % Create image
        [~, ~, RGB_image] = create_test_STORM_images_dv(param_struct, data_ch1, data_ch2, STORMvars);
        
        % Take only the green channel
        images_array{index} = RGB_image(:, :, 2);
        
        % Create subtitle and put in array
        text = ['S/N = ', num2str(param_struct.SN_ratio), '; ', sprintf('%1$3.0e', param_struct.number_events_ch1), ' events'];
        subtitle_array{index} = text;
end

% Save the datafiles because they can take a long time to generate
save([binary_path, 'NPIF_part_B_data.mat'], 'binary_path', 'images_array', 'subtitle_array', 'param_array');