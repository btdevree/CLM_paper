% Make figure part B, example images with S/N and completeness differences

% Define variable quantities for the S/N and event_number
SN_ratios = [1; 4; 10];
true_event_numbers = [1e4; 1e5; 1e6];

% Generate/define seeds
seeds = round(rand(3,3).*1e6); % Can specify to reproduce the same images

% Load in the image creation parameters as variable params
load 'params_part_B'

% Initialize cell array to hold the images and titles
images_array = cell(3, 3);
subtitle_array = cell(3, 3);

% Create image and subtitle for each combo of SN_ratio and event_number
% Loop through the two parameters
for SN_index = 1:size(SN_ratios, 1)
    SN_ratio = SN_ratios(SN_index);
    for event_num_index = 1:size(true_event_numbers, 1)
        event_number = true_event_numbers(event_num_index);
        
        % Update the parameter structure to create the correct image
        params.number_events_ch1 = event_number;
        params.number_background_events_ch1 = ciel(event_number/SN_ratio);
        
        % Get a new seed
        seed = seeds(SN_index, event_num_index);
        
        % Create data
        [data_ch1, data_ch2, ~, STORMvars] = create_test_data_dv(params, seed);
        
        % Create image
        [~, ~, RGB_image] = create_test_STORM_image_dv(params, data_ch1, data_ch2, STORMvars);
        
        % Take only the green channel
        images_array{SN_index, event_num_index} = RGB_image(:, :, 2);
        
        % Create subtitle
        text = ['S/N = ', num2str(SN_ratio), ', 