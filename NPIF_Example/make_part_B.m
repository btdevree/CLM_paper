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
        params.number_background_events_ch1 = ceil(event_number/SN_ratio);
        
        % Get a new seed
        seed = seeds(SN_index, event_num_index);
        
        % Create data
        [data_ch1, data_ch2, ~, STORMvars] = create_test_data_dv(params, seed);
        
        % Create image
        [~, ~, RGB_image] = create_test_STORM_images_dv(params, data_ch1, data_ch2, STORMvars);
        
        % Take only the green channel
        images_array{SN_index, event_num_index} = RGB_image(:, :, 2);
        
        % Create subtitle and put in array
        text = ['S/N = ', num2str(SN_ratio), '; ', sprintf('%1$3.0e', event_number), ' events'];
        subtitle_array{SN_index, event_num_index} = text;
    end
end



% Make new figure
hfig = figure;
set(gcf, 'Position', [100, 100, 1000, 1000]);

% Calculate indices
[num_rows, num_columns] = size(images_array);
num_images = num_rows * num_columns;

% Transpose image and title array because MATLAB uses column-major linear
% indices for subplot ONLY
images_array_t = images_array.';
subtitle_array_t = subtitle_array.';

% Make each subplot
for image_index = 1:num_images
    subplot(num_rows, num_columns, image_index);

    % Make image and set options
    imshow(images_array_t{image_index});
    title(subtitle_array_t{image_index});
    
    % Increase subplot size axes by a percentage
    fraction_increase = 0.25;
    subplot_position = get(gca, 'Position');
    subplot_position(1) = subplot_position(1) - 0.5 * fraction_increase * subplot_position(3);
    subplot_position(3) = subplot_position(3) + fraction_increase * subplot_position(3);
    subplot_position(2) = subplot_position(2) - 0.5 * fraction_increase * subplot_position(4);
    subplot_position(4) = subplot_position(4) + fraction_increase * subplot_position(4);
    set(gca, 'Position', subplot_position);
end

% Save figure
print(gcf, 'Figure_part_B.png', '-dpng');