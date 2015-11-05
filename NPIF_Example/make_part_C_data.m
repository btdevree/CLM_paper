% Make part C figure data

% Expects to find the cell arrays 'images_array' and 'param_array' in
% the workspace and the string 'binary_path'.

% NOTE: Working on RGB images, not the real STORM images. Need to change
% this.

% Get the ideal image
ideal_image = calculate_ideal_image(param_array{1}); % Assume same relevent parameters for the ideal image of all images

% Turn into RGB to match the images_array from Part B - NOT IDEAL, MUST CHANGE!!!
ideal_RGB_image = convert_STORM_image_to_RGB({ideal_image}, {[0, 1, 0]}, 3, 'channel_max');
ideal_image = double(ideal_RGB_image(:, :, 2));

% Get the surprisal value of a zero image
image_size = size(ideal_image);
zero_image_surprisal = calculate_surprisal(zeros(image_size), ideal_image, 'squared_distance');

% Loop through each image in the array and calculate the surprisal
surprisal_values = zeros(size(images_array));
for index = 1:length(images_array(:));
    
    % Calc surprisal
    surprisal_values(index) = calculate_surprisal(images_array{index}, ideal_image, 'squared_distance');    
end

% Take the ratio of the surprisal to the zero image surprisal
surprisal_ratios = surprisal_values ./ zero_image_surprisal;
