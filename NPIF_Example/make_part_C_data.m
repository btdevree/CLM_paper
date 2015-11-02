% Make part C figure data

% Expects to find the cell arrays 'images_array' and 'params_array' in
% the workspace and the string 'binary_path'.

% Get the surprisal value of a zero image
image_size = size(images_array{1});
zero_image_surprisal = calculate_surprisal(zeros(image_size), ideal_image, 'squared_distance');