% Make figure part B figure

% Expects to find the cell arrays 'images_array' and 'subtitle_array' in
% the workspace and the string 'binary_path'.

% Make new figure
hfig = figure;
set(gcf, 'Position', [100, 100, 1000, 1000]);
set(gcf, 'PaperPositionMode', 'auto');

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
    
    % Increase subplot size axes by a percentage
    fraction_increase = 0.25;
    subplot_position = get(gca, 'Position');
    subplot_position(1) = subplot_position(1) - 0.5 * fraction_increase * subplot_position(3);
    subplot_position(3) = subplot_position(3) + fraction_increase * subplot_position(3);
    subplot_position(2) = subplot_position(2) - 0.5 * fraction_increase * subplot_position(4);
    subplot_position(4) = subplot_position(4) + fraction_increase * subplot_position(4);
    set(gca, 'Position', subplot_position);
    
    % Set title
    title(subtitle_array_t{image_index}, 'Fontsize', 14);
end

% Save figure
print(gcf, [binary_path, 'Figure_part_B.png'],'-r150', '-dpng');