function make_parameteric_circle_plot(filename, image, center, radius, pixel_size)
%MAKE_PARAMETRIC_CIRCLE_PLOT Creates a .png image of the parametric fit to 
%   a simulated circular cell image 

% Assumes a Cartesian coordinate system, with the origin at the lower-left
%   corner of the lower-left pixel.
% Inputs:
%   filename - string, full filename and path of image to be created
%   image - image data matrix
%   center - center of circle, given in pixels.
%   radius - radius of circle, given in pixels.
%   pixel_size - size of pixel, in nanometers.
% Output:
%   Image saved with given filename, no return arguments.

% Get image size
[image_height, image_width] = size(image);

% Create new figure
hfig = figure('Units', 'pixels', 'Position', [100, 100, 700, 700], 'PaperPositionMode', 'auto', 'InvertHardcopy', 'off');

% Create clear axes for plotting shapes
haxes = axes('Units', 'pixels', 'Position', [0, 0, 700, 700]);
set(haxes, 'Xlim', [0, image_width], 'Ylim', [0, image_height], 'Color', 'none');
set(haxes, 'XTickLabel', [], 'YTickLabel', [], 'XTick', [], 'YTick', []);

% Create put image on backround axes
hbackaxes = axes('Units', 'pixels', 'Position', get(haxes, 'Position'));
imagesc((image), 'Parent', hbackaxes);
uistack(hbackaxes, 'bottom');% Move the background axes to the bottom
set(hbackaxes, 'XTickLabel', [], 'YTickLabel', [], 'XTick', [], 'YTick', []);
colormap(hbackaxes, gray);

% Plot a circle
rect_diameter = radius * 2;
rect_x = center(1) - radius;
rect_y = center(2) - radius;
rectangle('Position', [rect_x, rect_y, rect_diameter, rect_diameter], 'Curvature', [1,1],...
    'EdgeColor', [1, 0, 0], 'LineStyle', '--', 'LineWidth', 4, 'Parent', haxes);

% Plot a line for radius
line([center(1), center(1) + radius], [center(2), center(2)], 'Color', [1, 0, 0],...
    'LineStyle', '--', 'LineWidth', 4, 'Parent', haxes);

% Plot a filled circle for center
center_dot_radius = radius/20;
rect_diameter = center_dot_radius * 2;
rect_x = center(1) - center_dot_radius;
rect_y = center(2) - center_dot_radius;
rectangle('Position', [rect_x, rect_y, rect_diameter, rect_diameter], 'Curvature', [1,1],...
    'FaceColor', [1, 0, 0], 'EdgeColor', 'none', 'Parent', haxes);

% Print a scale bar
scalebar_position = [30, 30, 100, 20] ./ pixel_size;
rectangle('Position', scalebar_position, 'EdgeColor', 'none', 'FaceColor', [1, 1, 1], 'Parent', haxes);

% Save image
print(hfig, filename, '-dpng');

% Delete figure
delete(hfig);
end

