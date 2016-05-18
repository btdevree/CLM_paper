function make_parameteric_circle_plot(filename, image, center, radius, pixel_size)
%MAKE_PARAMETRIC_CIRCLE_PLOT Creates a .png image of the parametric fit to 
%   a simulated circular cell image 

function HIV_plots(HIV_filename, data, circle_radius)
% Plot the found HIV particles

% Default values
if nargin < 3; circle_radius = 5; end;

% Load image
image_data = flipud(imrotate(double(imread(HIV_filename)), 90));

% create new figure
figure('Units', 'pixels', 'Position', [50, 50, 900, 900]);

% Create main axes object and handle
haxes = axes('Units', 'pixels', 'Position', [50, 50, 800, 800]);
set(haxes, 'DataAspectRatio', [250,250,1], 'Xlim', [0,512], 'Ylim', [0,512], 'Color', 'none');

% Create the background image axes with a default all black background
hbackaxes = axes('Units', 'pixels', 'Position', get(haxes, 'Position'));
imagesc(flipud(image_data), 'Parent', hbackaxes); % autoflips image data, I think.
uistack(hbackaxes,'bottom');% Move the background axes to the bottom
set(hbackaxes, 'XTickLabel', [], 'YTickLabel', [], 'XTick', [], 'YTick', []);
colormap(hbackaxes, gray);

% Plot a circle around each found spot
rect_diameter = circle_radius .* 2;
rect_x = data.x - circle_radius;
rect_y = data.y - circle_radius;
for rect_ind = 1:size(rect_x, 1) 
    rectangle('Position', [rect_x(rect_ind) - .5, rect_y(rect_ind) - .5, rect_diameter, rect_diameter],...
        'Curvature', [1,1], 'EdgeColor', [1, 1, 0], 'Parent', haxes)
end


end

