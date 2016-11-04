function make_nonparameteric_plot(filename, image, pixel_size, no_line_flag, big_scalebar)
%MAKE_NONPARAMETRIC_PLOT Creates a .png image of a profile line on a 
%   simulated circular cell image 

% Assumes a Cartesian coordinate system, with the origin at the lower-left
%   corner of the lower-left pixel.
% Inputs:
%   filename - string, full filename and path of image to be created
%   image - image data matrix
%   pixel_size - size of pixel, in nanometers.
%   no_line_flag - boolean flag, if true, don't draw the profile line on
%       image. Default = false
%   big_scalebar - boolean flag, if true, the scale bar is 10x bigger than
%       normal
% Output:
%   Image saved with given filename, no return arguments.

% Set defaults
if nargin < 4; no_line_flag = false; end; 
if nargin < 5; big_scalebar = false; end; 

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
if max(image(:)) == 0 % For all zero images
    set(hbackaxes, 'Clim', [0, 1]);
end
uistack(hbackaxes, 'bottom');% Move the background axes to the bottom
set(hbackaxes, 'XTickLabel', [], 'YTickLabel', [], 'XTick', [], 'YTick', []);
colormap(hbackaxes, gray);

% Plot a line across the center
if ~no_line_flag
    line([0.1 * image_width, 0.9 * image_width], [0.5 * image_height, 0.5 * image_height], 'Color', [1, 0, 0],...
        'LineStyle', '-', 'LineWidth', 4, 'Parent', haxes);
end

% Print a scale bar
scalebar_position = [30, 30, 100, 20] ./ pixel_size;
if big_scalebar; scalebar_position = scalebar_position .* 10; end;
rectangle('Position', scalebar_position, 'EdgeColor', 'none', 'FaceColor', [1, 1, 1], 'Parent', haxes);

% Save image
print(hfig, filename, '-dpng');

% Delete figure
delete(hfig);
end

