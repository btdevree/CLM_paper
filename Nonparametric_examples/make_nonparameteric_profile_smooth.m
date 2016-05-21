function make_nonparameteric_profile_smooth(filename, image, pixel_size)
%MAKE_NONPARAMETRIC_PROFILE_SMOOTH Creates a .png image of the 1D image 
%   profile of a nonparametric circular region image

% Assumes a Cartesian coordinate system, with the origin at the lower-left
%   corner of the lower-left pixel.
% Inputs:
%   filename - string, full filename and path of image to be created
%   image - image data matrix
%   pixel_size - size of pixel, in nanometers.
% Output:
%   Image saved with given filename, no return arguments.

% Get image size
[image_height, image_width] = size(image);

% Normalize image
image = (image - min(image(:))) / (max(image(:)) - min(image(:)));

% Generate x coordinates
dist_vector = [floor(image_width * 0.1):ceil(image_width * 0.9)] * pixel_size;
intensity_vector = image(floor(image_height * 0.5), floor(image_width * 0.1):ceil(image_width * 0.9));

% Create new figure
hfig = figure('Units', 'pixels', 'Position', [100, 100, 600, 500], 'PaperPositionMode', 'auto', 'InvertHardcopy', 'off', 'Color', [1, 1, 1]);

% Create clear axes for plotting shapes
haxes = axes('Units', 'pixels', 'Position', [70, 50, 500, 400]);
plot(haxes, dist_vector, intensity_vector, 'LineWidth', 4, 'Color', [0, 0, 0]);
set(haxes, 'Xlim', [0, image_width * 0.8 * pixel_size], 'Ylim', [0, 1], 'Color', 'none');
title('1D Intensity Profile', 'FontSize', 16);
xlabel('Distance (nanometers)','FontSize', 12);
ylabel('Normalized Image Intensity','FontSize', 12);

% Save image
print(hfig, filename, '-dpng');

% Delete figure
delete(hfig);
end

