function make_nonparameteric_SSQ_profile(filename, STORM_image, ideal_image, image_scale_factor, pixel_size)
%MAKE_NONPARAMETRIC_SSQ_profile Creates a .png image of the 1D image 
%   profile of a nonparametric circular region image
%
% Assumes a Cartesian coordinate system, with the origin at the lower-left
%   corner of the lower-left pixel.
% Inputs:
%   filename - string, full filename and path of image to be created.
%   STORM_image/ideal_image - image data matrix.
%   image_scale_factor - scaling factor for the STORM factor to minimize
%       the SSQ discrepency distance.
%   pixel_size - size of pixel, in nanometers.
% Output:
%   Image saved with given filename, no return arguments.

% Get image size
[image_height, image_width] = size(STORM_image);

% Normalize images
STORM_min = min(STORM_image(:));
STORM_max = max(STORM_image(:));
STORM_image = STORM_image - STORM_min;
if STORM_max ~= 0 % Don't divide by zero!
    STORM_image = STORM_image / STORM_max;
end
STORM_image = STORM_image * image_scale_factor; 
ideal_min = min(ideal_image(:));
ideal_max = max(ideal_image(:));
ideal_image = ideal_image - ideal_min;
if ideal_max ~= 0 % Don't divide by zero!
    ideal_image = ideal_image / ideal_max;
end

% Generate x coordinates
dist_vector = [floor(image_width * 0.1):ceil(image_width * 0.9)] * pixel_size;
intensity_vector = STORM_image(floor(image_height * 0.5), floor(image_width * 0.1):ceil(image_width * 0.9));
ideal_dist_tempmatrix = repmat(dist_vector, 2, 1) + repmat([-0.5;0.5] * pixel_size, 1, length(dist_vector)); % upsample the ideal image in order to make a stepped line
ideal_dist_vector = ideal_dist_tempmatrix(:);
ideal_intensity_tempmatrix = repmat(ideal_image(floor(image_height * 0.5), floor(image_width * 0.1):ceil(image_width * 0.9)), 2, 1);
ideal_intensity_vector = ideal_intensity_tempmatrix(:);

% Create new figure
hfig = figure('Units', 'pixels', 'Position', [100, 100, 600, 500], 'PaperPositionMode', 'auto', 'InvertHardcopy', 'off', 'Color', [1, 1, 1]);

% Plot figure
STORM_upperlim = 1.1*max(intensity_vector(:));
if STORM_upperlim == 0 % For all zero images
    STORM_upperlim = 1.1; 
end
haxes = axes('Units', 'pixels', 'Position', [70, 50, 500, 400]);
hbackaxes = axes('Units', 'pixels', 'Position', get(haxes, 'Position'));
bar(hbackaxes, dist_vector, intensity_vector, 'FaceColor', [0, 0, 0], 'EdgeColor', 'none');
hlineplot = plot(haxes, ideal_dist_vector, ideal_intensity_vector, 'LineWidth', 2, 'Color', [0, 0, 0]);
set(haxes, 'Xlim', [min(dist_vector(:)), max(dist_vector(:))], 'Ylim', [0, STORM_upperlim], 'Color', 'none');
set(hbackaxes, 'Xlim', [min(dist_vector(:)), max(dist_vector(:))], 'Ylim', [0, STORM_upperlim], 'Color', 'none');

% Draw SSQ lines
for pixel_index = 1:length(dist_vector);
    line([dist_vector(pixel_index), dist_vector(pixel_index)], [intensity_vector(pixel_index), ideal_intensity_tempmatrix(1, pixel_index)],...
        'Color', [1, 0, 0], 'LineStyle', '-', 'LineWidth', 2, 'Parent', haxes);
end

% Bring line plot in front and bar plot in back
uistack(hlineplot, 'top');
uistack(hbackaxes, 'bottom');

% Titles and labels
title('1D Intensity Profile', 'FontSize', 16);
xlabel('Distance (nanometers)','FontSize', 12);
ylabel('Normalized Image Intensity','FontSize', 12);

% Save image
print(hfig, filename, '-dpng');

% Delete figure
delete(hfig);
end

