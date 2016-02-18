function image = create_binned_STORM_image(data, resolution, dims, sparse_output_flag)
% CREATE_STORM_IMAGE Creates a high-resolution image built from binning
%   event centers. 
%
% Inputs:
% data: Structure with fields 'x' and 'y'. Assumes a cartesian coordinate 
%   system with the orgin at the bottom left corner of the bottom left 
%   pixel. Given in original pixel values.
% resolution: Desired resolution in fractions of original pixel length.
% dims: Desired image size in original pixel length. Given as a matrix
%   [x_size, y_size].
% sparse_output_flag: logical, default = false. If set to true, image
%   will be retuned as a sparse array.
% 
% Ouptuts:
% image: Super-resolution image of the detected points. Each point has
%   a signal of 1 placed in the nearest pixel. Returned as a
%   sparse or full matrix of floating point doubles.

% Set defaults
if nargin < 4; sparse_output_flag = false; end 

% Calc image parameters
total_number_pixels_x = ceil(dims(1)/resolution);
total_number_pixels_y = ceil(dims(2)/resolution);

% Create STORM image
STORM_image = zeros(total_number_pixels_x, total_number_pixels_y);

% Populate STORM image with coordinates
for coord_index = 1:size(data.x, 1);

    % Calculate the index values for the given coordinates
    column_index = ceil(data.x(coord_index, 1) / resolution); % x coord
    row_index = ceil(total_number_pixels_y - (data.y(coord_index, 1) / resolution)); % y coord

    % Add to the STORM image
    STORM_image(row_index, column_index) = STORM_image(row_index, column_index) + 1;
end

if sparse_output_flag
    % Convert to sparse image to save space
    image = sparse(STORM_image);
else
    image = STORM_image;
end
end