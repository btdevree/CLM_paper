function [ ch1_STORM_image, ch2_STORM_image, STORM_RGB_image ] = create_test_STORM_images_dv( parameters_struct, ch1_data, ch2_data, STORM_variables_structure )
%CREATE_TEST_STORM_IMAGES_DV Creates a test movie with the given parameters. Assumes
%   a dualview configuration.

%   Creates a test movie with the given parameters.
%   Inputs:
%   parameters_struct: structure of parameters made according to the
%       template created by the test_movie_parameters function. Please see 
%       that function for more documentation.
%   ch1/ch2_data: matrix of event coordiantes calculated with 
%       create_test_data_dv
%   Outputs:
%   ch1/ch2_STORM_image: sparse double 
%   STORM_RGB_image: An RGB image (8 bit per color channel) for display of 
%       the STORM ground-truth image

% ----------- Create STORM images --------------

% Parameters needed for STORM image creating function
STORM_resolution = params.STORM_pixel_size; % data is given in nm, so resolution is just the nm per pixel
STORM_sigma = params.STORM_precision; % data is given in nm, so precision is in nm
STORM_dims = [max_y_bound - min_y_bound, max_x_bound - min_x_bound]; % data is given in nm, so bounds are in nm

% Create nested function to re-use code for each channel
function [STORM_image] = calc_STORM_image(STORM_resolution, STORM_sigma, STORM_dims, event_coords)

    % Convert points to a data structure
    data_struct = struct();
    data_struct.x = event_coords(:, 1);
    data_struct.y = event_coords(:, 2);

    % Call image generating function, use parallel processing if possible
    STORM_image = create_STORM_image(data_struct, STORM_resolution, STORM_sigma, STORM_dims, true); % note that output is spars
end

% Run STORM image function for each channel
ch1_STORM_image = calc_STORM_image(STORM_resolution, STORM_sigma, STORM_dims, ch1_event_coords);
ch2_STORM_image = calc_STORM_image(STORM_resolution, STORM_sigma, STORM_dims, ch2_event_coords);

% Get a two-color storm image
ch1_color = params.ch1_color;
ch2_color = params.ch2_color;

% Create a display RGB image
STORM_RGB_image = convert_STORM_image_to_RGB({ch1_STORM_image, ch2_STORM_image}, {ch1_color, ch2_color}, 3, 'channel_max');

% Save images as mat files
save('testmovie_STORM_images', 'ch1_STORM_image', 'ch2_STORM_image', 'STORM_RGB_image');

end
