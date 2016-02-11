function [ ch1_STORM_image, ch2_STORM_image, STORM_RGB_image ] = create_test_STORM_images_dv( parameters_structure, ch1_data, ch2_data,...
    passed_variables_structure, sparse_output_flag, parallel_image_creation_flag, use_MEX_flag )
%CREATE_TEST_STORM_IMAGES_DV Creates a test movie with the given parameters. Assumes
%   a dualview configuration.

%   Creates a test movie with the given parameters.
%   Inputs:
%   parameters_struct: structure of parameters made according to the
%       template created by the test_movie_parameters function. Please see 
%       that function for more documentation.
%   ch1/ch2_data: matrix of event coordiantes calculated with 
%       create_test_data_dv. Use empty matrix if no ch2 image is requested.
%   passed_variables_structure: Required variables from
%       test_movie_parameters_dv
%   sparse_output_flag: if true, output image is a sparse matrix. Default =
%       true.
%   parallel_image_creation_flag; if true, the image is created using all
%       available cores. Default = true.
%   use_MEX_flage: logical, default = false. If set to true, the core image 
%       creation will be done with the Gauss_STORM_image_MEX command, which
%       needs to be pre-compiled from the Gauss_STORM_image_MEX.cpp source.
%   Outputs:
%   ch1/ch2_STORM_image: sparse double 
%   STORM_RGB_image: An RGB image (8 bit per color channel) for display of 
%       the STORM ground-truth image
%   NOTE: can return only the ch1/ch2 images or just the ch1 image if
%   called like: 
%       [ch1_image, ch2_image] = create_test_STORM_images_dv(...) or 
%       [ch1_image] = create_test_STORM_images_dv(...), respectively.

% ------ Collect and rename variables ------------

% Set defaults
if nargin < 5, sparse_output_flag = true; end;
if nargin < 6, parallel_image_creation_flag = true; end;
if nargin < 7, use_MEX_flag = false; end;

% Rename the data that was generated with create_test_data_dv
ch1_event_coords = ch1_data;
ch2_event_coords = ch2_data;

% Unpack required info from passed_variables_structure
s = passed_variables_structure;
min_x_bound = s.min_x_bound;
min_y_bound = s.min_y_bound;
max_x_bound = s.max_y_bound;
max_y_bound = s.max_y_bound;

% Rename parameters_structure for convenience
params = parameters_structure; 

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

    % Call image generating function, use parallel processing if requested
    STORM_image = create_STORM_image(data_struct, STORM_resolution, STORM_sigma, STORM_dims, sparse_output_flag, parallel_image_creation_flag, use_MEX_flag);
end

% Run STORM image function 
ch1_STORM_image = calc_STORM_image(STORM_resolution, STORM_sigma, STORM_dims, ch1_event_coords);

% Run only if the second image channel is requested
if nargout > 1
    ch2_STORM_image = calc_STORM_image(STORM_resolution, STORM_sigma, STORM_dims, ch2_event_coords);

    % Run only if the RGB image is requested
    if nargout > 2
        
        % Get a two-color storm image
        ch1_color = params.ch1_color;
        ch2_color = params.ch2_color;

        % Create a display RGB image
        STORM_RGB_image = convert_STORM_image_to_RGB({ch1_STORM_image, ch2_STORM_image}, {ch1_color, ch2_color}, 3, 'channel_max');
    end
end

% Save images as mat files
if params.save_STORM_images_flag   
    if nargout == 3
        save('testmovie_STORM_images', 'ch1_STORM_image', 'ch2_STORM_image', 'STORM_RGB_image');
    elseif nargout == 2
        save('testmovie_STORM_images', 'ch1_STORM_image', 'ch2_STORM_image');
    elseif nargout == 1
        save('testmovie_STORM_images', 'ch1_STORM_image');
    end
end
end

