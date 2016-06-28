function make_test( parameter_structure, output_directory, test_number)
%MAKE_TEST Prepare a parametric STORM image test and save it into the 
% specified directory.
%
% Requires functions from Common, Actin_simulator, Border_simulator,
% Dots_simulator, and Region_simulator folders; make sure these are on the
% MATLAB path. 
%   
% Input:
%   parameter_struct: parameter structure made with
%       test_movie_parameters_dv.
%   output_directory: string, filepath for saving the test folder
%   test_number; integer, the number to assign to the test
% Output:
%   Test is created and saved to the output directory in a folder called
%   "Test_<test_number>". Inside the folder, there is an archived .mat file
%   called "test_archive.tar.bz" containing all the information about the
%   test and a file "test_files.mat" with the information needed by the
%   parametric STORM image testing GUI. 

% Rename parameter structure for convenience
params = parameter_structure;

% ---- Define the test images -----

% Define test characteristics
test_version = '0.1';

% Define important pdf and STORM image parameters just in case the given params don't have this set correctly
% NOTE: If you want to change these, you'll probalby have to change the code below to make sure the ideal image is still correct.
params.ch1_distribution_params{2} = 3; % pdf_map resolution, nanometers, should give a ~73MB map with a 25.6 x 25.6 um field
params.STORM_pixel_size = 21; % works fine for 25 nm STORM precision, ideal image is 1/7th the pdf map (0.024% alaised signal in Fourier domain)

% Regions
region_number_images = 6;
region_contrast_ratios = [0.2, 1, 4]; % background:additional_density
region_event_number_range = [3e2, 1e6];
region_number_vertices = [6];

% Dots
dots_number_images = 15;
dots_contrast_ratios = [1, 4, 10];
dots_sizes = [20, 50, 100, 200, 500]; % nanometers
dots_event_number_range = [3e2, 3e6];
dots_number_dots = [10];

% Actin lines
actin_number_images = 18;
actin_contrast_ratios = [1, 4, 10];
actin_event_number_range = [3e2, 1e7];
actin_line_types = {'line_segment', 'cubic', 'quadratic'};
actin_line_widths = [9, 26]; % actin and microtubules, respectivly
actin_number_lines = [6];

% Border
border_number_images = 15;
border_contrast_ratios = [0.2, 1, 4];
border_event_number_range = [1e3, 1e7];
border_roughness = [.35, .45, .55, .65, .75];
border_displacement_factor = [.4, .35, .3, .25, .2]; % Apply together with the above roughness factor, not for all combinations possible

% ---- Determine the required parameters for each image in the test ----

% Add up number of images
total_number_images = region_number_images + dots_number_images + actin_number_images + border_number_images;

% Shuffle the RNG
rng('shuffle');

% Create folder
mkdir(output_directory, ['Test_', num2str(test_number)]);

% Create test structure and fill with basic information
test_info = struct();
test_info.test_number = test_number;
test_info.test_version = test_version;
test_info.RNG_state = rng();
test_info.params = params;
test_info.image_info = cell(total_number_images, 1);

% Prep region info
image_contrast_ratios = choose_evenly(region_number_images, region_contrast_ratios);
image_event_numbers = choose_from_log_range(region_number_images, region_event_number_range);
image_number_vertices = choose_evenly(region_number_images, region_number_vertices);
info_cells = cell(region_number_images, 1);

% Write region info into image info cells
for image_index = 1:region_number_images
    info_struct = struct();
    info_struct.image_type = 'region';
    info_struct.contrast_ratio = image_contrast_ratios(image_index);
    info_struct.event_number = image_event_numbers(image_index);
    info_struct.number_vertices = image_number_vertices(image_index);
    info_cells{image_index} = info_struct;
end
test_info.image_info(1:region_number_images) = info_cells;
current_image_number = region_number_images;

% Prep dots info
image_contrast_ratios = choose_evenly(dots_number_images, dots_contrast_ratios);
image_event_numbers = chose_from_log_range(dots_number_image, dots_event_number_range);
image_dot_sizes = choose_evenly(dots_number_images, dots_sizes);
image_number_dots = choose_evenly(dots_number_images, dots_number_dots);
info_cells = cell(dots_number_images, 1);

% Write dots info into image info cells
for image_index = 1:dots_number_images
    info_struct = struct();
    info_struct.image_type = 'dots';
    info_struct.contrast_ratio = image_contrast_ratios(image_index);
    info_struct.event_number = image_event_numbers(image_index);
    info_struct.dot_sizes = image_dot_sizes(image_index);
    info_struct.number_dots = image_number_dots(image_index);
    info_cells{image_index} = info_struct;
end
test_info.image_info(current_image_number:current_image_number + dots_number_images) = info_cells;
current_image_number = current_image_number + dots_number_images;

% Prep actin info
image_contrast_ratios = choose_evenly(actin_number_images, actin_contrast_ratios);
image_event_numbers = choose_from_log_range(actin_number_image, actin_event_number_range);
image_line_widths = choose_evenly(actin_number_images, actin_line_widths);
image_line_types = choose_evenly(actin_number_images, actin_line_types);
image_number_lines = choose_evenly(actin_number_images, actin_number_lines);
info_cells = cell(actin_number_images, 1);

% Write actin info into image info cells
for image_index = 1:actin_number_images
    info_struct = struct();
    info_struct.image_type = 'actin';
    info_struct.contrast_ratio = image_contrast_ratios(image_index);
    info_struct.event_number = image_event_numbers(image_index);
    info_struct.line_width = image_line_widths(image_index);
    info_struct.line_type = image_line_types{image_index};
    info_struct.number_lines = image_number_lines(image_index);
    info_cells{image_index} = info_struct;
end
test_info.image_info(current_image_number:current_image_number + actin_number_images) = info_cells;
current_image_number = current_image_number + actin_number_images;

% Prep border info
image_contrast_ratios = choose_evenly(border_number_images, border_contrast_ratios);
image_event_numbers = choose_from_log_range(border_number_image, border_event_number_range);
[image_roughness, image_displacements] = choose_evenly(border_number_images, border_roughness, border_displacement_factor);
info_cells = cell(border_number_images, 1);

% Write border info into image info cells
for image_index = 1:border_number_images
    info_struct = struct();
    info_struct.image_type = 'border';
    info_struct.contrast_ratio = image_contrast_ratios(image_index);
    info_struct.event_number = image_event_numbers(image_index);
    info_struct.roughness = image_roughness(image_index);
    info_struct.displacement = image_displacements(image_index);
    info_cells{image_index} = info_struct;
end
test_info.image_info(current_image_number:current_image_number + border_number_images) = info_cells;

% ---- Create the images for testing and analysis ----

% Initalize cells for holding the images
test_info.ideal_images = cell(total_number_images, 1);
test_info.STORM_images = cell(total_number_images, 1);
test_info.event_data = cell(total_number_images, 1);
test_info.ground_truth_coords = cell(total_number_images, 1);

% Repeat image creation loop
for image_index = 1:total_number_images
    
    % Get and set parameters common to all images
    info = info_struct.image_info{image_index};
    params.number_events_ch1 = info.event_number;
    
    % Create the pdf image for sampling and reference coordinates
    if strcmp(info.image_type, 'region')
        [pdf_map, polygon_coords] = region_pdf_map(params, info.number_vertices, info.contrast_ratio);
        info.ground_truth_coords{image_index} = polygon_coords;
    elseif strcmp(info.image_type, 'dots')
        [pdf_map, dot_coords] = region_pdf_map(params, info.number_dots, info.dot_sizes, info.contrast_ratio);
        info.ground_truth_coords{image_index} = dot_coords;
    elseif strcmp(info.image_type, 'actin')
        [pdf_map, control_points_x, control_points_y] = lines_pdf_map(params, info.number_lines, info.line_width, info.contrast_ratio, info.line_type);
        info.ground_truth_coords{image_index} = cat(3, control_points_x, control_points_y);
    elseif strcmp(info.image_type, 'border')
        [pdf_map, border_coords] = border_pdf_map(params, info.displacement, info.roughness, info.contrast_ratio);
        info.ground_truth_coords{image_index} = border_coords;
    end
    
    % Generate event datapoint using pdf map
    [dataset_coords, ~, STORM_vars] = create_test_data_dv(params); % Uses pdf_map based on params setttings
    info.event_data{image_index} = dataset_coords;
    
    % Create STORM image
    STORM_image = create_test_STORM_images_dv(params, dataset_coords, [], STORM_vars, false, true, true);
    info.STORM_images{image_index} = STORM_image;
    
    % Calculate ideal image
    sigma = params.STORM_precision / params.STORM_pixel_size;
    filter = fspecial('gaussian', ceil(sigma*5), sigma);
    pdf_map = imfilter(pdf_map, filter, 0);
    ideal_image = pdf_map(4:7:end, 4:7:end); % Sample middle pixel of 7x7 field, maintains origin position.
    info.ideal_images{image_index} = ideal_image;
end
end

function [choices, extra_choices] = choose_evenly(number_choices, choice_vector, extra_vector)
% Local function for choosing from a list of choices as evenly as possible,
% with any remainder from a non-evenly divisible number of choices randomly
% assigned. Returned list is randomly permuted.

number_full_sets = floor(number_choices / length(choice_vector));
number_random_choices = mod(number_choices, length(choice_vector));
random_indices = randperm(length(choice_vector), number_random_choices);
choices = repmat(choice_vector(:), number_full_sets, 1);
choices = [choices; choice_vector(random_indices)'];
perm_indices = randperm(number_choices);
choices = choices(perm_indices);
if nargout > 1
    extra_choices = repmat(extra_vector(:), number_full_sets, 1);
    extra_choices = [extra_choices; extra_vector(random_indices)'];
    extra_choices = extra_choices(perm_indices);
end
end

function choices = choose_from_log_range(number_choices, choice_range)
% Local function for choosing a value randomly from a range with a 
% logrithmic distribution. 

log_min = log(choice_range(1));
log_max = log(choice_range(2));
log_choices = rand(number_choices, 1) * (log_max - log_min) + log_min;
choices = exp(log_choices);
end
