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
params = parameter_struct;

% ---- Define the test images -----

% Define test characteristics
test_version = '0.1';

% Regions
region_number_images = 6;
region_contrast_ratios = [0.2, 1, 4]; % background:additional_density
region_event_number_range = [3e2, 1e6];
region_number_points = [6];

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
image_event_numbers = choose_from_log_range(region_number_image, region_event_number_range);
image_number_points = choose_evenly(region_number_images, region_number_points);
info_cells = cell(region_number_images, 1);

% Write region info into image info cells
for image_index = 1:region_number_images
    info_struct = struct();
    info_struct.image_type = 'region';
    info_struct.contrast_ratio = image_contrast_ratios(image_index);
    info_struct.event_number = image_event_numbers(image_index);
    info_struct.region_points = image_number_points(image_index);
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
current_image_number = current_image_number + actin_number_images;

% Prep actin info
image_contrast_ratios = choose_evenly(actin_number_images, actin_contrast_ratios);
image_event_numbers = choose_from_log_range(actin_number_image, actin_event_number_range);
image_line_widths = choose_evenly(actin_number_images, actin_line_widths);
image_line_types = choose_evenly(actin_number_images, actin_line_types);
info_cells = cell(region_number_images, 1);

% Write actin info into image info cells
for image_index = 1:actin_number_images
    info_struct = struct();
    info_struct.image_type = 'actin';
    info_struct.contrast_ratio = image_contrast_ratios(image_index);
    info_struct.event_number = image_event_numbers(image_index);
    info_struct.line_width = image_line_widths(image_index);
    info_struct.line_type = image_line_types(image_index);
    info_cells{image_index} = info_struct;
end
test_info.image_info(current_image_number:current_image_number + actin_number_images) = info_cells;
current_image_number = current_image_number + actin_number_images;

% Prep border info
image_contrast_ratios = choose_evenly(border_number_images, border_contrast_ratios);
image_event_numbers = choose_from_log_range(border_number_image, border_event_number_range);
[image_roughness, image_displacements] = choose_evenly(border_number_images, border_roughness, border_displacement_factor);
info_cells = cell(region_number_images, 1);

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

test_info.ideal_images = cell(total_number_images, 1);
test_info.STORM_images = cell(total_number_images, 1);
test_info.event_data = cell(total_number_images, 1);

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
