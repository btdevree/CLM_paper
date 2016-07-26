function make_test(parameter_structure, output_directory, test_number)
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
%   called "test_archive.tar.gz" containing all the information about the
%   test and a file "test_files.mat" with the information needed by the
%   parametric STORM image testing GUI. 

% Rename parameter structure for convenience
params = parameter_structure;

% Start up a new parallel pool using default settings, if needed.
gcp

% ---- Define the test images -----

% Define test characteristics
test_version = '0.5';

% Change log: 
%   v0.2 - reduce maximum event number to 3e6 to help avoid out of memory errors
%   v0.3 - change parameter ranges for more targeted test
%   v0.4 - adjust parameter ranges for dots and actin, remove quatratic curves
%   v0.5 - adjust parameter ranges for dots and actin, remove cubic curves

% Define important pdf and STORM image parameters just in case the given params don't have this set correctly
% NOTE: If you want to change these, you'll probalby have to change the code below to make sure the ideal image is still correct.
params.ch1_distribution_params{2} = 3; % pdf_map resolution, nanometers, should give a ~73MB map with a 25.6 x 25.6 um field
params.STORM_pixel_size = 21; % works fine for 25 nm STORM precision, ideal image is 1/7th the pdf map (0.024% alaised signal in Fourier domain)

% Regions
region_number_images = 10;
region_contrast_ratios = [0.5; 1; 4]; % background:additional_density
region_event_number_range = [3e2, 1e6; 2e2, 5e5; 1e2, 2e5]; % Multiple rows are for each contrast ratio
region_event_number_focus = [3e2, 1e5; 2e2, 5e4; 1e2, 2e4];
region_number_vertices = [5; 6; 6; 7]; % Don't want people to always know exactly how many points they should be able to find

% Dots
dots_number_images = 15;
dots_sizes = [50; 100; 200]; % nanometers
dots_contrast_ratios = [4, 19; 2, 9; 1, 4]; % Multiple rows are for each dot size
dots_event_number_range = [1e3, 3e6];
dots_event_number_focus = [3e3, 1e6];
dots_number_dots = [8; 9; 9; 10; 10; 11; 11; 12]; % Don't want people to always know exactly how many points they should be able to find

% Actin lines
actin_number_images = 15;
actin_line_widths = [9; 26]; % actin and microtubules, respectivly
actin_contrast_ratios = [5, 19; 3, 9]; % Multiple rows are for each line width size
actin_event_number_range = [1e3, 3e6];
actin_event_number_focus = [3e3, 1e6];
actin_line_types = {'line_segment', 'quadratic'}; % 'cubic' also available, but not useful as it's too hard to match
actin_number_lines = [4; 5; 5; 6]; % Don't want people to always know exactly how many lines they should be able to find

% Border
border_number_images = 10;
border_contrast_ratios = 4;
border_event_number_choices = [3e3; 3e4; 3e5];
border_roughness = [.35, .4, .45, .475, .5, .525, .55, .575, .6, .65, .7, .75];
border_displacement_factor = [.36, .35, .34, .33, .32, .3, .28, .26, .24, .22, .19, .16]; % Apply together with the above roughness factor, not for all combinations possible

% ---- Determine the required parameters for each image in the test ----

% Add up number of images
total_number_images = region_number_images + dots_number_images + actin_number_images + border_number_images;

% Shuffle the RNG
rng('shuffle');

% Create test structure and fill with basic information
test_info = struct();
test_info.test_number = test_number;
test_info.test_version = test_version;
test_info.RNG_state = rng();
test_info.params = params;
test_info.image_info = cell(total_number_images, 1);

% Prep region info
image_contrast_ratios = choose_evenly(region_number_images, region_contrast_ratios);
image_event_numbers = zeros(region_number_images, 1);
for cr_index = 1:size(region_contrast_ratios, 1); % Set seperate event_number ranges for each contrast ratio
    cr_value = region_contrast_ratios(cr_index);
    selection_index_vector = image_contrast_ratios == cr_value;
    num_needed = sum(selection_index_vector, 1);
    event_number_values = round(choose_from_log_range(num_needed, region_event_number_range(cr_index, :), region_event_number_focus(cr_index, :)));
    image_event_numbers(selection_index_vector) = event_number_values;
end
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
image_dot_sizes = choose_evenly(dots_number_images, dots_sizes);
image_contrast_ratios = zeros(dots_number_images, 1);
for size_index = 1:size(dots_sizes, 1); % Set seperate constrast ratio ranges for each size of dot
    size_value = dots_sizes(size_index);
    selection_index_vector = image_dot_sizes == size_value;
    num_needed = sum(selection_index_vector, 1);
    cr_values = choose_evenly(num_needed, dots_contrast_ratios(size_index, :));
    image_contrast_ratios(selection_index_vector) = cr_values;
end
image_event_numbers = round(choose_from_log_range(dots_number_images, dots_event_number_range, dots_event_number_focus));
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
test_info.image_info(current_image_number + 1:current_image_number + dots_number_images) = info_cells;
current_image_number = current_image_number + dots_number_images;

% Prep actin info
image_line_widths = choose_evenly(actin_number_images, actin_line_widths);
image_contrast_ratios = zeros(actin_number_images, 1);
for width_index = 1:size(actin_line_widths, 1); % Set seperate constrast ratio ranges for each line width
    width_value = actin_line_widths(width_index);
    selection_index_vector = image_line_widths == width_value;
    num_needed = sum(selection_index_vector, 1);
    cr_values = choose_evenly(num_needed, actin_contrast_ratios(width_index, :));
    image_contrast_ratios(selection_index_vector) = cr_values;
end
image_event_numbers = round(choose_from_log_range(actin_number_images, actin_event_number_range, actin_event_number_focus));
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
test_info.image_info(current_image_number + 1:current_image_number + actin_number_images) = info_cells;
current_image_number = current_image_number + actin_number_images;

% Prep border info
image_contrast_ratios = choose_evenly(border_number_images, border_contrast_ratios);
image_event_numbers = choose_evenly(border_number_images, border_event_number_choices);
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
test_info.image_info(current_image_number + 1:current_image_number + border_number_images) = info_cells;

% ---- Create the images for testing and analysis ----

% Initalize cells for holding the images
test_info.ideal_images = cell(total_number_images, 1);
test_info.STORM_images = cell(total_number_images, 1);
test_info.event_data = cell(total_number_images, 1);
test_info.ground_truth_coords = cell(total_number_images, 1);

% Initalize structure and cells for holding GUI display images
GUI_info = struct();
GUI_info.STORM_images = cell(total_number_images, 1);

% Repeat image creation loop
fprintf('\nWorking on image number   ');
for image_index = 1:total_number_images
    
    % Feedback to the console
    fprintf('\b\b\b %2u', image_index);
    
    % Get and set parameters common to all images
    info = test_info.image_info{image_index};
    params.number_events_ch1 = info.event_number;
    
    % Create the pdf image for sampling and reference coordinates
    if strcmp(info.image_type, 'region')
        [pdf_map, polygon_coords] = region_pdf_map(params, info.number_vertices, info.contrast_ratio);
        test_info.ground_truth_coords{image_index} = polygon_coords;
    elseif strcmp(info.image_type, 'dots')
        [pdf_map, dot_coords] = dots_pdf_map(params, info.number_dots, info.dot_sizes, info.contrast_ratio);
        test_info.ground_truth_coords{image_index} = dot_coords;
    elseif strcmp(info.image_type, 'actin')
        [pdf_map, control_points_x, control_points_y] = lines_pdf_map(params, info.number_lines, info.line_width, info.contrast_ratio, info.line_type);
        test_info.ground_truth_coords{image_index} = cat(3, control_points_x, control_points_y);
    elseif strcmp(info.image_type, 'border')
        [pdf_map, border_coords] = border_pdf_map(params, info.displacement, info.roughness, info.contrast_ratio);
        test_info.ground_truth_coords{image_index} = border_coords;
    end
    
    % Generate event datapoint using pdf map
    [dataset_coords, ~, STORM_vars] = create_test_data_dv(params); % Uses pdf_map based on params setttings
    test_info.event_data{image_index} = single(dataset_coords); % Single precision is more than accurate enough
    
    % Create STORM image
    STORM_image = create_test_STORM_images_dv(params, dataset_coords, [], STORM_vars, false, true, true);
    test_info.STORM_images{image_index} = single(STORM_image); % Single precision is more than accurate enough for this work, would want doubles for fourier-domain math
    
    % Create display STORM image, there's no reason to keep more precision than 0-255 values.
    min_value = min(STORM_image(:));
    max_value = max(STORM_image(:));
    GUI_info.STORM_images{image_index} = uint8(256 * (STORM_image - min_value) / (max_value - min_value));
    
    % Calculate ideal image
    sigma = params.STORM_precision / params.ch1_distribution_params{2};
    filter = fspecial('gaussian', ceil(sigma*5), sigma);
    pdf_map = imfilter(pdf_map, filter, 'replicate'); 
    ideal_image = pdf_map(4:7:end, 4:7:end); % Sample middle pixel of 7x7 field, maintains origin position.
    test_info.ideal_images{image_index} = ideal_image;
end

% ---- Save files ----
fprintf('\n Image creation done, saving files.');

% Create folder
folder_name = ['Test_', num2str(test_number)];
if ~exist([output_directory, '/', folder_name], 'dir')
    mkdir(output_directory, folder_name);
end

% Copy the rest of the files for GUI program
GUI_info.test_number = test_info.test_number;
GUI_info.image_info = test_info.image_info;
GUI_info.params = test_info.params;
GUI_info.test_version = test_info.test_version;

% Save GUI info
pathname = [output_directory, '/Test_', num2str(test_number),'/test_files.mat'];
save(pathname, 'GUI_info');

% Save archive info
pathname = [output_directory, '/Test_', num2str(test_number),'/test_archive.mat'];
save(pathname, 'test_info');
tarfolder = [output_directory, '/Test_', num2str(test_number)];
tar([tarfolder, '/test_archive.tar.gz'], 'test_archive.mat', tarfolder);

% Delete archive matfile to save space
delete(pathname);

fprintf(' Done. \n');
end

function [choices, extra_choices] = choose_evenly(number_choices, choice_vector, extra_vector)
% Local function for choosing from a list of choices as evenly as possible,
% with any remainder from a non-evenly divisible number of choices randomly
% assigned. Returned list is randomly permuted.

number_full_sets = floor(number_choices / length(choice_vector));
number_random_choices = mod(number_choices, length(choice_vector));
random_indices = randperm(length(choice_vector), number_random_choices);
choices = repmat(choice_vector(:), number_full_sets, 1);
rand_choices = choice_vector(random_indices);
choices = [choices; rand_choices(:)];
perm_indices = randperm(number_choices);
choices = choices(perm_indices);
if nargout > 1
    extra_choices = repmat(extra_vector(:), number_full_sets, 1);
    extra_choices = [extra_choices; extra_vector(random_indices)'];
    extra_choices = extra_choices(perm_indices);
end
end

function choices = choose_from_log_range(number_choices, choice_range, focus_range)
% Local function for choosing a value randomly from a range with a 
% logrithmic distribution. Focus range gives option to concentrate on a
% smaller part of the range.

% Default focus is the same range
if nargin < 3; focus_range = choice_range; end;

% Split the choices into two groups - focus is doubled the rate of the full range
full_number_choices = ceil(number_choices / 2);
focus_number_choices = number_choices - full_number_choices;

% Set choices
log_min = log(choice_range(1));
log_max = log(choice_range(2));
log_choices = rand(full_number_choices, 1) * (log_max - log_min) + log_min;

% Set focus choices
focus_log_min = log(focus_range(1));
focus_log_max = log(focus_range(2));
focus_log_choices = rand(focus_number_choices, 1) * (focus_log_max - focus_log_min) + focus_log_min;

% Mix up choices
all_choices = [log_choices; focus_log_choices];
mixed_choices = all_choices(randperm(length(all_choices(:))));

% Exponentiate and return
choices = exp(mixed_choices);
end
