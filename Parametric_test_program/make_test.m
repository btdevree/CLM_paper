function make_test( parmeter_structure, output_directory, test_number)
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

% Shuffle the RNG
rng('shuffle');

% Define test characteristics
test_version = '0.1';

% Regions
region_number_images = 6;
region_SNratios = [1, 4, 10];
region_event_number_range = [3e2, 1e6];
region_number_points = [6];

% Dots
dots_number_images = 15;
dots_SNratios = [1, 4, 10];
dots_sizes = [20, 50, 100, 200, 500]; 
dots_event_number_range = [3e2, 3e6];
dots_number_dots = [10];

% Actin lines
actin_number_images = 18;
actin_SNratios = [1, 4, 10];
actin_event_number_range = [3e2, 1e7];
actin_line_types = {'line_segment', 'cubic', 'quadratic'};
actin_line_width = [9, 26]; % actin and microtubules, respectivly

% Border
border_number_images = 15;
border_SNratios = [1, 4, 10];
border_event_number_range = [1e3, 1e7];
border_roughness = [.35, .45, .55, .65, .75];
border_displacement_factor = [.4, .35, .3, .25, .2]; % Apply together with the above roughness factor, not for all combinations possible

% Add up number of images
total_number_images = region_number_images + dots_number_images + actin_number_images + border_number_images;

% Rename parameter structure for convenience
params = parameter_struct;

% Create folder
mkdir(output_directory, ['Test_', num2str(test_number)]);

% Create test structure and fill with basic information
test_info = struct();
test_info.test_number = test_number;
test_info.test_version = test_version;
test_info.image_info = cell(total_number_images, 1);


end

