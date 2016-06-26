function make_test(output_directory, test_number)
%MAKE_TEST Prepare a parametric STORM image test and save it into the 
% specified directory.
%
% Requires functions from Common, Actin_simulator, Border_simulator,
% Dots_simulator, and Region_simulator folders; make sure these are on the
% MATLAB path. 
%   
% Input: 
%   output_directory: string, filepath for saving the test folder
%   test_number; integer, the number to assign to the test
% Output:
%   Test is created and saved to the output directory in a folder called
%   "Test_<test_number>". Inside the folder, there is an archived .mat file
%   called "test_archive.tar.bz" containing all the information about the
%   test and a file "test_files.mat" with the information needed by the
%   parametric STORM image testing GUI. 

% Define test characteristics
test_version = '0.1';
actin_TCI_range = [


% Create folder
mkdir(output_directory, ['Test_', num2str(test_number)]);

% Create test structure and fill with basic information
test_info = struct();
test_info.test_number = test_number;
test_info.test_version = test_version;
test_info.image_info = cell(0);

end

