function [ image, data, parameter_structure ] = part_C_create_image_and_data( parameter_structure, seed )
%PART_C_CREATE_IMAGE_AND_DATA Function for parallel execution of image and
%   data creation 

% Create data and save it
[data, data_ch2, ~, STORMvars] = create_test_data_dv(parameter_structure, seed);

% Create image and save it in cell array, don't use parallel processing
% because the task is already distributed with parfeval in the calling
% function
[image, ~, ~] = create_test_STORM_images_dv(parameter_structure, data, data_ch2, STORMvars, false);

end

