function [ dilated_BWmap ] = standard_sequence_dilate( initial_BWmap, number_cycles )
% STANDARD_SEQUENCE_DILATE

% Dilates the image with a repeating sequence of square-, cross-, and 
% square-shaped kernels in order to help preserve the general shape of 
% the regions. Returns a new array with the dilated output.

% Prepare kernels
cross_kernel = strel('diamond', 1);
square_kernel = strel('square', 3);

% Divide the number of cycles by 3 to find the number of times to loop through the three-cycle sequence.
number_of_full_sequences = floor(number_cycles/3);
number_of_extra_cycles = mod(number_cycles, 3);

% Copy initial map into a working image
working_image = initial_BWmap;

% Dilate sequentually with square-, cross-, and square-shaped kernels for each three-cycle sequence
for index = 1:number_of_full_sequences
    working_image = imdilate(working_image, square_kernel);
    working_image = imdilate(working_image, cross_kernel);
    working_image = imdilate(working_image, square_kernel);
end

% Dilate once with a square kernel if remainder = 1
if number_of_extra_cycles == 1;
    working_image = imdilate(working_image, square_kernel);
end

% Dilate once with a square kernel and once with a cross kernel if remainder = 2
if number_of_extra_cycles == 2;
    working_image = imdilate(working_image, square_kernel);
    working_image = imdilate(working_image, cross_kernel);
end

% Copy working image into final returned product
dilated_BWmap = working_image;
end

