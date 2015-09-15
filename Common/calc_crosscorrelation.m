function [crosscorrelation] = calc_crosscorrelation(image1, image2, max_radius, binary_mask, twoD_window)
%CALC_CROSSCORRELATION Calculates the crosscorrelation between two STORM 
%   images.
%  
%   Calculates the normalized crosscorrelation of the STORM images.
%   Inputs:
%   image1, image2: STORM image given as full or sparse floating-point 
%       double matrices. Assumes these images are the same size as each 
%       other and as the binary mask.
%   max_radius: maximum radius between the center of the correlation and 
%       the image edge for the returned normalized correlation image. Be 
%       sure that the requested distances are actually achievable given the
%       size of the mask or else you will divide by zero. Given in pixels.
%   binary_mask: Binary image that marks the region of interest. Optional,
%       default mask is the whole imaged region.
%   twoD_window: A 2D window function that will be used to filter the
%       binary map of the region of interest in order to prevent spectral
%       leakage from the application of DFTs to the masking region. Optional, 
%       Default is a very simple kernel that makes a triangular window.
%   Outputs:
%   normalized_xcor: Centered crosscorrelation image normalized so that a 
%       random distribution is equal to 1.


% Set defaults
if nargin < 5,
    twoD_window = [1/12, 1/6, 1/12;
                    1/6,   0,  1/6;
                   1/12, 1/6, 1/12];
end
if nargin < 4,
    binary_mask = true(size(image1));
end

% Cut out the smallest possible rectangle that contains the masked regions
props = regionprops(binary_mask, 'BoundingBox'); % N.B. returned box is in stupid matlab image coordinates
bound_box = props.BoundingBox;
row_start = ceil(bound_box(1)); % imcrop will not accecpt sparse matrices, so we have to convert to normal indexes
row_end = row_start + bound_box(3)-1;
column_start = ceil(bound_box(2));
column_end = column_start + bound_box(4)-1;
image1_roi = image1(column_start:column_end, row_start:row_end);
image2_roi = image2(column_start:column_end, row_start:row_end);
mask_roi = binary_mask(column_start:column_end, row_start:row_end);

% Clear the images from memory, since they might take up a lot of space
clear image1 image2 binary_mask;

% Zero pad images to prevent circular convolution
added_zeros = max_radius + (length(twoD_window) - 1) / 2;
number_rows = size(mask_roi, 1) + added_zeros;
number_columns = size(mask_roi, 2) + added_zeros; 
image1_padded = zeros(number_rows, number_columns);
image2_padded = zeros(number_rows, number_columns);
mask_padded = false(number_rows, number_columns);

% Keep this code around for comparing doubling the DFT size to only minimal padding - integrate into options later
% image1_padded = zeros(size(image1_roi, 1)*2+1, size(image1_roi, 2)*2+1);
% image2_padded = zeros(size(image2_roi, 1)*2+1, size(image2_roi, 2)*2+1);
% mask_padded = false(size(mask_roi, 1)*2+1, size(mask_roi, 2)*2+1);

% Put the image values into the zero padded matrices. Convert sparse matrices to full if relevent
mask_padded(1:size(mask_roi, 1), 1:size(mask_roi, 2)) = mask_roi;
clear mask_roi;
if issparse(image1_roi)
    image1_roi = full(image1_roi);
end
image1_padded(1:size(image1_roi, 1), 1:size(image1_roi, 2)) = image1_roi;
clear image1_roi;
if issparse(image2_roi)
    image2_roi = full(image2_roi);
end
image2_padded(1:size(image2_roi, 1), 1:size(image2_roi, 2)) = image2_roi;
clear image2_roi

% Filter the mask and normalize
mask_padded = imfilter(double(mask_padded), twoD_window, 'circular');
mask_padded = mask_padded/max(mask_padded(:));

% Apply the mask to the images
image1_padded = image1_padded.*mask_padded;
image2_padded = image2_padded.*mask_padded;

% Calculate the normalized correlations
crosscorrelation = calc_normalized_correlation(image1_padded, image2_padded, mask_padded, max_radius);
end