function normalized_corr_image = calc_normalized_correlation(image1, image2, mask, max_radius)
% CALC_NORMALIZED_CORRELATION Calculates the 2D correlation of the two
%   images normalized by the given mask.
%
%   Calculate the normalized correlation between image1 and image 2
%   according to the normalization region given in the mask. 
%   Inputs:
%   image1 and image2: floating-point double image matrices. Assumed to be
%       zero padded and have the mask applied.
%   mask: floating-point double matrix with intensity range from 1 to 0. 
%       Assumed to be zero padded.
%   max_radius: maximum radius between the center of the correlation and 
%       the image edge for the returned normalized correlation image. Be 
%       sure that the requested distances are actually achievable given the
%       size of the mask or else you will divide by zero. Given in pixels.

% Calculate average intensities 
mask_sum = sum(mask(:));
aveval1 = sum(image1(:))/mask_sum;
aveval2 = sum(image2(:))/mask_sum;

% Calculate normalization factor image
norm_image_full = fftshift(ifft2(abs(fft2(mask)).^2));

% Cut out the requested region
center_pixel_column_index = ceil(size(norm_image_full, 1)/2);
center_pixel_row_index = ceil(size(norm_image_full, 2)/2); 
column_start = center_pixel_column_index - max_radius;
column_end = center_pixel_column_index + max_radius;
row_start = center_pixel_row_index - max_radius;
row_end = center_pixel_row_index + max_radius;
norm_image = norm_image_full(column_start:column_end, row_start:row_end);

% Clear big normalization image and mask as they might be large
clear mask norm_image_full

% Calculate 2D correlation
corr_raw_full = fftshift(ifft2(fft2(image1).*conj(fft2(image2))));
corr_raw = corr_raw_full(column_start:column_end, row_start:row_end);

% Clear images and big correlation as they might be large
clear image1 image2 corr_raw_full

% Normalize the final correlation image
normalized_corr_image = real((corr_raw./norm_image)/(aveval1*aveval2));
end

