function [ entropy_value ] = calc_image_joint_entropy( image_A, image_B, log_base, number_bins )
%CALC_IMAGE_JOINT_ENTROPY Calcualte the joint entropy of two images using 
%   the joint histogram of the intensity distributions.
%
%   Calculates the Shannon entropy of the joint intensity distributions.
% 
% Inputs:
%   image_A/image_B: possibly sparse arrays of floating-point doubles or 
%       uint16.
%   log_base: base of the logarithm used to calculate the Shannon entropy.
%       options are: 
%           '2' - entropy given in bits (default)
%           'e' - entropy given in nats
%           '10' - entropy given in Hartleys 
%   number_bins: number of bins used to form the histogram. Any number may 
%       be used, but default is 256 (same as calculated by  conversion to 
%       uint8 by MATLAB entropy() function).
% Outputs:
%   entropy_value: floating-point double value, given in units determined
%       by the log_base input

% Set defaults
if nargin < 4; number_bins = 256; end;
if nargin < 3; log_base = '2'; end;

% Make sure images are full matrices
if issparse(image_A)
    image_A = full(image_A);
end
if issparse(image_B)
    image_B = full(image_B);
end

% Find the minimum and maximum values of the images
min_value_A = min(image_A(:));
max_value_A = max(image_A(:));
min_value_B = min(image_B(:));
max_value_B = max(image_B(:));

% Get step sizes
step_A = (max_value_A - min_value_A) / number_bins;
step_B = (max_value_B - min_value_B) / number_bins;

% Divide images by the step and take the integer part
% Add 1 to make the indices on MATLAB-style 1-based indexing
if step_A ~= 0
    integer_A = floor((image_A - min_value_A) ./ step_A) + 1;
else % If step == 0, it means the image is flat. We'll assume that we want a flat, zero-valued image.
    integer_A = ones(size(image_A));
end
if step_B ~= 0
    integer_B = floor((image_B - min_value_B) ./ step_B) + 1;
else % If step == 0, it means the image is flat. We'll assume that we want a flat, zero-valued image.
    integer_B = ones(size(image_B));
end

% Count the number of each 2D integer pair
counts = accumarray([integer_A(:), integer_B(:)], 1, [number_bins + 1, number_bins + 1]);

% Add the last endpoint to the previous row or column
counts = [counts(1:end-2, :); counts(end-1, :) + counts(end, :)];
counts = [counts(:, 1:end-2), counts(:, end-1) + counts(:, end)];

% Normalize to a pdf
total_counts = sum(counts(:));
pdf = counts ./ total_counts;

% Calculate entropy
if strcmp(log_base, '2')
    entropy_values = -(pdf .* log2(pdf));
end
if strcmp(log_base, 'e')
    entropy_values = -(pdf .* log(pdf));
end
if strcmp(log_base, '10')
    entropy_values = -(pdf .* log10(pdf));
end
entropy_value = nansum(entropy_values(:));

end

