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

% Make bin ranges for use with imquantize (negative valued and reversed)
step_A = (max_value_A - min_value_A) / number_bins;
bin_edges_A = [-max_value_A:step_A:-min_value_A];
assert(length(bin_edges_A) == number_bins + 1); % double check that floating-point errors didn't throw us off
step_B = (max_value_B - min_value_B) / number_bins;
bin_edges_B = [-min_value_B:step_B:-max_value_B];
assert(length(bin_edges_B) == number_bins + 1);

% Get bin counts and add endpoint to last true bin
discrete_A = imquantize(-image_A, bin_edges_A); % level 1 = endpoint at max_value, level N = full bin  from min_value to level (N-1) cutoff

bincounts = histc(image(:), bin_edges);
bincounts = [bincounts(1:end-2); bincounts(end-1) + bincounts(end)];

% Normalize to a pdf
total_counts = sum(bincounts(:));
pdf = bincounts ./ total_counts;

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

