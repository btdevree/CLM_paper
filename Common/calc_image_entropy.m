function [ entropy_value ] = calc_image_entropy( image, log_base, number_bins )
%CALC_IMAGE_ENTROPY Calcualte the entropy of an image using the histogram
%   of the intensity distribution.
%
%   Calculates the Shannon entropy of the intensity distribution.
% 
% Inputs:
%   image: possibly sparse array of floating-point doubles or uint16.
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
if nargin < 3; number_bins = 256; end;
if nargin < 2; log_base = '2'; end;

% Make sure image is a full matrix
if issparse(image)
    image = full(image);
end

% Find the minimum and maximum values of the image
min_value = min(image(:));
max_value = max(image(:));

% Make bin ranges
step = (max_value - min_value) / number_bins;
bin_edges = [min_value:step:max_value];
assert(length(bin_edges) == number_bins + 1); % double check that floating-point errors didn't throw us off

% Get bin counts and add endpoint to last true bin
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

