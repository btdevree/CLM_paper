function [ file_size ] = calc_image_compression_size(image, number_bins)
%CALC_IMAGE_COMPRESSION_SIZE Calcualte the file size of a compressed image 
%   using the LZMA2 algorithm from p7zip.
% 
% Inputs:
%   image: possibly sparse array of floating-point doubles or uint16.
%   number_bins: number of discrete values used when writing the image 
%       files to be compressed. Default = []; no changes to the number of
%       values.
% Outputs:
%   entropy_value: floating-point double value, given in units determined
%       by the log_base input

% Set defaults
if nargin < 2; number_bins = []; end;

% Make sure image is a full matrix
if issparse(image)
    image = full(image);
end

% Make binned image, if requested 
if ~isempty(number_bins)

    % Convert the image to integer levels
    image = image - min(image(:));
    image = image / max(image(:));
    image = image * (number_bins - 1);
    image = round(image);
end



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

