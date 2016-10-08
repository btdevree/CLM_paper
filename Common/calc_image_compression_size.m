function [ file_size ] = calc_image_compression_size(image, number_bins)
%CALC_IMAGE_COMPRESSION_SIZE Calcualte the file size of a compressed image 
%   using the LZMA2 algorithm from p7zip.
%
%   Works on uint16 images for enhansed speed and accuracy of compression.
% 
% Inputs:
%   image: possibly sparse array of floating-point or uint16 values.
%   number_bins: number of discrete values used when writing the image 
%       files to be compressed. Default = []; no changes to the number of
%       values.
% Outputs:
%   file_size: floating-point double value, given in bytes

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
    image = uint16(round(image));

% Convert floating point images
elseif ~isa(image, 'uint16');
    
    % Convert to double before manipulating
    image = double(image);
    image = image - min(image(:));
    image = image / max(image(:));
    image = image * (2.^16 - 1);
    image = uint16(round(image));
end

% We assert that the image is now of uint16 type
assert(isa(image, 'uint16'));

% Create file to compress
filename = tempname;
fileID = fopen(filename, 'w');
fwrite(fileID, image, 'uint16');
fclose(fileID);

% Call p7zip to compress file
command = ['p7zip ', filename];
[status, cmdout] = system(command);

% Check status, print message if not 0
if ~status == 0
    cmdout
end

% Check file size after compression
file_info = dir([filename, '.7z']);
file_size = file_info.bytes;

% Delete file
delete([filename, '.7z']);
end

