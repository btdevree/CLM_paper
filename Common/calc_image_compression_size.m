function [ file_size ] = calc_image_compression_size(image, number_bins, command_flag)
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
%   command_flag: flag value telling function which commands to use. 0 =
%       use "p7zip" simple command from installed Ubuntu package p7zip. 1 = 
%       use locally installed 7za program, will probably need to custumize 
%       for each machine. Default = 0;
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
if command_flag == 0
    command = ['~/p7zip ', filename];
    [status, cmdout] = system(command);
elseif command_flag == 1
    command = ['~/p7zip_16.02/bin/7za a -sdel ', filename, ' ', filename, '.7z'];
    [status, cmdout] = system(command);
end

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

