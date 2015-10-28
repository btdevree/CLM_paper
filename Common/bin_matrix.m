function [M] = bin_matrix(M, number_binned_pixels)
%BIN_MATRIX Bins values in a matrix. 
%   Matrix columns and rows must be evenly divisible by the number of 
%   binned pixels.

% Determine final row and column sizes
number_rows = size(M, 1)/number_binned_pixels;
number_columns = size(M, 2)/number_binned_pixels;

% Break up columns into bin-size peices arranged in a row
M = reshape(M, number_binned_pixels,[]);

% Sum columns
M = sum(M, 1);

% Re-assemble binned columns and transpose
M = reshape(M, number_rows, []).'; 

% Break up columns (previously rows) into bin-size peices arranged in a row
M = reshape(M, number_binned_pixels, []);

% Sum columns
M = sum(M, 1);

% Re-assemble binned columns and transpose into original orientation
M = reshape(M, number_columns, []).'; 
end

