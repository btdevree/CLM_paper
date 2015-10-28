function [RGB_image] = convert_STORM_image_to_RGB(channel_pd_data, channel_color, number_binned_pixels, normalization_method, percent_saturated, channel_multipliers)
%CONVERT_STORM_PDF_TO_RGB Converts a set of per-channel probability density images to an 8bit RGB image for display.
%  
%   Accecpts a set of (possibly sparse) per-channel matrices and creates
%   an 8 bit RGB image for display purposes. 
%   
%   Inputs:
%   channel_pd_data: cell array of possibly sparse, registered matrices 
%       representing images of the probability distribution of the channel 
%       events.
%   channel_color: cell array of RBG color triplets that specify the color
%       of each channel in the final image.
%   number_binned_pixels: optional binning of the origonal pixels to make a
%       smaller image. Default = 1 (no binning). Extra pixels are truncated
%       from the bottom row and/or left column.
%   normalization_method: method used to normalize the image. Choices are:
%       'total' - image scaled to the overall data range across all
%       channels. Use with channel multipliers to manualy change the 
%       relative strength of each channel. Default method. 
%       'channel_max' - each channel is scaled according to its maximum
%       value.
%       'channel_nonzero_median' - each channel is scaled according to its median
%       value.
%       'channel_nonzero_mean' - each channel is scaled according to its mean 
%       value.
%   percent_saturated: percentage of saturated color in pixels of the final 
%       RGB image. Default = 0.
%   channel_multipliers: cell array of per-channel multiplier values to
%       manually change the channel strengths with the 'total'
%       normalization method. Default = 1 for each channel.
%   Outputs:
%   RGB_image: 8 bit RGB image, possibly downsampled with binning. Use for
%       display purposes, NOT for math.

% Set defaults
if nargin < 6,
    channel_multipliers = cell(size(channel_color));
    channel_multipliers(:) = {1};
end
if nargin < 5; percent_saturated = 0; end
if nargin < 4; normalization_method = 'total'; end
if nargin < 3; number_binned_pixels = 1; end

% Initialize a new set of sparse RGB matrices
[num_rows, num_columns] = size(channel_pd_data{1});
red_sparse = sparse(num_rows, num_columns);
green_sparse = sparse(num_rows, num_columns);
blue_sparse = sparse(num_rows, num_columns);

% Scale channels according to method
if strcmp(normalization_method, 'total'), % Total method
    
    % Apply the multiplier to the channel 
    for channel_index = 1:length(channel_pd_data),
        channel_pd_data{channel_index} = channel_pd_data{channel_index}.*channel_multipliers{channel_index};
    end
    
else % Per-channel methods
    
    % Determine divisors according to selected method
    divisors = cell(size(channel_color));
    
    % Per-channel maximum
    if strcmp(normalization_method, 'channel_max'), 
        
        % Scale the values in the channel to maximum = 1
        for channel_index = 1:length(channel_pd_data),
           divisors{channel_index} = max(channel_pd_data{channel_index}(:));
        end
        
    % Per-channel mean
    elseif strcmp(normalization_method, 'channel_nonzero_mean'),
        
        % Scale the values in the channel to mean = 1
        for channel_index = 1:length(channel_pd_data),
           divisors{channel_index} = mean(channel_pd_data{channel_index}(channel_pd_data{channel_index} > 0));
        end
        
    % Per-channel median
    elseif strcmp(normalization_method, 'channel_nonzero_median'),
        
        % Scale the values in the channel to median = 1
        for channel_index = 1:length(channel_pd_data),
           divisors{channel_index} = median(channel_pd_data{channel_index}(channel_pd_data{channel_index} > 0));
        end
        
    % No valid method
    else
        error('No valid normalization method given');
    end
    
    % Apply divisors
    for channel_index = 1:length(channel_pd_data),
        channel_pd_data{channel_index} = channel_pd_data{channel_index}./divisors{channel_index};
    end
end

% Multiply each channel by the desired color and add to the RGB matrices
for channel_index = 1:length(channel_pd_data),
    red_sparse = red_sparse + channel_pd_data{channel_index}.*channel_color{channel_index}(1);
    green_sparse = green_sparse + channel_pd_data{channel_index}.*channel_color{channel_index}(2);
    blue_sparse = blue_sparse + channel_pd_data{channel_index}.*channel_color{channel_index}(3);
end

% Avoid binning math if not requested
% If no binning is requested, just create a RGB array by filling out sparse arrays
if number_binned_pixels == 1,
    [num_rows, num_columns] = size(red_sparse);
    RGB_double = zeros(num_rows, num_columns, 3);
    RGB_double(:, :, 1) = full(red_sparse);
    RGB_double(:, :, 2) = full(green_sparse);
    RGB_double(:, :, 3) = full(blue_sparse);

% Otherwise, bin the matrices    
else   
    % Calculate the number of rows or columns that will be removed before binning
    [num_rows, num_columns] = size(red_sparse);
    num_rows_to_delete = mod(num_rows, number_binned_pixels);
    num_columns_to_delete = mod(num_columns, number_binned_pixels);

    % Initialize binned RGB double matrix
    num_rows_final = floor(num_rows/number_binned_pixels);
    num_columns_final = floor(num_columns/number_binned_pixels);
    RGB_double = zeros(num_rows_final, num_columns_final, 3);

    % Bin the matrices one-by-one to avoid having lot of full-size copies around

    % Cut off extra rows and columns and make a full red matrix
    red_full = full(red_sparse(1:end-num_rows_to_delete, 1:end-num_columns_to_delete));

    % Bin red and clear memory space
    RGB_double(:, :, 1) = bin_matrix(red_full, number_binned_pixels);
    clear red_full

    % Cut off extra rows and columns and make a full red matrix
    green_full = full(green_sparse(1:end-num_rows_to_delete, 1:end-num_columns_to_delete));

    % Bin green and clear memory space
    RGB_double(:, :, 2) = bin_matrix(green_full, number_binned_pixels);
    clear green_full

    % Cut off extra rows and columns and make a full red matrix
    blue_full = full(blue_sparse(1:end-num_rows_to_delete, 1:end-num_columns_to_delete));

    % Bin blue and clear memory space
    RGB_double(:, :, 3) = bin_matrix(blue_full, number_binned_pixels);
    clear blue_full
end

% Scale the RGB_double matrix between 0 and 255

% Avoid sorting the RGB matrix if no saturated pixels are wanted
if percent_saturated == 0,
    max_value = max(RGB_double(:));
    RGB_double = 255*(RGB_double/max_value);

% Otherwise, find the requested percentile value and scale to that value.
else
    scale_value = prctile(RGB_double(:), 100-percent_saturated);
    RGB_double = 255*(RGB_double/scale_value);
end

% Convert to an 8 bit matrix. 
RGB_image = uint8(RGB_double);
end
