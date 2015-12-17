function [ xy_samples ] = sample_2D_pdf( number_samples, pdf_image, resolution )
%SAMPLE_2D_PDF Samples 2D random variables from an arbitrary 2D pdf image.
%
% Uses linear interpolation to determine sampling rates between the pixel
%   centers. Assumes a Cartesian coordinate system with the origin in the 
%   bottom left corner of the bottom left pixel. Algorithm employs inverse 
%   sampling of a greyscale dilated image to sample points at a slightly 
%   higher rate than needed and then corrects the rate with rejection 
%   sampling. 
%
% Inputs:
%   number_samples: number of samples to return.
%   pdf_image: array of floating-point doubles, a discrete pdf to be 
%       sampled from. Does not have to be normalized.
%   resolution: optional, the number of measurment units per pixel. 
%       Default = 1.
% Outputs:
%   xy_samples: n by 2 array to floating point doubles. Given in the
%       measurement units supplied in the resolution parameter.

% Set defaults
if nargin < 3; resolution = 1; end;

% Get dimensions
image_n = size(pdf_image, 1);
image_m = size(pdf_image, 2);

% Create meshgrid variable for interpolation
x_vec = resolution * [0.5: 1: image_m - 0.5].';
y_vec = resolution * [image_n - 0.5: -1: 0.5].';
[ x_mesh, y_mesh ] = meshgrid(x_vec, y_vec);

% Create a dilated image
se = strel('square',3);
enclosing_pdf = imdilate(pdf_image, se);

% Average dilated and original pdfs to get the maximum sampling rate that will be needed for each pixel (given linear interpolation)
enclosing_pdf = (enclosing_pdf + pdf_image) / 2;

% Calculate the expected rejection rate in rejected points per created point
acceptance_rate = sum(pdf_image(:)) / sum(enclosing_pdf(:));
rejection_rate = 1 - acceptance_rate;

% Calculate the cumulative sum of the enclosing pdf image unraveled as a 1D array
enclosing_cdf = [0; cumsum(enclosing_pdf(:))];
enclosing_cdf = enclosing_cdf / enclosing_cdf(end); % normalize to range [0,1]

% Initialize variables
xy_samples = zeros(number_samples, 2);
num_completed_samples = 0;

% Attempt to make the required number of samples. Repeat if this didn't give enough.
while num_completed_samples < number_samples
   
    % Calculate the required number of trials to end up with enough accecpted samples
    % Use a 99.5% single-tailed CI to the Bernoulli process calculated with the Wald approximation
    % Solution is quadratic in form with n = x^2
    a = acceptance_rate;
    b = -2.576 * sqrt(acceptance_rate * rejection_rate);
    c = -(number_samples - num_completed_samples);
    num_new_samples = ceil(((-b + sqrt(b.^2 - 4*a*c)) / (2*a)).^2);
      
    % Get the pixel indices using inverse sampling
    new_rands = rand(num_new_samples);
    [~, new_pixel_indices] = histc(new_rands, enclosing_cdf);
    
    % The last two bins of the histogram should be added together (default for newer histcounts function)
    new_pixel_indices(new_pixel_indices == length(enclosing_cdf)) = length(enclosing_cdf) - 1;
    
    % Calculate the x coordinate
    new_rands = rand(num_new_samples);
    new_x = resolution * (floor((new_pixel_indices - 1) / image_n) + new_rands);
    
    % Calculate the y coordinate
    new_rands = rand(num_new_samples);
    new_y = resolution * (floor(mod(new_pixel_indices - 1, image_n)) + new_rands);
    
    % Interpolate the pdf map with the new points
    new_pdf_values = interp2(x_mesh, y_mesh, pdf_image, new_x, new_y);
    
    % Calculate the accecptance threshold for each point
    acceptance_threshold = new_pdf_values ./ enclosing_pdf(new_pixel_indices);
    
    % Decide whether to accecpt or reject the new point
    new_rands = rand(num_new_samples);
    selection_booleans = new_rands <= acceptance_threshold;
    num_accepted = sum(selection_booleans(:));
    
    % Put new, accepted points into the xy_samples array
    % Most of the time, it should work out that we have more samples accecpted than we need to fill the output array
    if num_accepted >= number_samples - num_completed_samples
        
        % Collect the new accepted points
        new_points = zeros(num_accepted, 2);
        new_points(:, 1) = new_x(selection_booleans);
        new_points(:, 2) = new_y(selection_booleans);
                
        % Copy the needed points from the beginning of the list of accepted points
        num_points_needed = number_samples - num_completed_samples;
        xy_samples(1 + num_completed_samples:number_samples, :) = new_points(1:num_points_needed, :);
        
        % Keep track of the number of points added
        num_completed_samples = num_completed_samples + num_points_needed;
    
    % Sometimes, there will not be enough samples
    elseif num_accepted < number_samples - num_completed_samples
        
        % Collect the new accepted points
        new_points = zeros(num_accepted, 2);
        new_points(:, 1) = new_x(selection_booleans);
        new_points(:, 2) = new_y(selection_booleans);
        
        % Copy the accepted points to the output list
        xy_samples(1 + num_completed_samples:num_accepted + num_completed_samples, :) = new_points;
        
        % Keep track of the number of points added
        num_completed_samples = num_completed_samples + num_accepted;
    end
end
end

