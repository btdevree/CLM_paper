function [ discrepency_value ] = calculate_discrepency( experimental_image, ideal_image, method, NVI_entropy_bins)
%CALCULATE_SURPRISAL Calculate a discrepency measure of the experimental
%   image to the ideal image according the the chosen method.
%
%   Experimental discrepency measures are normalized to the discrepency 
%   between the ideal image and a zero-valued image so that 0 = perfect
%   agreement between the ideal and experimental image and 1 = the same
%   dissimilarity as the discrepency between ideal and zero images. This
%   means that an experimental image could have a negative discrepency
%   (less similar than even a zero image), but probably shouldn't under 
%   normal imaging conditions. 
%
% Inputs:
%   experimental_image: image to be compared with the ideal, possibly 
%       sparse matrix of doubles
%   ideal_image: ideal image used as the expected result, possibly sparse 
%       matrix of doubles
%   method: method used for estimating the discrepency. Choices are:
%       'sum_of_squares': the sum of the sqared differences between the 
%           pixel values.
%       'discrimination': sum of (a*log(a/b) - a + b) where a is the
%           pixel value of the ideal image and b is the pixel value of the
%           approximation. Images cannot have zeros values.
%       'normalized_variation_of_information': 1 - variation of information
%           (i.e. the mutual information normalized by the joint entropy 
%           calculated from intensity histograms.
%   NVI_entropy_bins: number of bins to use in for the normalized variation
%       of information method. Integer, ignored for ssq or discrimination 
%       methods. Default = 100
% Outputs:
%   discrepency_value: floating-point double, normalized value of 
%       discrpency between the experimental and ideal image.

% Set defaults
if nargin < 4; NVI_entropy_bins = 100; end;

% Convert to full matrices, if needed
if issparse(experimental_image)
    experimental_image = full(experimental_image);
end
if issparse(ideal_image)
    ideal_image = full(ideal_image);
end

% Enforce double conversion and normalize images
experimental_image = double(experimental_image);
experimental_image = experimental_image - min(experimental_image(:));
if max(experimental_image(:)) > 0 % Don't divide by zero!
    experimental_image = experimental_image / max(experimental_image(:));
end
ideal_image = double(ideal_image);
ideal_image = ideal_image - min(ideal_image(:));
if max(ideal_image(:)) > 0 % Don't divide by zero!
    ideal_image = ideal_image / max(ideal_image(:));
end

% Calculate the discrepency as the sum of squared distance
if strcmp(method, 'sum_of_squares')
    
    % Calc zero image ssq value
    zero_image_squares = (ideal_image - zeros(size(experimental_image))).^2;
    zero_image_ssq = sum(zero_image_squares(:));
    
    % Calc experimental ssq value
    experimental_squares = (ideal_image - experimental_image).^2;
    experimental_ssq = sum(experimental_squares(:));
    
    % Get the normalized discrepency
    discrepency_value = experimental_ssq / zero_image_ssq;
end

% Calculate the discrepency as the sum of the discrimination
if strcmp(method, 'discrimination')
        
    % Calc zero image value
    zero_image_discrim = calc_discrim(ideal_image, zeros(size(experimental_image)));
   
    % Calc experimental value
    experimental_discrim = calc_discrim(ideal_image, experimental_image);
    
    % Get the normalized discrepency
    discrepency_value = experimental_discrim / zero_image_discrim;
end

% Calculate the discrepency with mutual information based Normalized Variation of Information
if strcmp(method, 'normalized_variation_of_information')
    
    % Calc discrepency
    discrepency_value = calc_NVI(ideal_image, experimental_image, '2', NVI_entropy_bins);
  
end
end

% Local function to calculate the discrimination
function [value] = calc_discrim(A, B)
    all_values = A .* log2(A ./ B) - A + B;
    value = sum(all_values(:));
end

% Local function to calculate the NVI
function [NVI] = calc_NVI(A, B, log_base, number_bins)
    H_A = calc_image_entropy(A, log_base, number_bins);
    H_B = calc_image_entropy(B, log_base, number_bins);
    H_AB = calc_image_joint_entropy(A, B, log_base, number_bins);
    I_AB = H_A + H_B - H_AB;
    NVI = 1 - (I_AB/H_AB); % Bounded between 1 and 0; 0 = perfect identity
end


