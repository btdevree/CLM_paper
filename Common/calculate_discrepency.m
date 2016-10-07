function [ discrepency_value, optimized_scale ] = calculate_discrepency( experimental_image, ideal_image, method, optimize_scaling, number_bins)
%CALCULATE_DISCREPENCY Calculate a discrepency measure of the experimental
%   image to the ideal image according the the chosen method.
%
% Inputs:
%   experimental_image: image to be compared with the ideal, possibly 
%       sparse matrix of doubles
%   ideal_image: ideal image used as the expected result, possibly sparse 
%       matrix of doubles
%   method: method used for estimating the discrepency. Choices are:
%       'sum_of_squares': the sum of the sqared differences between the 
%           pixel values.
%       'l2_norm': the sum of the l^2 or Euclidean distances between the 
%           two pixel values, also know as sum of absolute distances. 
%       'discrimination': sum of (a*log(a/b) - a + b) where a is the
%           pixel value of the ideal image and b is the pixel value of the
%           approximation. Images cannot have zeros values.
%       'normalized_variation_of_information': 1 - variation of information
%           (i.e. the mutual information normalized by the joint entropy 
%           calculated from intensity histograms.
%       'normalized_compression_distance': 1 - normalized compression 
%           distance, based on the 7zip LZMA2 algorithm.
%   optimize_scaling: boolean, if true the experimental image is scaled in
%       order to minimize the returned discrepency value. Default = true.
%   number_bins: number of bins to use in for the normalized variation
%       of information method or the normalized_compression_distance. 
%       Integer, ignored for other methods. Default = 100 for NVI, 256 for
%       NCD.
% Outputs:
%   discrepency_value: floating-point double, value of discrepency between 
%       the experimental and ideal image.
%   optimized_scale: optimized scaling factor for the experimental image

% Set defaults
if nargin < 4; optimize_scaling = true; end;
if nargin < 5
    if strcmp(method, 'normalized_variation_of_information')
        number_bins = 100;
    elseif strcmp(method, 'normalized_compression_distance')
        number_bins = 256;
    end
end

% Convert to full matrices, if needed
if issparse(experimental_image)
    experimental_image = full(experimental_image);
end
if issparse(ideal_image)
    ideal_image = full(ideal_image);
end

% Enforce double conversion and normalize images
experimental_image = double(experimental_image);
exp_min = min(experimental_image(:));
exp_max = max(experimental_image(:));
experimental_image = experimental_image - exp_min;
if exp_max ~= 0 % Don't divide by zero!
    experimental_image = experimental_image / exp_max;
end
ideal_image = double(ideal_image);
ideal_min = min(ideal_image(:));
ideal_max = max(ideal_image(:));
ideal_image = ideal_image - ideal_min;
if ideal_max ~= 0 % Don't divide by zero!
    ideal_image = ideal_image / ideal_max;
end

% Start with optimized scale = 1
optimized_scale = 1;

% Optimize scale_factor, if requested
initial_scale = 1;
if optimize_scaling
    
    % Calculate the discrepency as the sum of squared distance
    if strcmp(method, 'sum_of_squares')
        % Optimize with ssq
        f = @(scale_factor)calc_ssq(scale_factor, ideal_image, experimental_image);
    end

    % Calculate the discrepency as the sum of l2 distance
    if strcmp(method, 'l2_norm')
        % Optimize with ssq then l2
        f = @(scale_factor)calc_l2_norm(scale_factor, ideal_image, experimental_image);
    end

    % Calculate the discrepency as the sum of the discrimination
    if strcmp(method, 'discrimination')
        % Calc experimental value
        f = @(scale_factor)calc_discrim(scale_factor, ideal_image, experimental_image);
    end
    
    % Calculate the discrepency with mutual information based Normalized Variation of Information
    if strcmp(method, 'normalized_variation_of_information')
        % Calc discrepency
        f = @(scale_factor)calc_NVI(scale_factor, ideal_image, experimental_image, '2', number_bins);
    end
    
    % Calculate the discrepency with mutual information based Normalized Variation of Information
    if strcmp(method, 'normalized_compression_distance')
        % Calc discrepency
        f = @(scale_factor)calc_NCD(scale_factor, ideal_image, experimental_image, number_bins);
    end
    
    % Run the optimization
    options = optimoptions('fminunc', 'Algorithm', 'quasi-newton', 'Display', 'none');
    if ~strcmp(method, 'sum_of_squares')
        % First optimize on ssq; more stable
        [initial_scale] = fminunc(@(scale_factor)calc_ssq(scale_factor, ideal_image, experimental_image), 1, options);
    end
    % Now optimize on the objective function
    [optimized_scale] = fminunc(f, initial_scale, options);
end

% Get the discrepency for return
% Calculate the discrepency as the sum of squared distance
if strcmp(method, 'sum_of_squares')
    [discrepency_value] = calc_ssq(optimized_scale, ideal_image, experimental_image);
end

% Calculate the discrepency as the sum of l2 distance
if strcmp(method, 'l2_norm')
    [discrepency_value] = calc_l2_norm(optimized_scale, ideal_image, experimental_image);
end

% Calculate the discrepency as the sum of the discrimination
if strcmp(method, 'discrimination')
    [discrepency_value] = calc_discrim(optimized_scale, ideal_image, experimental_image);
end

% Calculate the discrepency with mutual information based Normalized Variation of Information
if strcmp(method, 'normalized_variation_of_information')
    [discrepency_value] = calc_NVI(optimized_scale, ideal_image, experimental_image, '2', NVI_entropy_bins);
end
end

% Local function to calculate the discrimination
function [value] = calc_discrim(scale_factor, A, B)
    sB = scale_factor .* B;
    all_values = A .* log2(A ./ sB) - A + sB;
    value = sum(all_values(:));
end

% Local function to calculate the NVI
function [NVI] = calc_NVI(scale_factor, A, B, log_base, number_bins)
    sB = scale_factor .* B;
    H_A = calc_image_entropy(A, log_base, number_bins);
    H_B = calc_image_entropy(sB, log_base, number_bins);
    H_AB = calc_image_joint_entropy(A, sB, log_base, number_bins);
    I_AB = H_A + H_B - H_AB;
    NVI = 1 - (I_AB/H_AB); % Bounded between 1 and 0; 0 = perfect identity
end

% Local funciton to calculate the SSQ error
function [error] = calc_ssq(scale_factor, A, B)
    squares = (A - scale_factor .* B).^2;
    error = sum(squares(:));
end

% Local funciton to calculate the L2 norm error
function [error] = calc_l2_norm(scale_factor, A, B)
    squares = (A - scale_factor .* B).^2;
    error = sum(sqrt(squares(:)));
end

% Local function to calculate the NCD
function [NCD] = calc_NCD(scale_factor, A, B, number_bins)
    sB = scale_factor .* B;
    AB = [A, sB];
    C_a = calc_image_compression_size(A, number_bins);
    C_b = calc_image_compression_size(sB, number_bins);
    C_ab = calc_image_compression_size(AB, number_bins);  
    NCD = (C_ab - min(C_a, C_b)) / max(C_a, C_b); % Nominally bounded between 1 and 0; 0 = perfect identity, can get to values of 1.1 due to compressor inefficiencies
end
