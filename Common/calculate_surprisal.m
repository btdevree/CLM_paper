function [ surprisal_value ] = calculate_surprisal( experimental_image, ideal_image, method )
%CALCULATE_SURPRISAL Calculate the surprisal, or self-information of the 
%   experimental image according the the chosen method.

% Inputs:
%   experimental_image: image to be compared with the ideal, possibly 
%       sparse matrix of doubles
%   ideal_image: ideal image used as the expected result, possibly sparse 
%       matrix of doubles
%   method: method used for estimating the discrepency. Choices are
%       'Entropy' (real-valued functions) and 'Discrepency' (positive real-
%       valued functions). FOR NOW: 'squared_distance'
%
% Outputs:
%   surprisal_value: self-information or surprisal value, given in bits.

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

% Calculate the discrepency as the squared distance
if strcmp(method, 'squared_distance')
    surprisal_image = (experimental_image - ideal_image).^2;
    surprisal_value = sum(surprisal_image(:));
end
end

