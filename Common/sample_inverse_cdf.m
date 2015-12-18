function [ indices ] = sample_inverse_cdf(sample_rands, cdf)
%SAMPLE_INVERSE_CDF Convert uniform [0-1] random numbers into samples from
%   a cdf function.
%
%   Calculates the indices as the highest index value where the random 
%       number is greater than or equal to the cdf value at that index.  
%
% Inputs:
%   sample_rands: column vector of floating-point doubles, array of random 
%       numbers between 0 and 1.
%   cdf: column_vector of floating-point doubles, monotonically
%       non-decreasing values between 0 and 1.
% Output:
%   indices: 32 bit integers, the index values of the cdf entry that
%       corrospond to the values of the random samples.

% Initialize indices matrix
indices = zeros(size(sample_rands), 'int32');

% Get cdf max index
max_index = size(cdf, 1);

% Loop through each sample
for rand_index = 1:size(sample_rands, 1)
    rand_value = sample_rands(rand_index);
    
    % Reset the upper and lower index limits
    upper_limit = max_index;
    lower_limit = 1; 
    
    % Search for the corrosponding cdf value by repeated comparison at halfway between the min and max limits
    while upper_limit - lower_limit > 1
        
        % Average the upper and lower limits and get the corrosponding cdf value
        test_index = round((upper_limit + lower_limit)/2);
        test_value = cdf(test_index);
        
        % If the random value is smaller than the test value, move the upper limit to the test index
        if rand_value < test_value
            upper_limit = test_index;
            
        % If the random value is equal to or larger than the test value, move the lower limit to the test index
        else
            lower_limit = test_index;
        end
    end
    
    % Write the lower limit into the output matrix
    indices(rand_index) = int32(lower_limit);
end     
end

