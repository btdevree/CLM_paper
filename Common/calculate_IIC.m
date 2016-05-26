function [info_improvement_data] = calculate_IIC(params, dataset, fraction_vector, number_pseudoreplicates, discrepency_method, optimize_flag)
%CALCULATE_IIC Calculate the information improvement curve for a STORM
% image generated with the given parameters

% Inputs: 
%   params: parameters structure for image generation
%   dataset: n by 2 matrix of event centers used to create the images 
%   fraction_vector: column vector of the fractions of data used, include 
%       full and zero fractions.
%   number_pseudoreplicates: number of times to regenerate each fractional
%       data split and recalculate the information improvement. Optional, 
%       default = 1.
%   discrepency_method: string, method to use for the discrepency 
%       measurement. See function calculate_discrepency for options.
%       Default = 'sum_of_squares'.
%   optimize_flag: boolean, when set to true the scaling of the partial 
%       image is optimized to minimize the discrepency. Default = true.
% Output:
%   info_improvement_data: matrix of information improvement fractions,
%       arranged as seperate columns for each pseudoreplicate.

% Set defaults
if nargin < 4; number_pseudoreplicates = 1; end;
if nargin < 5; discrepency_method = 'sum_of_squares'; end;
if nargin < 6; optimize_flag = true; end;

% Get a copy of the STORMvars structure and channel 2 data created with the given parameter structures
testparams = params; % param_array from loading in 'FD_data_SNxx'
testparams.number_events_ch1 = 10;
testparams.number_background_events_ch1 = 10;
[ ~, data_ch2, STORMvars] = create_test_data_dv(testparams, 10);

% Get the full image and maximum discrepency value
full_image = create_test_STORM_images_dv(params, dataset, data_ch2, STORMvars, false, true, true);
max_discrepency = calculate_discrepency(zeros(size(full_image)), full_image, discrepency_method, optimize_flag);

% Initialize results vector
info_improvement_data = zeros(length(fraction_vector), number_pseudoreplicates); 

% Loop through each fraction and calculate the info improvement
for fraction_index = 1:length(fraction_vector)
    fraction = fraction_vector(fraction_index);
    
    % Loop through each pseudoreplicate
    for pseudorep_index = 1:number_pseudoreplicates
            
        % Generate a logical index array for splitting the data
        number_datapoints = size(dataset, 1);
        number_datapoints_taken = round(number_datapoints * fraction);
        bool_selection = false(size(dataset));
        index_integers = randsample(number_datapoints, number_datapoints_taken);                
        bool_selection(index_integers) = true;
        bool_selection(index_integers + number_datapoints) = true; % select all the y points for corrosponding x in linear representation

        % Split the dataset into larger and smaller datasets
        partial_dataset = reshape(dataset(bool_selection), [number_datapoints_taken, 2]);

        % Create images for the larger and smaller datasets
        [partial_image] = create_test_STORM_images_dv(params, partial_dataset, data_ch2, STORMvars, false, true, true);

        % Get the discrepency and information improvement 
        discrepency = calculate_discrepency(partial_image, full_image, discrepency_method, optimize_flag);
        info_improvement_data(fraction_index, pseudorep_index) = 1 - (discrepency / max_discrepency);
    end     
end    
end

