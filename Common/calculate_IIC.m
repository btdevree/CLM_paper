function [info_improvement_data, exp_discrepency_data, ideal_discrepency_data] = calculate_IIC(params, dataset, fraction_vector,...
    number_pseudoreplicates, discrepency_method, ideal_image, optimize_flag)
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
%   discrepency_method: string or cell array of strings, method to use for 
%       the discrepency measurement. See function calculate_discrepency for 
%       options. Default = 'sum_of_squares'.
%   ideal_image: ideal image to measure the ideal, theoretical discrepency.
%       Requried to output the ideal discrepency matrix. Default = [] (no 
%       ideal calculations).
%   optimize_flag: boolean, when set to true the scaling of the partial 
%       image is optimized to minimize the discrepency. Default = true.
% Output:
%   info_improvement_data: matrix or cell array of matricies of information
%       improvement fractions, arranged as seperate columns for each 
%       pseudoreplicate.
%   experimental_discrepency: raw experimental image discrepency values,
%       arranged in the same order as the info_improvement_data.
%   ideal_discrepency: raw ideal image discrepency values, arranged in the 
%       same order as the IIC_results. Must supply ideal image in order to
%       calculate.



% Set defaults
if nargin < 4; number_pseudoreplicates = 1; end;
if nargin < 5; discrepency_method = 'sum_of_squares'; end;
if nargin < 6; ideal_image = []; end;
if nargin < 7; optimize_flag = true; end;

% If the discrepency_method is a simple string, put it into a cell array
if ischar(discrepency_method)
    discrepency_method = {discrepency_method};
end

% Get a copy of the STORMvars structure and channel 2 data created with the given parameter structures
testparams = params; % param_array from loading in 'FD_data_SNxx'
testparams.number_events_ch1 = 10;
testparams.number_background_events_ch1 = 10;
[ ~, data_ch2, STORMvars] = create_test_data_dv(testparams, 10);

% Get the full image
full_image = create_test_STORM_images_dv(params, dataset, data_ch2, STORMvars, false, true, true);

% Loop through each discrepency method
exp_discrepency_data = cell(0);
ideal_discrepency_data = cell(0);
max_discrepency = cell(0);
for result_cell_index = 1:length(discrepency_method)

    % Initialize results matricies
    exp_discrepency_data{result_cell_index} = zeros(length(fraction_vector), number_pseudoreplicates);
    if ~isempty(ideal_image)
        ideal_discrepency_data{result_cell_index} = zeros(length(fraction_vector), number_pseudoreplicates);
    end
    
    % Get the maximum discrepency for each method
    max_discrepency{result_cell_index} = calculate_discrepency(zeros(size(full_image)), full_image, discrepency_method{result_cell_index}, optimize_flag);
end

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

        % Reduce the dataset into a partial dataset
        partial_dataset = reshape(dataset(bool_selection), [number_datapoints_taken, 2]);

        % Create images for the larger and smaller datasets
        [partial_image] = create_test_STORM_images_dv(params, partial_dataset, data_ch2, STORMvars, false, true, true);
         
        % Loop through each discrepency method
        for result_cell_index = 1:length(discrepency_method)
            
            % Record the discrepency 
            exp_discrepency_data{result_cell_index}(fraction_index, pseudorep_index) = calculate_discrepency(partial_image, full_image,...
                discrepency_method{result_cell_index}, optimize_flag);
            if ~isempty(ideal_image)
                ideal_discrepency_data{result_cell_index}(fraction_index, pseudorep_index) = calculate_discrepency(partial_image,...
                    ideal_image, discrepency_method{result_cell_index}, optimize_flag);
            end
        end
    end
end

% Loop through each discrepency method
info_improvement_data = cell(0);
for result_cell_index = 1:length(discrepency_method)

    % Transform to information improvement
    info_improvement_data{result_cell_index} = 1 - (exp_discrepency_data{result_cell_index} / max_discrepency{result_cell_index});
end

% If only one discrepency method is requested, return just the results and not the cell array
if length(discrepency_method) == 1
    info_improvement_data = info_improvement_data{1};
    exp_discrepency_data = exp_discrepency_data{1};
    ideal_discrepency_data = ideal_discrepency_data{1};
end
end

