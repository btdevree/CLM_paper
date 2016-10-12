function [IIC_results, experimental_discrepency, ideal_discrepency] = generate_IIC_curve(params, fraction_vector, number_events_vector, SNratio,...
    number_replicates, number_pseudoreplicates, discrepency_method, ideal_image, optimize_flag, verbose_flag)
%GENERATE_IIC_CURVE Generate information improvement characteristic curves 
% for test images created with the create_test_data_dv data generator.

% Input:
%   params: parameters structure for image generation
%   fraction_vector: column vector of the fractions of data used, include 
%       full and zero fractions.
%   number_events_vector: column vector of the number of events used to
%       make each full image.
%   SNratio: signal to noise ratio, number of true points / number of 
%       randomly distributed points. 
%   number_replicates: number of times to regenerate each dataset.
%       Optional, default = 1;
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
%   verbose_flag: boolean, when set to true the function outputs progress
%       to the console. Default = false;
% Outputs:
%   IIC_results: information improvement characteristic data, arranged with
%       the fraction of events used on the first dimention, number of 
%       events on the second dimension, and replicates on the third 
%       dimension.
%   experimental_discrepency: raw experimental image discrepency values,
%       arranged in the same order as the IIC_results.
%   ideal_discrepency: raw ideal image discrepency values, arranged in the 
%       same order as the IIC_results. Must supply ideal image in order to
%       calculate.

% Set defaults
if nargin < 5; number_replicates = 1; end;
if nargin < 6; number_pseudoreplicates = 1; end;
if nargin < 7; discrepency_method = 'sum_of_squares'; end;
if nargin < 8; ideal_image = []; end;
if nargin < 9; optimize_flag = true; end;
if nargin < 10; verbose_flag = false; end;

% If the discrepency_method is a simple string, put it into a cell array
if ischar(discrepency_method)
    discrepency_method = {discrepency_method};
end

% Initialize results matricies
IIC_results = cell(0);
experimental_discrepency = cell(0);
ideal_discrepency = cell(0);
for result_cell_index = 1:length(discrepency_method)
    IIC_results{result_cell_index} = zeros(length(fraction_vector), length(number_events_vector), number_replicates);
    experimental_discrepency{result_cell_index} = zeros(length(fraction_vector), length(number_events_vector), number_replicates);
    if ~isempty(ideal_image)
        ideal_discrepency{result_cell_index} = zeros(length(fraction_vector), length(number_events_vector), number_replicates);
    else
        ideal_discrepency{result_cell_index} = [];
    end
end

% start timer
elapsed_time = 0;
tic

% Loop through each set of parameters to make an IIC curve
for num_events_index = 1:length(number_events_vector)
    number_events = number_events_vector(num_events_index);
    for replicate_index = 1:number_replicates  

        % Report
        if verbose_flag
            if replicate_index == 1
                fprintf('\n Working on event number %1.1E replicate  1 t = %3d min', number_events, round(elapsed_time/60));
            else
                elapsed_time = toc;
                fprintf('\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b %2d t = %3d min', replicate_index, round(elapsed_time/60));
            end
        end

        % Edit parameters
        true_events = round(number_events / (1 + 1 / SNratio));
        params.number_events_ch1 = true_events;
        params.number_background_events_ch1 = number_events - true_events;

        % Generate new data
        seed = randi(1e6);
        [dataset] = create_test_data_dv(params, seed);

        % Get IIC curve
        [IIC_pesudoreps, exp_pesudoreps, ideal_pesudoreps] = calculate_IIC(params, dataset, fraction_vector,...
            number_pseudoreplicates, discrepency_method, ideal_image, optimize_flag, true);
         
        % Take the mean values of the pseudoreplicates and put into result matrices
        for result_cell_index = 1:length(discrepency_method)
            IIC_results{result_cell_index}(:, num_events_index, replicate_index) = mean(IIC_pesudoreps{result_cell_index}, 2);
            experimental_discrepency{result_cell_index}(:, num_events_index, replicate_index) = mean(exp_pesudoreps{result_cell_index}, 2);
            if ~isempty(ideal_image)
                ideal_discrepency{result_cell_index}(:, num_events_index, replicate_index) = mean(ideal_pesudoreps{result_cell_index}, 2);
            end
        end
    end
end

% If only one discrepency method, return just the results and not the cell array
if length(discrepency_method) == 1
    IIC_results = IIC_results{1};
    experimental_discrepency = experimental_discrepency{1};
    ideal_discrepency = ideal_discrepency{1};
end

% Report
if verbose_flag; fprintf('\n'); end;
end

