% Make part C figure data

% Calculates the completeness index data for the part C figure. Assumes we
% can find the images and other data at ..CLM_figures_and_data under the 
% names 'NPIF_part_C_data' and 'NPIF_part_C_images_SN[SN ratio value]'.

% We assume that we're in the CLM_paper repository, and we want to save the
% big binary figure and data files to CLM_figures_and_data folder, not on 
% the repository but on the same file level as CLM_paper.
binary_path_parts = strsplit(pwd, 'CLM_paper');
binary_path = [binary_path_parts{1}, 'CLM_figures_and_data/'];

% Read in the data file
load([binary_path, 'NPIF_part_C_data.mat']);

% Define number of pseudoreplicates for the approximation completeness index, get number of replicates
approx_pseudoreplicates = 5;
replicates = size(param_array, 3);

% Initialize a matrix for the completeness index data
ideal_completeness_index = zeros(SN_ratios, true_event_numbers, size(param_array, 3));
approx_completeness_index = zeros(SN_ratios, true_event_numbers, size(param_array, 3), approx_pseudoreplicates);

% Loop through each group of images
for SN_index = 1:size(SN_ratios, 1);
    SN_ratio = SN_ratios(SN_index);
    
    % Load in the group of images
    load([binary_path, 'NPIF_part_C_images_SN', num2str(SN_ratio), '.mat'], 'images_array');
    
    % Loop through each number of events and replicates
    for event_num_index = 1:size(true_event_numbers, 1)
        for replicate_index = 1:replicates
            
            % Calculate the ideal completeness index
            ideal_discrepency = calc_discrepency(images_array{1, event_num_index, replicate_index},...
                ideal_image, 'normalized_variation_of_information', 100, false);
            ideal_completeness_index(SN_index, event_num_index, replicate_index) = 1-ideal_discrepency;
            
            % Calculate the 
            
        end
    end
    
    % Calc surprisal
    surprisal_values(index) = calculate_surprisal(images_array{index}, ideal_image, 'squared_distance');    
end

% Take the ratio of the surprisal to the zero image surprisal
surprisal_ratios = surprisal_values ./ zero_image_surprisal;

real_completeness_index = zeros(SN_ratios, true_event_numbers, size(param_array, 3), );