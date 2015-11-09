% Make part C figure data

% We assume that we're in the CLM_paper repository, and we want to save the
% big binary figure and data files to CLM_figures_and_data folder, not on 
% the repository but on the same file level as CLM_paper.
binary_path_parts = strsplit(pwd, 'CLM_paper');
binary_path = [binary_path_parts{1}, 'CLM_figures_and_data/'];

% Define variable quantities for the S/N and event_number
SN_ratios = [1; 4; 10];
true_event_numbers = [1e3; 2e4];% 5e5; 1e4; 2e4; 5e4; 1e5; 2e5; 5e5; 1e6; 2e6; 5e6; 1e7];
replicates = 2;
num_ratios = size(SN_ratios, 1);
num_eventnums = size(true_event_numbers, 1);

% Generate seeds 
seeds = round(rand(num_ratios, num_eventnums, replicates).*1e6); % Can specify to reproduce the same images

% Initialize cell arrays to hold the datasets and parameters arrays
datasets = cell(num_ratios, num_eventnums, replicates);
param_array = cell(num_ratios, num_eventnums, replicates);

% Load in the image creation parameters as variable params
load 'params_part_B'

% Copy the parameter structure into a cell array so that we can make the images in parallel
% Loop through all SN_ratio and event_number options
for SN_index = 1:size(SN_ratios, 1)
    SN_ratio = SN_ratios(SN_index);
    for event_num_index = 1:size(true_event_numbers, 1)
        event_number = true_event_numbers(event_num_index);
        
        % Edit the parameter structure and make a copy in the parameter array
        params.number_events_ch1 = event_number;
        params.number_background_events_ch1 = ceil(event_number/SN_ratio);
        params.SN_ratio = SN_ratio;
        params.binary_path = binary_path;
        for replicate_index = 1:replicates
            params.replicate = replicate_index;
            params.output_indices = [SN_index, event_num_index, replicate_index];
            param_array{SN_index, event_num_index, replicate_index} = params;
        end
    end
end

% Make a test image to see what size the results will be
test_params = param_array{1};
test_params.number_events_ch1 = 10;
test_params.number_background_events_ch1 = 10;
[test_image, ~, ~] = part_C_create_image_and_data(test_params, 10);
image_height = size(test_image, 1);
image_width = size(test_image, 2);

% Initialize array to hold images directly on disc
image_file = matfile([binary_path, 'NPIF_part_C_images.mat'], 'Writable', true);
image_file.images_array = zeros(image_height, image_width, num_ratios, num_eventnums, replicates);

% Calculate the images in parallel; run and save each SN ration seperate
number_results = length(param_array(:));
for index = 1:number_results

    % Evaluate the function asynchronously
    future_results(index) = parfeval(@part_C_create_image_and_data, 3, param_array{index}, seeds(index)); 
end

% Collect and save results
for index = 1:number_results

  % Get next available result, fetchNext blocks until next results are available.
  [completedIdx, image, data, params] = fetchNext(future_results);
  
  % Must use n-dimensional indices to save image directly to disc
  output_indices = params.output_indices;
  image_file.images_array(:, :, output_indices(1), output_indices(2), output_indices(3)) = image;
  
  % save data file
  datasets{output_indices(1), output_indices(2), output_indices(3)} = data;
  
  % Let user know about progress
    disp(['finished S/N ratio = ', num2str(params.SN_ratio), ', event number = ',...
        num2str(params.number_events_ch1), ', replicate = ', num2str(params.replicate)]);

end  
    
% Get the ideal image
ideal_image = calculate_ideal_image(param_array{1}); % Assume same relevent parameters for the ideal image of all images

% Save the datafiles because they can take a long time to generate
save([binary_path, 'NPIF_part_C_data.mat'], 'binary_path', 'param_array', 'SN_ratios', 'true_event_numbers', 'ideal_image', 'datasets', '-v7.3');

