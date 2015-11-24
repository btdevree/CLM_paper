% Make part C figure data

% Calculates the completeness index data for the part C figure. Assumes we
% can find the images and other data at ..CLM_figures_and_data under the 
% names 'NPIF_part_C_data.mat' and 'NPIF_part_C_images.h5'.

% We assume that we're in the CLM_paper repository, and we want to save the
% big binary figure and data files to CLM_figures_and_data folder, not on 
% the repository but on the same file level as CLM_paper.
%binary_path_parts = strsplit(pwd, 'CLM_paper');
%binary_path = [binary_path_parts{1}, 'CLM_figures_and_data/'];
binary_path = '/home/btdevree/large_file_temp/'; % Network drive is just too slow and causes process to get killed

% Re-define SN-ratios
SN_ratios = [1; 4; 10];

% These files are pretty big, so we'll have to make each SNratio seperately
for SN_index = 1:length(SN_ratios)
    SN_ratio = SN_ratios(SN_index);
   
    % Initialize hdf5 array to hold images directly on disc
    image_filepath = [binary_path, 'NPIF_part_C_images_SN', num2str(SN_ratio), '.h5'];
    data_filepath = [binary_path, 'NPIF_part_C_data_SN', num2str(SN_ratio), '.mat']; 
   
    % Read in the data file
    load(data_filepath);
    % Re-define SN-ratios
    SN_ratios = [1; 4; 10];

    % Define number of pseudoreplicates for the approximation completeness index, get number of replicates
    approx_pseudoreplicates = 5;

    % Get a copy of the STORMvars structure and channel 2 data created with the given parameter structures
    testparams = param_array{1, 1, 1};
    testparams.number_events_ch1 = 10;
    testparams.number_background_events_ch1 = 10;
    [ ~, data_ch2, ~, STORMvars] = create_test_data_dv(testparams, 10);

    % Generate seeds 
    seeds = round(rand(size(sliced_param_array)).*1e6); % Can specify to reproduce the same images

    % Get image size
    image_info = h5info(image_filepath, '/images_array');
    image_dims = image_info.Dataspace.Size; % Note that Matlab translates this to row-major standard
    image_height = image_dims(1);
    image_width = image_dims(2);

    % Initialize a matrix for the completeness index data
    ideal_CI_method1 = zeros(size(sliced_param_array));
    ideal_CI_method2 = zeros(size(sliced_param_array));
    approx_CI_method1 = zeros([size(sliced_param_array), approx_pseudoreplicates]);
    approx_CI_method2 = zeros([size(sliced_param_array), approx_pseudoreplicates]);

    % Loop through each parameter struture and corrosponding image
    for param_index = 1:length(sliced_param_array(:));
        params = sliced_param_array{param_index};

        % Give user feedback
        disp(['Working on, S/N ratio = ', num2str(params.SN_ratio), ', event number = ',...
            num2str(params.number_events_ch1), ', replicate = ', num2str(params.replicate)]);
        
        % Channge the random number generator seed
        rng(seeds(param_index));

        % Load in the image
        image = h5read(image_filepath, '/images_array', [1, 1, 1, params.output_indices(2:3)], [image_height, image_width, 1, 1, 1]);

        % Calculate the ideal completeness index
        ideal_discrepency_method1 = calculate_discrepency(image, ideal_image, 'sum_of_squares');
        ideal_discrepency_method2 = calculate_discrepency(image, ideal_image, 'l2_norm');
        ideal_CI_method1(param_index) = 1-ideal_discrepency_method1;
        ideal_CI_method2(param_index) = 1-ideal_discrepency_method2;

        % Get the full dataset
        data_100pct = datasets{param_index};

        % Calculate the approximated completeness index
        % Repeat multiple times for pseudo-replicates 
        for pseudorep_index = 1:approx_pseudoreplicates

            % Generate a logical index array for splitting the data
            number_datapoints_100pct = size(data_100pct, 1);
            number_datapoints_1pct = ceil(number_datapoints_100pct/100);
            number_datapoints_99pct = number_datapoints_100pct - number_datapoints_1pct;
            bool_selection = false(size(data_100pct));
            integers_1pct = randsample(number_datapoints_100pct, number_datapoints_1pct);
            for int_index = 1:length(integers_1pct);
                integer_value = integers_1pct(int_index);
                bool_selection(integer_value, :) = true;
            end

            % Split the dataset into 1% and 99% datasets
            data_1pct = reshape(data_100pct(bool_selection), [number_datapoints_1pct, 2]);
            data_99pct = reshape(data_100pct(~bool_selection), [number_datapoints_99pct, 2]);

            % Create images for the 1% and 99% datasets
            [image_1pct, ~, ~] = create_test_STORM_images_dv(params, data_1pct, data_ch2, STORMvars, false);
            [image_99pct, ~, ~] = create_test_STORM_images_dv(params, data_99pct, data_ch2, STORMvars, false);
            
            % Get the discrepency in a blank image
            zero_image = zeros(size(image));
            discrepency_0pct_method1 = calculate_discrepency(zero_image, image_99pct, 'sum_of_squares');
            discrepency_0pct_method2 = calculate_discrepency(zero_image, image_99pct, 'l2_norm');
            
            % Get the discrepency in the 1% image
            discrepency_1pct_method1 = calculate_discrepency(image_1pct, image_99pct, 'sum_of_squares');
            discrepency_1pct_method2 = calculate_discrepency(image_1pct, image_99pct, 'l2_norm');

            % Get the discrepency increase by adding the 1% image to the 99% image
            discrepency_100pct_method1 = calculate_discrepency(image, image_99pct, 'sum_of_squares');
            discrepency_100pct_method2 = calculate_discrepency(image, image_99pct, 'l2_norm');

            % Calculate completeness approximation
            CI_method1 = discrepency_100pct_method1 / (discrepency_0pct_method1 - discrepency_1pct_method1);
            CI_method2 = discrepency_100pct_method2 / (discrepency_0pct_method2 - discrepency_1pct_method2);
            
            out_ind = params.output_indices;
            approx_CI_method1(1, out_ind(2), out_ind(3), pseudorep_index) = CI_method1;
            approx_CI_method2(1, out_ind(2), out_ind(3), pseudorep_index) = CI_method2;
        end
    end
    
    % Save CI info
    CI_filename = ['NPIF_part_C_CIdata_ssq99_SN', num2str(SN_ratio), '.mat'];
    save(CI_filename, 'ideal_CI_method1', 'ideal_CI_method2', 'approx_CI_method1', 'approx_CI_method2',...
        'param_array', 'SN_ratios', 'true_event_numbers', 'seeds');
end