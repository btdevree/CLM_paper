% Script to generate aliaised cross-correlations 

% Set directory information
path_parts = strsplit(pwd, 'CLM_paper');
filepath = [path_parts{1}, 'CLM_figures_and_data/'];
data_name = 'aliasing_part_B.mat';

% Load a parameter structure from current dirctory with deisred settings
load('aliasing_params.mat') % loads as 'params'

% Set parameters
number_replicates = 30;
max_correlation_radius = 500; %nm
pixel_lengths = [60; 55; 50; 45; 40; 36; 32; 28; 25; 22; 19; 16; 14; 12; 10; 8; 7]; % nm
analytical_calc_cutoff = 10; % sigma
radial_average_sampling = 5; %nm

% Initialize the results variables
xcor_curves =  cell(length(pixel_lengths), number_replicates);
residuals = cell(length(pixel_lengths), number_replicates);
analytical_curves = cell(length(pixel_lengths), 1);
dist_vectors = cell(length(pixel_lengths), 1);

% Loop through the different pixel lengths
for pixel_ind = 1:length(pixel_lengths)
    pixel_length = pixel_lengths(pixel_ind);
    
    % Set pixel size
    params.STORM_pixel_size = pixel_length;
    
    % Loop through each replicate
    for repeat_ind = 1:number_replicates
        
        % Get event data for the image
        [ch1_data, ch2_data, passed_vars ] = create_test_data_dv(params);
        
        % Create a STORM image
        [ch1_STORM_image, ch2_STORM_image] = create_test_STORM_images_dv(params, ch1_data, ch2_data, passed_vars, false, true, true);
        
        % Create an image mask
        passed_vars.cell_radius = passed_vars.cell_radius + 2*params.ch2_crosscor_params; % Add 2x xcor radius to mask 
        STORM_mask = create_test_cell_STORM_mask(params, passed_vars);
        
        % Run the cross-correlation
        max_radius_px = ceil(max_correlation_radius/params.STORM_pixel_size);
        xcor_image = calc_crosscorrelation(ch1_STORM_image, ch2_STORM_image, max_radius_px, STORM_mask);
        
        % Radially average the xcor image
        radial_sampling_px = radial_average_sampling/pixel_length;
        [distance_vector_px, mean_vector] = radial_average_2D_correlation(xcor_image, radial_sampling_px);
        
        % Calculate the expected result and record the distance and analytical results, first cycle only
        if repeat_ind == 1
           
            % Record the distance vector in nanometers
            distance_vector = pixel_length * distance_vector_px;
            dist_vectors{pixel_ind, 1} = distance_vector;
            
            % Add required info to the passed_vars structure
            passed_vars.analytical_arc_length_sampling = radial_average_sampling/3;
            passed_vars.analytical_calc_cutoff = analytical_calc_cutoff; % sigma
            
            % Get analytical solution
            expected_radav_xcor = alaising_test_analytical_solution(distance_vector, params, passed_vars);
        end
        
        % Record the results
        xcor_curves{pixel_ind, repeat_ind} = mean_vector;
        
        % Record residuals
        residuals{pixel_ind, repeat_ind} = mean_vector - expected_radav_xcor;
        
        % Report finishing a curve
        fprintf('.')
    end
    % Report finishing a set of curves
    fprintf('\n')
end

% Save data
data_filename = [filepath, data_name];
save(data_filename, 'xcor_curves', 'residuals', 'analytical_curves', 'dist_vectors', 'pixel_lengths');
