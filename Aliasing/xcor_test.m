% Script to generate aliaised cross-correlations 

% Set directory information
path_parts = strsplit(pwd, 'CLM_paper');
filepath = [path_parts{1}, 'CLM_figures_and_data/'];
data_name = 'aliasing_part_B_pdf_200nm.mat';

% Load a parameter structure from current dirctory with deisred settings
load('aliasing_params.mat') % loads as 'params'

% Set parameters
number_replicates = 30;
max_correlation_radius = 1000; %nm
pixel_lengths = logspace(log10(75), log10(5), 30); % nm
analytical_calc_cutoff = 10; % sigma
radial_average_sampling = 5; %nm
use_bins = false;
test_radius = 200;

if use_bins
    params.ch2_crosscor = 'Gaussian';
    params.ch2_crosscor_params = [test_radius, sqrt(2 * params.STORM_precision.^2)];
else
    params.ch2_crosscor = 'exact';
    params.ch2_crosscor_params = test_radius;
end
    
% Initialize the results variables
xcor_curves =  cell(length(pixel_lengths), number_replicates);
residuals = cell(length(pixel_lengths), number_replicates);
analytical_curves = cell(length(pixel_lengths), 1);
dist_vectors = cell(length(pixel_lengths), 1);

% Loop through each replicate
for repeat_ind = 1:number_replicates
    
    % Get event data for the image
    [ch1_data, ch2_data, passed_vars] = create_test_data_dv(params);
        
    % Loop through the different pixel lengths
    for pixel_ind = 1:length(pixel_lengths)
    pixel_length = pixel_lengths(pixel_ind);
    
        % Set pixel size
        params.STORM_pixel_size = pixel_length;
       
        % Create a STORM image
        if use_bins
            [ch1_STORM_image, ch2_STORM_image] = create_test_STORM_images_dv(params, ch1_data, ch2_data, passed_vars, false, true, true, true);
        else
            [ch1_STORM_image, ch2_STORM_image] = create_test_STORM_images_dv(params, ch1_data, ch2_data, passed_vars, false, true, true);
        end
            
        % Create an image mask
        modified_passed_vars = passed_vars;
        modified_passed_vars.cell_radius = passed_vars.cell_radius - 2*test_radius; % subtract 2x xcor radius to mask 
        STORM_mask = create_test_cell_STORM_mask(params, modified_passed_vars);
        
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
            expected_radav_xcor = aliasing_test_analytical_solution(distance_vector, params, passed_vars);
            analytical_curves{pixel_ind} = expected_radav_xcor;
        end
        
        % Record the results
        xcor_curves{pixel_ind, repeat_ind} = mean_vector;
        
        % Record residuals
        residuals{pixel_ind, repeat_ind} = mean_vector - analytical_curves{pixel_ind};
        
        % Report finishing a curve
        fprintf('.')
    end
    % Report finishing a set of curves
    fprintf('\n')
end

% Save data
data_filename = [filepath, data_name];
save(data_filename, 'xcor_curves', 'residuals', 'analytical_curves', 'dist_vectors', 'pixel_lengths');

