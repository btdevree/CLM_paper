function [timing_vector, whos_struct_cell_vector] = calc_Fourier_timing(number_points_vector, method)
%CALC_FOURIER_TIMING Calculate the run time and memory use for a Fourier 
%   calculation method of autocorrelations 
%
%   Calculates a radial pairwise autocorrelation curve from points
%   uniformaly distributed in a unit square and records the time of
%   execution and the memory use for the calculation.
%   Input:
%   number_points_vector: Column vector of doubles; the number of points to 
%       use in the autocorrelation
%   method: string denoting the method used for calculation. Choices are
%       'binning' and 'Gaussian_pdf'.
%   Output:
%   timing_vector: Time of execution for calculating the correlation and
%       the radial average. Best of 10 repetitions for number_of_points <= 
%       1e4, best of 3 repetitions for number_of_points <= 1e6, and a 
%       single repetition for all else. Column vector of doubles.
%   whos_struct_cell_vector: output of the whos command at the end of each
%       requested number of points. Column vector of cells containing 
%       strucures.

% Define repetition limits
repeat_limits = [5e3, 5e5];

% Define pixel resolution and max length
full_image_size = 5e4; % nanometers
STORM_pixel_resolution = 7; % nanometers
max_length = 1000; % nanometersm

% Loop through each number_of_points value
for number_points_index = 1:size(number_points_vector, 1);
    number_points = number_points_vector(number_points_index);

    % Create list of points
    points = rand(number_points, 2);

    % Determine the number of repeats
    if number_points <= repeat_limits(1)
        number_repeats = 1;
    elseif number_points <= repeat_limits(2)
        number_repeats = 1;
    else
        number_repeats = 1;
    end

    % Repeat the specified number of times
    for repeat_index = 1:number_repeats

        % Start timer
        tic

        % Create STORM image with  binning coordinates, a la Jennifer Lippincott-Schwartz
        if strcmp(method, 'binning')

            % Create STORM image
            image_size = ceil(full_image_size / STORM_pixel_resolution);
            STORM_image = zeros(image_size);

            % Populate STORM image with coordinates
            for coord_index = 1:size(points, 1);
                
                % Calculate the index values for the given coordinates
                column_index = floor(points(coord_index, 1) * image_size); % x coord
                row_index = ceil(image_size - points(coord_index, 2) * image_size); % y coord
                
                % Add to the STORM image
                STORM_image(row_index, column_index) = STORM_image(row_index, column_index) + 1;
            end
            
            % Calculate the autocorrelation of the image
            autocorrelation = calc_crosscorrelation(STORM_image, STORM_image, ceil(max_length/STORM_pixel_resolution)); 
                
            % Calculate radial average with binning
            [distance_vector, mean_vector, stdev_vector, sem_vector] = radial_average_2D_correlation_binning(autocorrelation);
                
        % Create STORM image with a sampled Gaussian psf
        elseif strcmp(method, 'Gaussian_psf')
            
            % Convert points to a data structure
            data_struct = struct();
            data_struct.x = points(:, 1);
            data_struct.y = points(:, 2);
            
            % Call the STORM image creating function; pretend an origonal pixel size = 70 nm
            original_pixel_size = 70;
            resolution = STORM_pixel_resolution/original_pixel_size; % fraction of original pixel, 7nm
            sigma = 25/original_pixel_size; % fraction of original pixel, 25nm
            dims = repmat(ceil(full_image_size / original_pixel_size), [1, 2]); % in original pixels
            STORM_image = create_STORM_image(data_struct, resolution, sigma, dims); % note that output is sparse
            
            % Convert to dense image
            STORM_image = full(STORM_image);

            % Calculate the autocorrelation of the image
            autocorrelation = calc_crosscorrelation(STORM_image, STORM_image, ceil(max_length/STORM_pixel_resolution));

            % Calculate the radial average of the autocorrelation
            [distance_vector, mean_vector, stdev_vector, sem_vector] = radial_average_2D_correlation(autocorrelation);
        end
        
        % We're done, get the timing
        timing_repeats(repeat_index) = toc;
    end
    
    % Take the best timing
    timing_vector(number_points_index) = min(timing_repeats);
    
    % Record the whos vector
    whos_struct_cell_vector{number_points_index} = whos;
end

% Convert outputs to column vectors
timing_vector = timing_vector.';
whos_struct_cell_vector = whos_struct_cell_vector.';
end


