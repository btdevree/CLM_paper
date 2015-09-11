% X-cor test

% 2-24-15 Brian DeVree

% Variables to test
use_new_DFT = true;
use_pdf_flag = false;
pdf_sigma = 4;
pixel_lengths = [50, 40, 30, 25, 20, 16, 12, 10, 8]; % nm
xcor_test_lengths = [5, 10, 15, 20, 25, 30, 35, 40, 50, 60, 80, 100, 150]; % nm
% pixel_lengths = [50, 25, 10]; % Shortened for testing
% xcor_test_lengths = [5, 15, 30, 50];
number_repeats = 10;
number_spots_per_image = 100;
density_of_test_image = .001;
max_xcor_test_length = 250; % nm Handles out to 5 sigma at psf_sigma = 30 nm
analytical_resolution = 2; % nm
analytical_calc_cutoff = 7; % sigma
analytical_radial_average_sampling = 5; %nm

% Initialize the results variables
full_xcor_radii = cell(length(pixel_lengths), length(xcor_test_lengths)); % nm
full_xcor_curves = cell(length(pixel_lengths), length(xcor_test_lengths));
full_xcor_stdev = cell(length(pixel_lengths), length(xcor_test_lengths));
full_xcor_sem = cell(length(pixel_lengths), length(xcor_test_lengths));
expected_xcor_curves = cell(length(pixel_lengths), length(xcor_test_lengths));
residuals = cell(length(pixel_lengths), length(xcor_test_lengths));

% Calculate how big the image needs to be in nm
area_covered_by_spots = number_spots_per_image * pi * pdf_sigma^2; % nm^2
area_of_image = area_covered_by_spots / density_of_test_image; % nm^2
length_of_image = sqrt(area_of_image); % nm

% Loop through the different test lengths
for test_ind = 1:length(xcor_test_lengths)
    xcor_test_length = xcor_test_lengths(test_ind);
         
    % Loop through the different pixel lengths
    for pixel_ind = 1:length(pixel_lengths)
        pixel_length = pixel_lengths(pixel_ind);
    
        % Initialize results matricies (don't know exaclty how many radial data points each will have)
        results_independent = [];
        results_dependent = [];
        
        % Calculate the number of pixels needed per side of image at this resolution
        number_pixels_per_side = ceil(length_of_image / pixel_length);
        
        % Loop through each replicate
        for repeat_ind = 1:number_repeats
                                                
            % Calculate the points used to make the channel 1 image
            new_spots_x = (rand(1, number_spots_per_image) .* number_pixels_per_side); % image generator wants x and y in pixels as row vectors
            new_spots_y = (rand(1, number_spots_per_image) .* number_pixels_per_side);
            data_ch1 = struct('x', new_spots_x, 'y', new_spots_y);
            
            % Calculate the points used to make the channel 2 image
            offest_rho = xcor_test_length / pixel_length; % offest in pixels
            offset_theta = rand(1, number_spots_per_image) .* 2 .* pi;
            [offset_x, offset_y] = pol2cart(offset_theta, offest_rho);
            data_ch2 = struct('x', new_spots_x + offset_x, 'y', new_spots_y + offset_y);
            
            % Make images
            if use_pdf_flag
                covar_matrix = [pdf_sigma.^2, 0; 0, pdf_sigma.^2];
                covar_inv = inv(covar_matrix);
                covar_det = det(covar_matrix);

                % Calc mesh x and y values
                x_vector = [0:1:number_pixels_per_side];
                y_vector = [number_pixels_per_side:-1:0];
                [x_mesh, y_mesh] = meshgrid(x_vector, y_vector);
                for data_ind = 1:length(data_ch1)
                    
                    % Calc pdf 
                    X = [x_mesh(:).';y_mesh(:).'];
                    X_meanshifted = X - repmat([data_ch1(data_ind).x; data_ch1(data_ind).y], 1, length(X)); 
                    pdf = (2 .* pi .* sqrt(covar_det))^-1 .* exp(-0.5 .* sum(X_testshift.' * covar_inv .* X_testshift.', 2));
                    
                    Add to total
                end
                
            elseif use_new_DFT
                % Method from STORM_analyzer_dualview circa October 2014
                % Modifed to use fully zero-padded DFTs for the convolution
                % Note: using original = new pixel length
                res = 1; % desired resolution, from user, in fractions of original pixel length
                s = pdf_sigma / pixel_length; % desired psf width, from user, in fractions of new pixel length
                dims = [number_pixels_per_side, number_pixels_per_side]; % desired image size, from user, in original pixels
                rmax = ceil(max_xcor_test_length/pixel_length); % in new pixels
                [~, image_ch1] = generate_STORM_image_MBS(data_ch1, 1, s, dims);
                [~, image_ch2] = generate_STORM_image_MBS(data_ch2, 1, s, dims);
                                                
                aveval1_raw = sum(sum(image_ch1))/sum(sum(ones(size(image_ch1))));
                factor1_raw = 1/aveval1_raw;
                aveval2_raw = sum(sum(image_ch2))/sum(sum(ones(size(image_ch2))));
                factor2_raw = 1/aveval2_raw;

                column_size = size(image_ch1, 1) + size(image_ch2, 1) - 1;
                row_size = size(image_ch1, 2) + size(image_ch2, 2) - 1;
                N = fftshift(ifft2(abs(fft2(ones(size(image_ch1)), column_size, row_size)).^2));
                C = factor1_raw*factor2_raw*real(fftshift( ...
                    ifft2(fft2(image_ch1, column_size, row_size).* ...
                    conj(fft2(image_ch2, column_size, row_size)))))./N;

                [radii, rad_av_xcor, rad_av_xcor_err] = radial_average_EM(C, rmax); % radii in pixels

                results_independent = pixel_length * radii(2:end); % ignore r = 0
                results_dependent = [results_dependent; rad_av_xcor(2:end)]; % Stack of row vectors with results 
                
            else 
                % Method from STORM_analyzer_dualview circa October 2014
                % Note: using original = new pixel length
                res = 1; % desired resolution, from user, in fractions of original pixel length
                s = pdf_sigma / pixel_length; % desired psf width, from user, in fractions of new pixel length
                dims = [number_pixels_per_side, number_pixels_per_side]; % desired image size, from user, in original pixels
                rmax = ceil(max_xcor_test_length/pixel_length); % in new pixels
                [~, image_ch1] = generate_STORM_image_MBS(data_ch1, 1, s, dims);
                [~, image_ch2] = generate_STORM_image_MBS(data_ch2, 1, s, dims);
                
                aveval1_raw = sum(sum(image_ch1))/sum(sum(ones(size(image_ch1))));
                factor1_raw = 1/aveval1_raw;
                aveval2_raw = sum(sum(image_ch2))/sum(sum(ones(size(image_ch2))));
                factor2_raw = 1/aveval2_raw;


                N = fftshift(ifft2(abs(fft2(ones(size(image_ch1)), size(image_ch1, 1)+rmax-1, size(image_ch1, 2)+rmax-1)).^2));
                C = factor1_raw*factor2_raw*real(fftshift( ...
                    ifft2(fft2(image_ch1, size(image_ch1, 1)+rmax-1, size(image_ch1, 2)+rmax-1).* ...
                    conj(fft2(image_ch2, size(image_ch2, 1)+rmax-1, size(image_ch2, 2)+rmax-1)))))./N;

                [radii, rad_av_xcor, rad_av_xcor_err] = radial_average_EM(C, rmax); % radii in pixels

                results_independent = pixel_length * radii(2:end); % ignore r = 0
                results_dependent = [results_dependent; rad_av_xcor(2:end)]; % Stack of row vectors with results 
            end
        end
        
        % Average results
        av_results_dependent = mean(results_dependent, 1);
        stdev_results_dependent = std(results_dependent, 1);
        
        % Record results
        full_xcor_radii{pixel_ind, test_ind} = results_independent;
        full_xcor_curves{pixel_ind, test_ind} = av_results_dependent;
        full_xcor_stdev{pixel_ind, test_ind} = stdev_results_dependent;
        full_xcor_sem{pixel_ind, test_ind} = stdev_results_dependent./sqrt(number_repeats - 1); % Unbiased estimate
        
        % Calculate analytical solution (work in nanaometers)
        % Calc average spot density
        ch1_total_density = number_spots_per_image / (number_pixels_per_side*pixel_length)^2; % spots/nm^2
        ch2_total_density = number_spots_per_image / (number_pixels_per_side*pixel_length)^2; % spots/nm^2
        ch1_rand_density = (number_spots_per_image - 1) / (number_pixels_per_side*pixel_length)^2; % spots/nm^2 Density of all other spots in the channel besides the one modeled in the psf
        ch2_rand_density = (number_spots_per_image - 1) / (number_pixels_per_side*pixel_length)^2; % spots/nm^2 Density of all other spots in the channel besides the one modeled in the psf
        covar_matrix = [pdf_sigma.^2, 0; 0, pdf_sigma.^2];
        covar_inv = inv(covar_matrix);
        covar_det = det(covar_matrix);
        
        % Calc mesh x and y values
        x_vector = [-max_xcor_test_length:analytical_resolution:max_xcor_test_length];
        y_vector = [-max_xcor_test_length:analytical_resolution:max_xcor_test_length];
        [x_mesh, y_mesh] = meshgrid(x_vector, y_vector);
        
        % Calc psf of ch2
        X = [x_mesh(:).';y_mesh(:).'];
        X_testshift = X - repmat([xcor_test_length; 0], 1, length(X)); 
        ch2_pdf = (2 .* pi .* sqrt(covar_det))^-1 .* exp(-0.5 .* sum(X_testshift.' * covar_inv .* X_testshift.', 2));
        
        
        % Radially average to get the expected c(r) curve
        % Loop through the desired indices
        expected_radav_xcor = zeros(size(results_independent));
        for rad_ind = 1:length(results_independent) 
            radius = results_independent(rad_ind);
            
            % Spread out sampling grid evenly along circumference 
            % from 0 to pi, only calculate one half of the circle because the 
            % other half is symmetric. Don't calculate when centers of the
            % pdf are farther away than the analytical_calc_cultoff
            total_circ_length = 2 .* pi .* radius; % circumfrence in nm
            half_num_points = floor((total_circ_length./2)./analytical_radial_average_sampling) + 1;
            total_num_points = (2 .* half_num_points) - 2;
            actual_interval = total_circ_length ./ total_num_points;

            % Decide which points to average.
            % Compute law of cosines expression a^2 + b^2 - c^2 / 2ab
            cosine_law_product = (radius.^2 + xcor_test_length.^2 -...
                (analytical_calc_cutoff .* pdf_sigma).^2) ./ (2 .* radius .* xcor_test_length);
            
            % If the product is <-1, the two pdf centers are always within
            % the cutoff length and we need to calculate the entire half circle
            if cosine_law_product <= -1
                half_circ_lengths = linspace(0, total_circ_length ./ 2, half_num_points);
                repeat_last_circ_point = false;
                num_uneval_points = 0;
               
            % If the product is between -1 and 1, there exists an angle where the two centers
            % are past the cutoff and we need to go from 0 radians to the closest point to that angle
            elseif -1 < cosine_law_product && cosine_law_product < 1
                max_angle = acos(cosine_law_product);
                max_circumference = max_angle .* radius;
                num_half_circ_points = floor(max_circumference ./ actual_interval);
                half_circ_lengths = linspace(0, num_half_circ_points .* actual_interval, num_half_circ_points);
                repeat_last_circ_point = true;
                num_uneval_points = total_num_points - num_half_circ_points .* 2 - 1;
            
            % If the product is >1, there is no angle where the two centers
            % are less than the cutoff and we don't need to evaluate anything
            elseif 1 < cosine_law_product
                half_circ_lengths = [];
                num_uneval_points = total_num_points;
                repeat_last_circ_point = false; % Doesn't matter for an empty matrix
            end
            
            % Calculate radians and convert to x,y coordinates
            radians = half_circ_lengths./radius;
            [x_circ, y_circ] = pol2cart(radians, radius);
                                  
            % Calculate psf of ch1 at the specified radial positions
            xcor_circ = zeros(length(x_circ), 1);
            for circ_ind = 1:length(x_circ)
                xshift = x_circ(circ_ind);
                yshift = y_circ(circ_ind);
                shift = [xshift; yshift];
                X_radshift = X - repmat(shift, 1, length(X)); 
                ch1_pdf = (2 .* pi .* sqrt(covar_det))^-1 .* exp(-0.5 .* sum(X_radshift.' * covar_inv .* X_radshift.', 2));
            
                % Calculate the xcor at the specified radial positions
                xcor = sum(ch1_pdf .* ch2_pdf) ./ sum(ch1_pdf .* ch2_total_density) +...
                    (ch1_rand_density .* ch2_rand_density) ./ (ch1_total_density .* ch2_total_density);
               
                % Add to results vector
                xcor_circ(circ_ind) = xcor;
            end
            
            % Duplicate the 2nd half of the curve to save computation
            if repeat_last_circ_point
                otherhalf_xcor = xcor_circ(2:end); % Don't repeat 0 radians
            elseif ~repeat_last_circ_point
                otherhalf_xcor = xcor_circ(2:end-1); % Don't repeat 0 and pi radians
            end
            
            % Add both halves together
            xcor_circ = [xcor_circ; otherhalf_xcor];
            
            % Average the values
            eval_sum = sum(xcor_circ); % Gives zero if xcor_circ is empty
            uneval_sum = num_uneval_points .* (ch1_rand_density .* ch2_rand_density) ./ (ch1_total_density .* ch2_total_density);
            expected_radav_xcor(rad_ind) = (eval_sum + uneval_sum) ./ (length(xcor_circ) + num_uneval_points);                
        end

        % Calc residuals
        xcor_residuals = av_results_dependent - expected_radav_xcor;
        
        % Record the values
        expected_xcor_curves{pixel_ind, test_ind} = expected_radav_xcor;
        residuals{pixel_ind, test_ind} = xcor_residuals;
        fprintf('.')
    end
    fprintf('\n')
end

% Save data
save('data.mat', 'full_xcor_radii', 'full_xcor_curves', 'full_xcor_stdev', 'full_xcor_sem',...
    'expected_xcor_curves', 'residuals', 'pixel_lengths', 'xcor_test_lengths');
