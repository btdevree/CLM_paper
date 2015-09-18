function sparse_image = create_STORM_image(data, resolution, sigma, dims, parallel_flag, calc_cutoff_sigmas, zero_cutoff)
    % CREATE_STORM_IMAGE Creates a high-resolution image built from 
    %   Gaussian pdfs. 
    %
    % Inputs:
    % data: Structure with fields 'x', 'y', and, optionally, 'sigma' for 
    %   each data point. Assumes a coordinate system with the orgin at the 
    %   top left corner of the top left pixel, and values increase with 
    %   value as one goes downward or rightward in the image. Given in 
    %   original pixel values.
    % resolution: Desired resolution in fractions of original pixel length.
    % sigma: Desired psf width in fractions of original pixel length. Can
    % be a scalar or 'per-point'.
    % dims: Desired image size in original pixel length. Given as a matrix
    %   [y_size, x_size].
    % parallel_flag: logical, default = false. If set to true, image
    %   creation will take place across multiple cores. Parallel pool must
    %   already be activated.
    % calculation_cutoff_sigma: maximum number of sigmas that the psf
    %   calculation will be made to. Default = 5 
    % zero_cutoff: minimum probability density value requried to not be 
    %   counted as a zero in the final map. Default = 1E-6.
    %
    % Ouptuts:
    % image: Super-resolution image of the detected points. Each point has
    %   a total integrated signal of 1 and is represented by the value of
    %   the pdf at the center of each pixel (NOT the integrated value over 
    %   the area of the pixel; choose small enough pixels that this 
    %   distinction does not matter. In practice, this is about 1/3 the 
    %   sigma value for relative errors that are maxiamlly 5% off. By 1/10 
    %   the sigma value, the difference between the sampled and integrated 
    %   pdf is vanishinly small for all real purposes.) Returned as a
    %   sparse matrix.
    
    % Set defaults
    if nargin < 7; zero_cutoff = 1e-6; end 
    if nargin < 6; calc_cutoff_sigmas = 5; end
    if nargin < 5; parallel_flag = false; end 
        
    % Calc parameters for psfs of all the same size
    if isscalar(sigma)
        covar_matrix = [sigma.^2, 0; 0, sigma.^2];
        covar_inv = inv(covar_matrix);
        covar_det = det(covar_matrix);
        calc_cutoff_pixels = ceil((calc_cutoff_sigmas * sigma)/resolution);
    else % value does not matter, but we want these to exit
        covar_matrix = [];
        covar_inv = [];
        covar_det = [];
    end
   
    % Calc image parameters
    total_number_pixels_x = ceil(dims(2)/resolution);
    total_number_pixels_y = ceil(dims(1)/resolution);
    
    % Create nested function for adding a pdf to the image
    function full_image = make_image(data, covar_inv, covar_det, calc_cutoff_pixels)
        
        % Initalize image array and loop through each point
        full_image = zeros(total_number_pixels_y, total_number_pixels_x);
        for data_ind = 1:length(data.x)
        
            % Recalc the psf parameters if each sigma value is given independently
            if strcmp(sigma, 'per-point');
                covar_matrix = [data.sigma(data_ind).^2, 0; 0, data.sigma(data_ind).^2];
                covar_inv = inv(covar_matrix);
                covar_det = det(covar_matrix);
                calc_cutoff_pixels = ceil((calc_cutoff_sigmas * sigma)/resolution);
            end
            
            % Add pdf to the image
            full_image = add_pdf_to_image(full_image, data.x(data_ind), data.y(data_ind), covar_inv, covar_det, calc_cutoff_pixels);
        end
        % Calc index range for the point
        center_row = round(total_number_pixels_y - y/resolution + 0.5);
        center_column = round(x/resolution + 0.5);
        min_row = center_row - calc_cutoff_pixels;
        max_row = center_row + calc_cutoff_pixels;
        min_column = center_column - calc_cutoff_pixels;
        max_column = center_column + calc_cutoff_pixels;

        % Make sure that the point and the index range does not go off the 
        % edges of the image.
        if center_row < 1 || center_row > total_number_pixels_y; return; end
        if center_column < 1 || center_column > total_number_pixels_x; return; end
        if min_row < 1; min_row = 1; end
        if max_row > total_number_pixels_y; max_row = total_number_pixels_y; end
        if min_column < 1; min_column = 1; end
        if max_column > total_number_pixels_x; max_column = total_number_pixels_x; end

        % Calc mesh x and y values in original pixel units
        x_vector = resolution .* [min_column-0.5:1:max_column-0.5];
        y_vector = resolution .* [total_number_pixels_y-min_row+0.5:-1:total_number_pixels_y-max_row+0.5];
        [x_mesh, y_mesh] = meshgrid(x_vector, y_vector);
        X = [x_mesh(:).'; y_mesh(:).'];

        % Calc pdf
        X_shifted = X - repmat([x; y], 1, size(X, 2)); 
        pdf = (2 .* pi .* sqrt(covar_det))^-1 .* exp(-0.5 .* sum(X_shifted.' * covar_inv .* X_shifted.', 2));

        % Add to total
        full_image(min_row:max_row, min_column:max_column) =...
            full_image(min_row:max_row, min_column:max_column) + reshape(pdf, size(x_mesh));
    end
    
    % Evaluate on one core
    if ~parallel_flag
        
       
    
    % Evaluate in parallel
    elseif parallel_flag
        
        % Determine indices needed to split data list into equal parts for each cluster
        cluster_info = parcluster('local');
        num_workers = cluster_info.NumWorkers;
        total_datapoints = length(data.x);        
        current_value = 0;
        start_index = zeros(num_workers, 1);
        end_index = zeros(num_workers, 1);
        for worker_index = 1:num_workers
            start_index(worker_index) = current_value + 1;
            worker_datapoints = floor(total_datapoints/num_workers);
            if worker_index <= mod(total_datapoints, num_workers)
                worker_datapoints = worker_datapoints + 1;
            end
            end_index(worker_index) = current_value + worker_datapoints;
            current_value = end_index(worker_index);
        end
        
        % Split the data vector
        data_chunks = cell(num_workers);
        for chunk_index = 1:num_workers
            newdata = struct();
            newdata.x = data.x(start_index(chunk_index), end_index(chunk_index));
            newdata.y = data.y(start_index(chunk_index), end_index(chunk_index));
            if strcmp(sigma, 'per-point')
                 newdata.sigma = data.sigma(start_index(chunk_index), end_index(chunk_index));
            end
            data_chunks{chunk_index} = newdata;
        end
        
        % create nested function to evaluate a chunk        
        for data_ind = 1:length(data.x)
            F(data_ind) = parfeval(@rand,1,1,10); % One row each future
end
R = fetchOutputs(F); % 10-by-10 concatenated output
        
        % Initalize image array and loop through each point
        full_image = zeros(total_number_pixels_y, total_number_pixels_x);
        
        
            % Recalc the psf parameters if each sigma value is given independently
            if strcmp(sigma, 'per-point');
                sigma = data.sigma(data_ind);
                covar_matrix = [sigma.^2, 0; 0, sigma.^2];
                covar_inv = inv(covar_matrix);
                covar_det = det(covar_matrix);
                calc_cutoff_pixels = ceil((calc_cutoff_sigmas * sigma)/resolution);
            end
            
            % Add pdf to the image
            full_image = add_pdf_to_image(full_image, data.x(data_ind), data.y(data_ind), covar_inv, covar_det, calc_cutoff_pixels);
        end
        
    % Convert very low values to zero
    dense_image(dense_image < zero_cutoff) = 0;
    
    % Convert to sparse image to save space
    sparse_image = sparse(dense_image);
end
