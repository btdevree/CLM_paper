function image = create_STORM_image(data, resolution, sigma, dims, sparse_output_flag, parallel_flag, use_MEX_flag, calc_cutoff_sigmas, zero_cutoff)
% CREATE_STORM_IMAGE Creates a high-resolution image built from 
%   Gaussian pdfs. 
%
% Inputs:
% data: Structure with fields 'x' and 'y'. Assumes a coordinate system with
%   the orgin at the top left corner of the top left pixel, and values 
%   increase with value as one goes downward or rightward in the image. 
%   Given in original pixel values.
% resolution: Desired resolution in fractions of original pixel length.
% sigma: Desired psf width in fractions of original pixel length. 
% dims: Desired image size in original pixel length. Given as a matrix
%   [y_size, x_size].
% parallel_flag: logical, default = false. If set to true, image
%   creation will take place across multiple cores. Parallel pool must
%   already be activated.
% use_MEX_flage: logical, default = false. If set to true, the core image 
%   creation will be done with the Gauss_STORM_image_MEX command, which
%   needs to be pre-compiled from the Gauss_STORM_image_MEX.cpp source.
% calculation_cutoff_sigma: maximum number of sigmas that the psf
%   calculation will be made to. Default = 5 
% zero_cutoff: minimum probability value requried to not be 
%   counted as a zero in the final map. Default = 1E-9.
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
%   sparse or full matrix of floating point doubles.

% Set defaults
if nargin < 9; zero_cutoff = 1e-9; end 
if nargin < 8; calc_cutoff_sigmas = 5; end
if nargin < 7; use_MEX_flag = false; end
if nargin < 6; parallel_flag = false; end 
if nargin < 5; sparse_output_flag = true; end 

% Calc parameters for psf
covar_matrix = [sigma.^2, 0; 0, sigma.^2];
covar_inv = inv(covar_matrix);
covar_det = det(covar_matrix);
calc_cutoff_pixels = ceil((calc_cutoff_sigmas * sigma)/resolution);

% Calc image parameters
total_number_pixels_x = ceil(dims(2)/resolution);
total_number_pixels_y = ceil(dims(1)/resolution);

% Evaluate on one core
if ~parallel_flag

    % make_image defined below
    full_image = make_image(data, covar_inv, covar_det, calc_cutoff_pixels, sigma, resolution,...
        total_number_pixels_y, total_number_pixels_x, dims, use_MEX_flag);

% Evaluate in parallel
% NOTE - TO DO: Should change this to use spmd command.
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
    data_chunks = cell(num_workers, 1);
    for chunk_index = 1:num_workers
        newdata = struct();
        newdata.x = data.x(start_index(chunk_index):end_index(chunk_index));
        newdata.y = data.y(start_index(chunk_index):end_index(chunk_index));
        if strcmp(sigma, 'per-point')
             newdata.sigma = data.sigma(start_index(chunk_index):end_index(chunk_index));
        end
        data_chunks{chunk_index} = newdata;
    end

    % Evaluate the chunks in parallel
    image_chunks = cell(num_workers, 1);
    parfor chunk_index = 1:num_workers
        image_chunks{chunk_index} = make_image(data_chunks{chunk_index}, covar_inv, covar_det,...
            calc_cutoff_pixels, sigma, resolution, total_number_pixels_y, total_number_pixels_x, dims, use_MEX_flag);
    end

    % Add results together
    full_image = zeros(total_number_pixels_y, total_number_pixels_x);
    for chunk_index = 1:num_workers
        full_image = full_image + image_chunks{chunk_index};
    end
end

% Convert very low values to zero
full_image(full_image < zero_cutoff) = 0;

if sparse_output_flag
    % Convert to sparse image to save space
    image = sparse(full_image);
else
    image = full_image;
end
end

% Create function for adding a pdf to the image
function full_image = make_image(data, covar_inv, covar_det, calc_cutoff_pixels,...
    resolution, total_number_pixels_y, total_number_pixels_x, dims, use_MEX_flag)

% Translate data into Cartesian system and create xy_data matrix for the Gauss_STORM_image function 
data.y = dims(1) - data.y;
xy_data = [data.x, data.y];

% Create x and y pixel coordinate vectors
x_vector = [0.5 : 1 : total_number_pixels_x - 0.5] .* resolution;
y_vector = [0.5 : 1 : total_number_pixels_y - 0.5] .* resolution;

% Assume square cutoff box
calc_cutoff_pixels_y = calc_cutoff_pixels;
calc_cutoff_pixels_x = calc_cutoff_pixels;

% Get the full image from Gauss_STORM_image
if ~use_MEX_flag
    full_image = Gauss_STORM_image(xy_data, resolution, covar_inv, covar_det, calc_cutoff_pixels_x, calc_cutoff_pixels_y, x_vector, y_vector);
else
    % Type checking to make sure we don't end up giving the MEX file bad inputs and cause a seg-fault
    assert(ismatrix(xy_data) && size(xy_data, 2) == 2 && isdouble(xy_data));
    assert(isnumeric(resolution) && isscalar(resolution));
    assert(size(covar_inv) == [2, 2] && isdouble(covar_inv));
    assert(isnumeric(covar_det) && isscalar(covar_det));
    assert(isnumeric(calc_cutoff_pixels_x) && isscalar(calc_cutoff_pixels_x) && calc_cutoff_pixels_x >= 0);
    assert(isnumeric(calc_cutoff_pixels_y) && isscalar(calc_cutoff_pixels_y) && calc_cutoff_pixels_y >= 0);
    assert(ismatrix(x_vector) && isdouble(x_vector) && size(xy_data, 2) == 1);
    
    % Get the image with verified inputs
    full_image = Gauss_STORM_image_MEX(xy_data, resolution, covar_inv, covar_det, calc_cutoff_pixels_x, calc_cutoff_pixels_y, x_vector, y_vector);
end
end
