%function [test_movie, test_timestamp, test_STORM_image_ch1, test_STORM_image_ch2] = create_test_movie_dv(parameters_struct, seed)
function create_test_movie_dv(parameters_struct, seed)
%CREATE_TEST_MOVIE Creates a test movie with the given parameters. Assumes
%   a dualview configuration.

%   Creates a test movie with the given parameters.
%   Inputs:
%   parameters_struct: structure of parameters made according to the
%       template created by the test_movie_parameters function. Please see 
%       that function for more documentation.
%   seed: random seed to start image creation from. Optional. If seed is 
%       empty or not given, matlab's native random seeding method will be 
%       used. 
%   Outputs:
%   test_movie: Tiff stack of frames of the test movie.
%   test_timestamp: Veatch lab simulated timestamp file for the movie.
%   test_STORM_image: A ground-truth image of what the reconstructed STORM
%       image should look like. 

% Set defaults
if nargin < 2; seed = []; end;

% Initialize the random number generator
if ~isempty(seed)
    rng(seed);
else
    rng('default');
end

% Rename parameters for convenience
params = parameters_struct; 

% Determine coordinate bounds
min_x_bound = params.bounds(1);
min_y_bound = params.bounds(2);
max_x_bound = params.bounds(3);
max_y_bound = params.bounds(4);

% ----- Create event distribution in cell------

% Determine cell center
if strcmp(params.cell_center, 'centered')
    cell_center_x = (max_x_bound - min_x_bound)/2;
    cell_center_y = (max_y_bound - min_y_bound)/2;
else
    cell_center_x = params.cell_center(1);
    cell_center_y = params.cell_center(2);
end

% Determine cell radius
if isnumeric(params.cell_radius)
    cell_radius = params.cell_radius;
else
    tokens = regexp(params.cell_radius, '(\d+)([%])', 'tokens'); % Look for a percentage value
    if ~isempty(tokens)
        max_radius = min([max_x_bound - cell_center_x, cell_center_x - min_x_bound,...
            max_y_bound - cell_center_y, cell_center_y - min_y_bound]);
        cell_radius = (str2num(tokens{1}{1})/100) * max_radius;
    end
end

% Place channel one events 
% Random distribution
if strcmp(params.ch1_autocor, 'random')
     ch1_event_coords = random_events_in_cell(params.number_events_ch1, [cell_center_x, cell_center_y], cell_radius);
end

% Assign channel two events to channel one events
% Initialize assignment vector; any assignments to remaining NaN values will be randomly placed in the cell
ch2_assignments = NaN(params.number_events_ch2, 1); 

% Random assignment
if strcmp(params.ch2_distribution, 'random')
    
    % Read parameters
    fraction_assigned = params.ch2_distribution_params(1);
    
    % Loop through each ch2 assignment value
    for index = 1:length(ch2_assignments)
        
        % Assign the value to a ch1 event at the requested probability
        test_value = rand;
        if test_value <= fraction_assigned
            ch2_assignments(index) = ceil(rand*size(ch1_event_coords, 1));
        end
    end

% Even assignment
elseif strcmp(params.ch2_distribution, 'evenly_distributed')
    
    % Get parameters
    max_assignment = params.ch2_distribution_params(1);
    number_ch1 = params.number_events_ch1;
    number_ch2 = params.number_events_ch2;
    
    % Calculate how many full and partial assignments need to be made
    number_assignments = min([max_assignment, number_ch2/number_ch1]);
    number_full_cycles = floor(number_assignments);
    fractional_cycle = number_assignments - number_full_cycles;
        
    % Assign full cycles
    for index = 1:number_full_cycles
        ch2_assignments((index-1)*number_ch1+1:index*number_ch1) = randperm(number_ch1);
    end
    
    % Assign partial cycles
    number_partial_assignments = round(fractional_cycle*number_ch1);
    ch2_assignments(number_ch1*number_full_cycles+1:number_ch1*number_full_cycles+number_partial_assignments) =...
        randperm(number_ch1, number_partial_assignments);
end

% Place ch2 events with specified crosscorrelation to ch1 events
% Create ch2 coordinate matrix
ch2_event_coords = zeros(params.number_events_ch2, 2);

% Random assignment
if strcmp(params.ch2_crosscor, 'random')
    % Reassign all vales to NaN so they will be randomly placed at the end of the crosscorrelation code block
    ch2_assignments = NaN(params.number_events_ch2, 1);

% Exact placement
elseif strcmp(params.ch2_crosscor, 'exact')
    
    % Read parameters and initialize 
    delta_r = params.ch2_distribution_params(1); % Distance between ch1 and ch2 event
    
    % Loop through each assignment value and create the appropriate event coordinate
    for index = 1:size(ch2_assignments, 1)
        
        % Test for NaN value, skip assignment if true
        if isnan(ch2_assignments(index))
            continue
        else
            
        % Generate a random angle of displacement and calculate coordinate shift
        angle = rand*2*pi;
        delta_x = cos(angle)*delta_r;
        delta_y = sin(angle)*delta_r;
        
        % Write down shifted coordinates
        parent_coord = ch1_event_coords(ch2_assignments(index), :);
        ch2_event_coords(index, :) = parent_coord + [delta_x, delta_y];
        end
    end
    
% Gaussian placement
elseif strcmp(params.ch2_crosscor, 'Gaussian')
         
    % Read parameters and initialize 
    mu = params.ch2_distribution_params(1); % Average distance between ch1 and ch2 event
    sigma = params.ch2_distribution_params(2); % Standard deviation of ch2 event distribution from the mean distance 
    
    % Loop through each assignment value and create the appropriate event coordinate
    for index = 1:size(ch2_assignments, 1)
        
        % Test for NaN value, skip assignment if true
        if isnan(ch2_assignments(index))
            continue
        else
            
        % Generate a random angle of displacement and calculate coordinate shift
        angle = rand*2*pi;
        delta_x = cos(angle)*mu;
        delta_y = sin(angle)*mu;
        
        % Write down shifted coordinates
        parent_coord = ch1_event_coords(ch2_assignments(index), :);
        ch2_event_coords(index, :) = parent_coord + [delta_x, delta_y] + mvnrnd([0,0], [sigma^2, 0; 0, sigma^2]);
        end
    end
end

% Randomly place any coords with NaN assignments
% Generate random coords
number_ch2_random_coords = sum(isnan(ch2_assignments));
ch2_random_coords = random_events_in_cell(number_ch2_random_coords, [cell_center_x, cell_center_y], cell_radius);

% Place random coords in appropriate slots
random_index = 1;
for coord_index = 1:size(ch2_assignments, 1)
    
    % Test for NaN value, assign if true and increment counter index
    if isnan(ch2_assignments(coord_index))
        ch2_event_coords(coord_index, :) = ch2_random_coords(random_index, :);
        random_index = random_index + 1;
    end
end

% --------- Create temporal length distributions ----------- 

% Calculate start and stop times  
% Ch1 and ch2 are the same, create nested function for re-use.
function [event_times] = calc_event_times(event_coords_size)
    
    % Initialize start and stop time matrices (nx2 doubles; each row = [start, stop])
    event_times = zeros(event_coords_size);
    
    % Random start times from t=-(expected event length) to t=(end of movie) in frames
    if strcmp(params.event_length_model, 'random')
        event_times(:, 1) = rand(size(event_times, 1), 1) .*...
            (params.number_frames + params.event_length_model_params(1)) - params.event_length_model_params(1);

    % Discrete start times from t=-(expected event length) to t=(end of movie - 1) in frames
    elseif strcmp(params.event_length_model, 'discrete')
        event_times(:, 1) = floor(rand(size(event_times, 1), 1) .*...
            (params.number_frames + params.event_length_model_params(1)) - params.event_length_model_params(1));

    % Discrete start times from t=0 to t=(end of movie - 1) in frames
    elseif strcmp(params.event_length_model, 'one_frame')
        event_times(:, 1) = floor(rand(size(event_times, 1), 1) .* params.number_frames);
    end

    % Calculate end times
    % Random end times from an exponential distribution
    if strcmp(params.event_length_model, 'random')
        lengths = exprnd(params.event_length_model_params(1), size(event_times, 1), 1);
        event_times(:, 2) = event_times(:, 1) + lengths;
        
    % Discrete end times from a Poission distribution
    elseif strcmp(params.event_length_model, 'discrete')
        lengths = poissrnd(params.event_length_model_params(1), size(event_times, 1), 1);
        event_times(:, 2) = event_times(:, 1) + lengths;
        
    % One frame only
    elseif strcmp(params.event_length_model, 'one_frame_only')
        event_times(:, 2) = event_times(:, 1) + 1;
    end
end

% Run the nested function to get the start and stop times
ch1_event_times = calc_event_times(size(ch1_event_coords));
ch2_event_times = calc_event_times(size(ch2_event_coords));

% ------ Convert to per-frame lists of positions and photons --------

% Calculate required constants
 if strcmp(params.event_length_model, 'random') || strcmp(params.event_length_model, 'discrete')
    ch1_photons_per_frame = params.ch1_photon_budget_params(1) / params.event_length_model_params(1);
    ch2_photons_per_frame = params.ch2_photon_budget_params(1) / params.event_length_model_params(1);
 elseif strcmp(params.event_length_model, 'one_frame')
     ch1_photons_per_frame = params.ch1_photon_budget_params(1);
     ch2_photons_per_frame = params.ch2_photon_budget_params(1);
 end

% Loop through each event and split it into the per-frame information
% Ch1 and ch2 are the same, create nested function for re-use.
function [frame_list] = split_events_to_frames(event_coords, event_times, photons_per_frame) 

    % Initialize frame lists
    frame_list = cell(params.number_frames, 1);
   
    % Loop through each event
    for event_index = 1:size(event_coords, 1)
        
        % Initialize frame and illumination fraction matrices
        frame_numbers = [];
        fractions = [];
        
        % Get frame numbers and illumination fraction 
        start_time = event_times(event_index, 1);
        stop_time = event_times(event_index, 2);
        start_frame = max(floor(start_time), 1);
        stop_frame = min(floor(stop_time), params.number_frames);
        for frame_index = start_frame:stop_frame
            frame_numbers(end + 1, 1) = frame_index;
            fractions(end + 1, 1) = min([frame_index + 1 - start_time, stop_time - frame_index, stop_time - start_time, 1]);
        end
        
        % Add to frame_list
        for frame_num_index = 1:size(frame_numbers, 1)
            
            % Get/calculate the required per-frame info 
            x = event_coords(event_index, 1);
            y = event_coords(event_index, 2);
            photons = photons_per_frame * fractions(frame_num_index);
            
            % Don't add anything if there are no photons
            if photons == 0
                continue
                
            % Make a [x, y, photons] matrix if there is nothing
            elseif isempty(frame_list{frame_numbers(frame_num_index)})
                frame_list{frame_numbers(frame_num_index)} = [x, y, photons];
                
            % Add to the [x, y, photons] matrix if exists
            else
                info_matrix = frame_list{frame_numbers(frame_num_index)};
                info_matrix(end + 1, :) = [x, y, photons];
                frame_list{frame_numbers(frame_num_index)} = info_matrix;
            end
        end
    end
end

% Run nested function on both channels
frame_event_list_ch1 = split_events_to_frames(ch1_event_coords, ch1_event_times, ch1_photons_per_frame);
frame_event_list_ch2 = split_events_to_frames(ch2_event_coords, ch2_event_times, ch2_photons_per_frame);

% ------------- Create video and add background noise --------------

% Allocate empty matrix
video_column_height = ceil((params.bounds(4) - params.bounds(2))/params.movie_pixel_size);
video_row_width = ceil((params.bounds(3) - params.bounds(1))/params.movie_pixel_size);
video_data_ch1 = zeros(video_column_height, video_row_width, params.number_frames, 'uint16');
video_data_ch2 = zeros(video_column_height, video_row_width, params.number_frames, 'uint16');

% Add background to video
if strcmp(params.background_noise_model, 'Gaussian')
    
    % Create random numbers; casting to uint16 automatically zeros any negative numbers
    noise_ch1 = uint16(normrnd(params.background_noise_params_ch1(1), params.background_noise_params_ch1(2), size(video_data_ch1)));
    noise_ch2 = uint16(normrnd(params.background_noise_params_ch2(1), params.background_noise_params_ch2(2), size(video_data_ch2)));
    
    % Add noise to video data
    video_data_ch1 = video_data_ch1 + noise_ch1;
    video_data_ch2 = video_data_ch2 + noise_ch2;
    clear noise_ch1 noise_ch2 % Free up memory

elseif strcmp(params.background_noise_model, 'Poisson')
    
    % Create random numbers
    noise_ch1 = uint16(poissrnd(params.background_noise_params_ch1(1), size(video_data_ch1)));
    noise_ch2 = uint16(poissrnd(params.background_noise_params_ch2(1), size(video_data_ch2)));
    
    % Add noise to video data
    video_data_ch1 = video_data_ch1 + noise_ch1;
    video_data_ch2 = video_data_ch2 + noise_ch2;
    clear noise_ch1 noise_ch2 % Free up memory
end

% Make a cell-sized mask
x_vector = [0.5:1:size(video_data_ch1, 2)-0.5] * params.movie_pixel_size;
y_vector = [size(video_data_ch1, 1)-0.5:-1:0.5] * params.movie_pixel_size; 
[x_mesh, y_mesh] = meshgrid(x_vector, y_vector);
x_cellshift = x_mesh - cell_center_x;
y_cellshift = y_mesh - cell_center_y;
cell_mask = uint16(sqrt(x_cellshift.^2 + y_cellshift.^2) <= cell_radius);

% Add cellular background to video
if strcmp(params.cell_noise_model, 'Gaussian')
    
    % Create random numbers; casting to uint16 automatically zeros any negative numbers
    noise_ch1 = uint16(normrnd(params.cell_noise_params_ch1(1), params.cell_noise_params_ch1(2), size(video_data_ch1)));
    noise_ch2 = uint16(normrnd(params.cell_noise_params_ch2(1), params.cell_noise_params_ch2(2), size(video_data_ch2)));
    
    % Multiply by the mask to get the cell noise
    noise_ch1 = noise_ch1 .* repmat(cell_mask, 1, 1, size(video_data_ch1, 3));
    noise_ch2 = noise_ch2 .* repmat(cell_mask, 1, 1, size(video_data_ch2, 3));
    
    % Add noise to video data
    video_data_ch1 = video_data_ch1 + noise_ch1;
    video_data_ch2 = video_data_ch2 + noise_ch2;
    clear noise_ch1 noise_ch2 % Free up memory
    
elseif strcmp(params.cell_noise_model, 'Poisson')
    
    % Create random numbers; casting to uint16 automatically zeros any negative numbers
    noise_ch1 = uint16(poissrnd(params.cell_noise_params_ch1(1), size(video_data_ch1)));
    noise_ch2 = uint16(poissrnd(params.cell_noise_params_ch2(1), size(video_data_ch2)));
    
    % Multiply by the mask to get the cell noise
    noise_ch1 = noise_ch1 .* repmat(cell_mask, 1, 1, size(video_data_ch1, 3));
    noise_ch2 = noise_ch2 .* repmat(cell_mask, 1, 1, size(video_data_ch2, 3));
    
    % Add noise to video data
    video_data_ch1 = video_data_ch1 + noise_ch1;
    video_data_ch2 = video_data_ch2 + noise_ch2;
    clear noise_ch1 noise_ch2 % Free up memory
end

% ----------- Add events to the video --------------

% Create nested function to re-use code for each channel
function [video_data] = add_events_to_video(frame_event_list, video_data, psf_params)
    
    % Initialize a double matrix to hold a frame of the video
    current_frame = zeros(size(video_data, 1), size(video_data, 2));
    
    % Determine the size of the area that needs to be calculated; both Gaussian and Airy methods use the same cutoff calculation
    calculation_pixel_radius = ceil((psf_params(1) * psf_params(2))/params.movie_pixel_size);
    
    % Loop through each frame
    for frame_index = 1:params.number_frames

        % Zero out the current_frame
        current_frame(:) = 0;
        
        % Rename event for convenience
        event_data = frame_event_list{frame_index};
        
        % Loop through each column in the frame event list
        for event_index = 1:size(event_data, 1);
            
            % Rename event for convenience 
            current_event = event_data(event_index, :);
           
            % Get the index values for the region to be calculated
            center_row_ind = size(current_frame, 1) + 1 - ceil(current_event(2)/params.movie_pixel_size); % need to flip axis to line up with indexing
            center_column_ind = ceil(current_event(1)/params.movie_pixel_size); 
            top_row_ind = center_row_ind - calculation_pixel_radius;
            bottom_row_ind = center_row_ind + calculation_pixel_radius;
            left_column_ind = center_column_ind - calculation_pixel_radius;
            right_column_ind = center_column_ind + calculation_pixel_radius;
            
            % Make sure index values don't go off the edges of the frame
            top_row_ind = max(1, top_row_ind);
            bottom_row_ind = min(size(current_frame, 1), bottom_row_ind);
            left_column_ind = max(1, left_column_ind);
            right_column_ind = min(size(current_frame, 2), right_column_ind);
            
            % Cut out appropreate pieces of meshgrid coordinates
            x_coords = x_mesh(top_row_ind:bottom_row_ind, left_column_ind:right_column_ind);
            y_coords = y_mesh(top_row_ind:bottom_row_ind, left_column_ind:right_column_ind);
            
            % Function call to calc matrix of pixel centers.
            if strcmp(params.psf_model, 'Gaussian')
                psf_values = calc_gaussian_psf(psf_params(1), current_event, x_coords, y_coords, params.movie_pixel_size, params.resampling_model, params.resampling_model_params);
            elseif strcmp(prams.psf_model, 'Airy')
                psf_values = calc_airy_psf(psf_params(1), current_event, x_coords, y_coords);
            end
            
            % Add psf to the image frame
            current_frame(top_row_ind:bottom_row_ind, left_column_ind:right_column_ind) =...
                current_frame(top_row_ind:bottom_row_ind, left_column_ind:right_column_ind) + psf_values;
        end
        
        % Add frame to video data
        video_data(:, :, frame_index) = video_data(:, :, frame_index) + uint16(current_frame);
    end
end

% Call function to add event data
video_data_ch1 = add_events_to_video(frame_event_list_ch1, video_data_ch1, params.psf_model_params_ch1);
video_data_ch2 = add_events_to_video(frame_event_list_ch2, video_data_ch2, params.psf_model_params_ch2);
     
% ----------- Add detection and camera effects to image -----------------

% Convert photon values in video to electrons entering the A/D converter
% Create nested function to re-use code for each channel
function [video_data] = add_EM_gain(video_data)

    % Loop frame by frame to avoid creating a full-size video in double precision floats (would be about 1/2 a gigabyte each)
    for frame_index = 1:params.number_frames          
    
        % Convert the frame to doubles
        current_frame = double(video_data(:, :, frame_index));
        
        % Multiply signal according to model
        if strcmp(params.EM_gain_model, 'Gaussian')
       
            % Calculate parameters
            mu = params.EM_gain .* current_frame;
            sigma = sqrt(2) .* params.EM_gain .* sqrt(current_frame); % sqrt(2) = excess noise factor at large (>10) gain settings
            
            % Sample the distribution to get the electrons after amplification
            current_frame = normrnd(mu, sigma);
            
        elseif strcmp(params.EM_gain_model, 'Gamma')
            
            % Sample the distribution to get the electrons after amplification
            current_frame = gamrnd(current_frame, params.EM_gain); 
        end
            
        % Round to integers and replace frame in video data
        video_data(:, :, frame_index) = uint16(current_frame);       
    end
end

% Run EM gain function for each channel
video_data_ch1 = add_EM_gain(video_data_ch1);
video_data_ch2 = add_EM_gain(video_data_ch2);

% Convert electron counts into ADU (grey levels)
% Create nested function to re-use code for each channel
function [video_data] = convert_to_ADU(video_data)

    % Loop frame by frame to avoid creating a full-size video
    for frame_index = 1:params.number_frames          
    
        % Convert the frame to doubles
        current_frame = double(video_data(:, :, frame_index));
        
        % Add read noise to the frame
        current_frame = normrnd(current_frame, params.A_to_D_noise);
        
        % Convert to ADU
        current_frame = current_frame ./ params.e_per_ADU;
        
        % Round to integers and replace frame in video data
        video_data(:, :, frame_index) = uint16(current_frame);     
    end
end

% Run ADU conversion function for each channel
video_data_ch1 = convert_to_ADU(video_data_ch1);
video_data_ch2 = convert_to_ADU(video_data_ch2);

% Add baseline to all pixels
video_data_ch1 = video_data_ch1 + uint16(params.baseline_clamp);
video_data_ch2 = video_data_ch2 + uint16(params.baseline_clamp);

% -------------- Save movie ----------------

% Call save_movie local function
save_movie(video_data_ch1, video_data_ch2, 'testmovie');

% Free up memory
clear video_data_ch1 video_data_ch2

% ----------- Create STORM images --------------

% Parameters needed for STORM image creating function
STORM_resolution = params.STORM_pixel_size; % data is given in nm, so resolution is just the nm per pixel
STORM_sigma = params.STORM_precision; % data is given in nm, so precision is in nm
STORM_dims = [max_y_bound - min_y_bound, max_x_bound - min_x_bound]; % data is given in nm, so bounds are in nm

% Create nested function to re-use code for each channel
function [STORM_image] = calc_STORM_image(STORM_resolution, STORM_sigma, STORM_dims, event_coords)

    % Convert points to a data structure
    data_struct = struct();
    data_struct.x = event_coords(:, 1);
    data_struct.y = event_coords(:, 2);

    % Call image generating function, use parallel processing if possible
    STORM_image = create_STORM_image(data_struct, STORM_resolution, STORM_sigma, STORM_dims, true); % note that output is spars
end

% Run STORM image function for each channel
ch1_STORM_image = calc_STORM_image(STORM_resolution, STORM_sigma, STORM_dims, ch1_event_coords);
ch2_STORM_image = calc_STORM_image(STORM_resolution, STORM_sigma, STORM_dims, ch2_event_coords);

% Get a two-color storm image
ch1_color = params.ch1_color;
ch2_color = params.ch2_color;

% Create a display RGB image
STORM_RGB_image = convert_STORM_image_to_RGB({ch1_STORM_image, ch2_STORM_image}, {ch1_color, ch2_color}, 3, 'channel_max');

% Save images as mat files
save('testmovie_STORM_images', 'ch1_STORM_image', 'ch2_STORM_image', 'STORM_RGB_image');

end

function [coords] = random_events_in_cell(number_events, cell_center, cell_radius)
% Creates the specified number of randomly places events inside the cell

%   Inupts:
%   number_events: positive integer, the number of events to generate.
%   cell_center: 1x2 floating-point double, center coordinate of the cell
%       given as [x, y] in nanometers.
%   cell_radius: floating-point double, radius of the cell in nanometers.
%   Output:
%   coords: nx2 floating-point double array of coordinates.

% Initialize variables
coords = zeros(number_events, 2);
event_counter = 0;

% Repeat loop until enough random events inside the circle are generated
while event_counter < number_events 
    event_coord = (rand(1,2) * cell_radius * 2) - cell_radius;
    if sqrt(sum(event_coord.^2)) <= cell_radius
       event_counter = event_counter + 1;
       coords(event_counter, :) = event_coord + cell_center;
    end
end
end

function [psf_values] = calc_gaussian_psf(psf_sigma, current_event, centered_x_coords, centered_y_coords, movie_pixel_size, resampling_model, resampling_model_params)
% Calculates a Gaussian psf and adds resampling noise

%   Inputs:
%   psf_sigma: zero-mean Gaussian parameters, given as [standard deviation(sigma)] in nanometers.
%   current_event: double vector given as [x, y, photons] with x and y in nanometers
%   centered_x_coords/centered_y_cords: coordinates in meshgrid format at the center of each pixel, given in nanometers. 
%   movie_pixel_size: the resolution of the pixels in the movie, given in nanometers.
%   resampling_model: string containing the choice of resampling model
%   resampling_model_params: matrix containing the required parameters for the resampling model
%   Outputs:
%   psf_values: the psf evaluated at each pixel, with resampling noise added

% Shift meshgrid values down by 0.5 pixels so they represent edges of pixels, not centers
x_coords = zeros(size(centered_x_coords)+1);
x_coords(2:end, 1:end-1) = centered_x_coords - 0.5 * movie_pixel_size;
y_coords = zeros(size(centered_y_coords)+1);
y_coords(2:end, 1:end-1) = centered_y_coords - 0.5 * movie_pixel_size;

% Add right and top edge to edge coordinate matrix
x_coords(1, :) = x_coords(2, :);
x_coords(:, end) = x_coords(:, end-1) + movie_pixel_size;
y_coords(1, :) = y_coords(2, :) + movie_pixel_size;
y_coords(:, end) = y_coords(:, end-1);

% Shift center coords if using Gaussian displacement for resampling
if strcmp(resampling_model, 'Gaussian_displacement')
    current_event(1) = current_event(1) + normrnd(0, resampling_model_params);
    current_event(2) = current_event(2) + normrnd(0, resampling_model_params);
end

% Evaluate the Gaussian cdf at each coordinate
cdf_values = reshape(mvncdf([x_coords(:), y_coords(:)], current_event(1:2), [psf_sigma.^2, 0; 0, psf_sigma.^2]), size(x_coords));

% Calculate the integrated gaussian over each pixel area
a = cdf_values(1:end-1, 1:end-1);
b = cdf_values(1:end-1, 2:end);
c = cdf_values(2:end, 1:end-1);
d = cdf_values(2:end, 2:end);
psf_values = (b + c) - (a + d);

% Multiply by the number of total detected photons for the frame
psf_values = psf_values * current_event(3);

% Corrupt with Poisson shot noise for resampling
if strcmp(resampling_model, 'shot_noise')
    psf_values = poissrnd(psf_values);
end
end

function save_movie(data_ch1, data_ch2, filename)
% SAVE_MOVIE Save the movie and timestamp info

%   Inputs:
%   data_ch1/data_ch2: 3D arrays of channel 1/2 data in uint16 format.
%   filename: string containing the filename that the move is saved under.
%       Do not include the .tiff file specification.
%   Outputs: 
%   Movie is saved in current directory under the indicated filename. The
%   associated timestamp file is saved as an m file with the name filename 
%       plus '_timestamp.m'

% Set up Tiff tag information
tagstruct = struct();
tagstruct.Photometric = Tiff.Photometric.MinIsBlack;
tagstruct.ImageLength = size(data_ch1, 1) + size(data_ch2, 1); 
tagstruct.ImageWidth = size(data_ch1, 2);
tagstruct.SamplesPerPixel = 1;
tagstruct.RowsPerStrip = size(data_ch1, 1) + size(data_ch2, 1);
tagstruct.PlanarConfiguration = Tiff.PlanarConfiguration.Chunky;
tagstruct.Compression = Tiff.Compression.None;
tagstruct.SampleFormat = Tiff.SampleFormat.UInt;
tagstruct.BitsPerSample = 16;

% Create Tiff object
tiff_obj = Tiff([filename, '.tif'], 'w'); 

% Loop through each frame
depth = size(data_ch1, 3);
for frame_index = 1:depth
    current_frame = [data_ch1(:, :, frame_index); data_ch2(:, :, frame_index)];
    tiff_obj.setTag(tagstruct);
    tiff_obj.write(current_frame);
    if frame_index ~= depth;
        tiff_obj.writeDirectory();
    end
end

% close Tiff object
tiff_obj.close();

% Setup timestamp info
start_time = clock;
Nframes = depth;
cropdims = [1, size(data_ch1, 2), 1, size(data_ch1, 1) + size(data_ch2, 1)];

% Save timestamp
save([filename, '_timestamp.mat'], 'start_time', 'Nframes', 'cropdims');
end

