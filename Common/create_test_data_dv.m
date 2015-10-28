function [ch1_data, ch2_data, movie_variables_structure ] = create_test_data_dv(parameters_struct, seed)
%CREATE_TEST_DATA_DV Creates data points for a test movie with the given 
%   parameters. Assumes a dualview configuration.

%   Creates a test movie with the given parameters.
%   Inputs:
%   parameters_struct: structure of parameters made according to the
%       template created by the test_movie_parameters function. Please see 
%       that function for more documentation.
%   seed: random seed to start image creation from. Optional. If seed is 
%       empty or not given, matlab's native random seeding method will be 
%       used. 
%   Output:
%   ch1_data/ch2_data: datapoints for use in making a test movie and/or 
%       test STORM images

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

% ------- Prepare output data ------------

% Rename event coords
ch1_data = ch1_event_coords;
ch2_data = ch2_event_coords;

% Collect required variables for movie
s = struct();
s.cell_center_x = cell_center_x;
s.cell_center_y = cell_center_y; 
s.cell_radius = cell_radius;
movie_variables_structure = s;
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