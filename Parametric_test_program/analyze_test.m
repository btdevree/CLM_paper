function analyze_test(test_directory, output_directory)
%ANALYZE_TEST Compare responses of a parametric STORM image test to the 
% ground truth and save the summary statistics.
%
% Requires functions from Common, Actin_simulator, Border_simulator,
% Dots_simulator, and Region_simulator folders; make sure these are on the
% MATLAB path. 
%   
% Input:
%   test_directory: string, filepath for the test folder where the test
%       archive and response are kept. Optional, default = current 
%       directory.
%   output_directory: string, filepath for saving the test results.
%       Optional, default = test_directory.
% Output:
%   Summary is saved to the output directory as "test_summary.mat"

% Set defaults
if nargin < 1; test_directory = ''; end;
if nargin < 2; output_directory = test_directory; end;    

% Start up a new parallel pool using default settings, if needed.
gcp

% Expand and load the test, load the answers
untar([test_directory, '/test_archive.tar.gz'], output_directory);
load([test_directory, '/test_archive.mat']); % Loads test_info
load([test_directory, '/response_info.mat']); % Load responses

% Get required constants and variables
number_images = size(test_info.image_info, 1);
params = test_info.params;

% Determine coordinate bounds
min_x_bound = params.bounds(1);
min_y_bound = params.bounds(2);
max_x_bound = params.bounds(3);
max_y_bound = params.bounds(4);
x_length = max_x_bound - min_x_bound;
y_length = max_y_bound - min_y_bound;

% Initalize variables
ts = struct(); % For test_summary, shorten for convenience
ts.TCI_mean = cell(number_images, 1);
ts.ECI_mean = cell(number_images, 1);
ts.TCI_stdev = cell(number_images, 1);
ts.ECI_stdev = cell(number_images, 1);

% A template structure for the shared attributes that all the image type need to keep track of. 
template_struct = struct();
template_struct.contrast_ratios = [];
template_struct.ECI = [];
template_struct.TCI = [];

% ------Calculate the TCI and ECI of each image-------

% Define fraction of events, number of events, and replicates
fraction_vector = [0; .005; .01; .02; .04; .06; .08; .1; .15; .2; .25; .3; .35; .4; .5; .6; .7; .8; .9; 1];
number_pseudoreplicates = 30; % Number of times to split up each dataset

% Run IIC curve script for each image
for image_index = 1:number_images
    
    % Get the ideal image and dataset
    ideal_image = test_info.ideal_images{image_index};
    dataset = test_info.event_data{image_index};

    % Get the IIC and ideal discrepency results
    [IIC_results, ~, ideal_discrepency_results] = calculate_IIC(params, dataset, fraction_vector, number_pseudoreplicates, 'sum_of_squares', ideal_image, true);

    % Calculate the AOC and ECI of the sum of squares IIC curves
    dFrac = fraction_vector(2:end) - fraction_vector(1:end-1);
    midpoint_II = (IIC_results(2:end, :, :) + IIC_results(1:end-1, :, :))/2;
    AOC = sum(repmat(dFrac, 1, size(midpoint_II, 2), size(midpoint_II, 3)) .* midpoint_II, 1);
    ECI_data = (2 * AOC - 1);
    ts.ECI_mean{image_index} = squeeze(mean(ECI_data, 3))';
    ts.ECI_stdev{image_index} = squeeze(std(ECI_data, 0, 3))';

    % Calculate the corrosponding TCI 
    TCI_data = 1 - (ideal_discrepency_results(end, :, :) ./ ideal_discrepency_results(1, :, :));
    ts.TCI_mean{image_index} = squeeze(mean(TCI_data, 3))';
    ts.TCI_stdev{image_index} = squeeze(std(TCI_data, 0, 3))';
end

% Find indices for each type of image
region_indices = [];
dots_indices = [];
actin_indices = [];
border_indices = [];
for image_index = 1:number_images
    info = test_info.ideal_images{image_index};
    if strcmp(info.image_type, 'region')
        region_indices = [region_indices; image_index];
    elseif strcmp(info.image_type, 'dots')
        dots_indices = [dots_indices; image_index];
    elseif strcmp(info.image_type, 'actin')
        actin_indices = [actin_indices; image_index];
    elseif strcmp(info.image_type, 'border')
        border_indices = [border_indices; image_index];
    end
end

% -----Region analysis------

%Calc needed number of pixels
area_map_resolution = 3; % nanometers
num_pixels_x = ceil(x_length / map_resolution);
num_pixels_y = ceil(y_length / map_resolution);

 % Calc meshgrid
x_vec = single([0.5: 1: num_pixels_x - 0.5] .* area_map_resolution + min_x_bound); % add grid vector and origin coords
y_vec = single([num_pixels_y - 0.5: -1: 0.5] .* area_map_resolution + min_y_bound);
[xmesh, ymesh] = meshgrid(x_vec, y_vec);

% Initalize results variables
found_regions = template_struct;
found_regions.net_area = [];
found_regions.abs_area = [];
found_regions.true_area = [];
missing_regions = template_struct;

for image_index = region_indices;
    info = test_info.image_info{image_index};
    ground_truth = test_info.ground_truth_coords;
       
    if ~isempty(responses.x{image_index});
   
        % Create area_maps
        true_area_map = inpolygon(xmesh, ymesh, ground_truth(:, 1), ground_truth(:, 2));
        response_area_map = inpolygon(xmesh, ymesh, responses.x{image_index}, responses.y{image_index});
        extra_area_map = response_area_map & ~true_area_map;
        missing_area_map = true_area_map & ~response_area_map;
        
        % Calculate areas
        true_area = sum(true_area_map(:)) * 9; % nanometers squared
        extra_area = sum(extra_area_map(:)) * 9;
        missing_area = sum(missing_area_map(:)) * 9;      

        % Record areas and related information for the found regions
        found_regions.true_area = [found_regions.true_area; true_area];
        found_regions.abs_area = [found_regions.abs_area; extra_area + missing_area];
        found_regions.net_area = [found_regions.net_area; extra_area - missing_area];
        found_regions.contrast_ratios = [found_regions.contrast_ratios; info.contrast_ratio];
        found_regions.ECI = [found_regions.ECI; ts.ECI{image_index}];
        found_regions.TCI = [found_regions.TCI; ts.TCI{image_index}];
    
    else % No response, record information for the too-hard problems
        missing_regions.contrast_ratios = [missing_regions.contrast_ratios; info.contrast_ratio];
        missing_regions.ECI = [missing_regions.ECI; ts.ECI{image_index}];
        missing_regions.TCI = [missing_regions.TCI; ts.TCI{image_index}];
    end
end

% Copy into summary structure
ts.region.found_regions = found_regions;
ts.region.missing_regions = missing_regions;

% -------Dots analysis------------------

% Initalize results variables
found_points = template_struct;
found_points.distances = [];
missing_points = template_struct;
extra_points = template_struct;

dots_distance_cutoff = 500; 
for image_index = dots_indices;
    info = test_info.image_info{image_index};
    ground_truth = test_info.ground_truth_coords;
    
    % Get the distance between the response points and the nearest true coordinate
    if ~isempty(responses.x{image_index});
        distances = distance_to_bezier(ground_truth, responses.x{image_index}, responses.y{image_index}, true);
    else
        distances = [];
    end
    
    % If the distance is beyond the cutoff, throw it away and count it as an extra point
    extra_point_selector = distances > dots_distance_cutoff;
    distances = distances(~extra_point_selector);
    number_extra_points = sum(extra_point_selector);
   
    % Record distances and related information for the found points and the extra points
    found_points.distances = [found_points.distances; distances];
    found_points.contrast_ratios = [found_points.contrast_ratios; repmat([info.contrast_ratio], size(distances, 1), 1)];
    found_points.ECI = [found_points.ECI; repmat([ts.ECI{image_index}], size(distances, 1), 1)];
    found_points.TCI = [found_points.TCI; repmat([ts.TCI{image_index}], size(distances, 1), 1)];
    extra_points.contrast_ratios = [extra_points.contrast_ratios; repmat([info.contrast_ratio], number_extra_points, 1)];
    extra_points.ECI = [extra_points.ECI; repmat([ts.ECI{image_index}], number_extra_points, 1)];
    extra_points.TCI = [extra_points.TCI; repmat([ts.TCI{image_index}], number_extra_points, 1)];
    
    % Get the distance between the true coordinates and the nearest response to look for missed detections
    if ~isempty(responses.x{image_index});
        distances = distance_to_bezier([responses.x{image_index}, responses.y{image_index}], ground_truth(:, 1), ground_truth(:, 2), true);
        missing_point_selector = distances > dots_distance_cutoff;
        number_missing_points = sum(missing_point_selector);
    else
        number_missing_points = size(ground_truth, 1);
    end
    
    % Record distances and related information for the found points and the extra points
    missing_points.contrast_ratios = [missing_points.contrast_ratios; repmat([info.contrast_ratio], number_missing_points, 1)];
    missing_points.ECI = [missing_points.ECI; repmat([ts.ECI{image_index}], number_missing_points, 1)];
    missing_points.TCI = [missing_points.TCI; repmat([ts.TCI{image_index}], number_missing_points, 1)];
end

% Copy into summary structure
ts.dots.found_points = found_points;
ts.dots.extra_points = extra_points;
ts.dots.missing_points = missing_points;

% -------Actin lines analysis------------------

% NOTE: We are going to estimate the areas by taking a sum of rectangular 
% blocks as long as the minimum distance between one curve and any points
% on the the other curve. The estimate is the average with the same
% strategy performed starting with the other curve. Endpoint areas are 
% estimated as a triangle between the last point on one curve, the endpoint
% of the other curve, and the point on the other curve closest to the last
% point on the one curve. 

% Initalize results variables
found_lines = template_struct;
found_lines.area = [];
found_lines.true_length = [];
found_lines.type = [];
found_lines.size = [];
missing_lines = template_struct;
missing_lines.true_length = [];
missing_lines.type = [];
missing_lines.size = [];
extra_lines = template_struct;
extra_lines.true_length = [];
extra_lines.type = [];
extra_lines.size = [];

matching_area_cutoff = 100; % nanometers squared per nanometermeter length 
for image_index = actin_indices;
    info = test_info.image_info{image_index};
    ground_truth_x = test_info.ground_truth_coords(:, :, 1);
    ground_truth_y = test_info.ground_truth_coords(:, :, 2);
    
    % Get number of responses
    number_response_curves = size(response.x, 1);
    number_true_curves = size(ground_truth_x, 1);
    
    % Get the length and a bunch of points on each response curve
    response_arc_length = zeros(number_response_curves, 1);
    response_curve_points = cell(number_response_curves, 1);
    for curve_index = 1:number_response_curves

        % Construct the control_points matrix
        control_points = [response.x(curve_index, :)', response.y(curve_index, :)'];

        % Approximate the arc length
        approx_curve_points = calc_bezier_line(control_points, 25);
        approx_length = sum(sqrt(sum((approx_curve_points(2:end, :) - approx_curve_points(1:end-1, :)).^2 , 2)), 1);

        % Get points along the bezier curve 
        curve_points = calc_bezier_line(control_points, round(approx_length/3)); % Get curves with about 3 nm spaced points
        response_curve_points{curve_index} = curve_points;
        response_arc_length(curve_index) = sum(sqrt(sum((curve_points(2:end, :) - curve_points(1:end-1, :)).^2 , 2)), 1);
    end
    
    % Get the length and a bunch of points on each ground truth curve
    true_arc_length = zeros(number_true_curves, 1);
    true_curve_points = cell(number_true_curves, 1);
    for curve_index = 1:number_true_curves

        % Construct the control_points matrix
        control_points = [ground_truth_x(curve_index, :)', ground_truth_y(curve_index, :)'];

        % Approximate the arc length
        approx_curve_points = calc_bezier_line(control_points, 25);
        approx_length = sum(sqrt(sum((approx_curve_points(2:end, :) - approx_curve_points(1:end-1, :)).^2 , 2)), 1);

        % Get points along the bezier curve 
        curve_points = calc_bezier_line(control_points, round(approx_length/3)); % Get curves with about 3 nm spaced points
        true_curve_points{curve_index} = curve_points;
        true_arc_length(curve_index) = sum(sqrt(sum((curve_points(2:end, :) - curve_points(1:end-1, :)).^2 , 2)), 1);
    end
    
    % Get the distance between the response points and the nearest true coordinate
    if number_response_curves > 0 % No need to do any of this if there is no response

        % Get the closest point on the other curve to the endpoints of each
        % curve
        
        
    if ~isempty(responses.x{image_index});
        distances = distance_to_bezier(ground_truth, responses.x{image_index}, responses.y{image_index}, true);
    else
        distances = [];
    end
    
    % If the distance is beyond the cutoff, throw it away and count it as an extra point
    extra_point_selector = distances > dots_distance_cutoff;
    distances = distances(~extra_point_selector);
    number_extra_points = sum(extra_point_selector);
   
    % Record distances and related information for the found points and the extra points
    found_points.distances = [found_points.distances; distances];
    found_points.contrast_ratios = [found_points.contrast_ratios; repmat([info.contrast_ratio], size(distances, 1), 1)];
    found_points.ECI = [found_points.ECI; repmat([ts.ECI{image_index}], size(distances, 1), 1)];
    found_points.TCI = [found_points.TCI; repmat([ts.TCI{image_index}], size(distances, 1), 1)];
    extra_points.contrast_ratios = [extra_points.contrast_ratios; repmat([info.contrast_ratio], number_extra_points, 1)];
    extra_points.ECI = [extra_points.ECI; repmat([ts.ECI{image_index}], number_extra_points, 1)];
    extra_points.TCI = [extra_points.TCI; repmat([ts.TCI{image_index}], number_extra_points, 1)];
    
    % Get the distance between the true coordinates and the nearest response to look for missed detections
    if ~isempty(responses.x{image_index});
        distances = distance_to_bezier([responses.x{image_index}, responses.y{image_index}], ground_truth(:, 1), ground_truth(:, 2), true);
        missing_point_selector = distances > dots_distance_cutoff;
        number_missing_points = sum(missing_point_selector);
    else
        number_missing_points = size(ground_truth, 1);
    end
    
    % Record distances and related information for the found points and the extra points
    missing_points.contrast_ratios = [missing_points.contrast_ratios; repmat([info.contrast_ratio], number_missing_points, 1)];
    missing_points.ECI = [missing_points.ECI; repmat([ts.ECI{image_index}], number_missing_points, 1)];
    missing_points.TCI = [missing_points.TCI; repmat([ts.TCI{image_index}], number_missing_points, 1)];
end

% Copy into summary structure
ts.dots.found_points = found_points;
ts.dots.extra_points = extra_points;
ts.dots.missing_points = missing_points;


end
