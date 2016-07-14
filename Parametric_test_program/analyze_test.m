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
%untar([test_directory, '/test_archive.tar.gz'], output_directory);
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

fraction_vector = [0; .02; .06; .15; .3; .6; 1];
number_pseudoreplicates = 3; % Number of times to split up each dataset

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
% point on the one curve. No attempt is made to combine several ground
% truth curves in the case that intersecting curves leads to tracing along
% multiple different tracks in the same response. 

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
matching_endpoint_cutoff = 500; % nanometers
for image_index = actin_indices;
    info = test_info.image_info{image_index};
    ground_truth_x = test_info.ground_truth_coords(:, :, 1);
    ground_truth_y = test_info.ground_truth_coords(:, :, 2);
    found_true_curve_indices = [];
    
    % Get number of responses
    number_response_curves = size(response.x{image_index}, 1);
    number_true_curves = size(ground_truth_x, 1);
    
    % Get the length and a bunch of points on each response curve
    response_arc_length = zeros(number_response_curves, 1);
    response_curve_points = cell(number_response_curves, 1);
    for curve_index = 1:number_response_curves

        % Construct the control_points matrix
        control_points = [response.x{image_index}(curve_index, :)', response.y{image_index}(curve_index, :)'];

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
        
        % Compute areas and endpoint distances for each response-truth curve pair.
        area_estimates = inf(number_response_curves, number_true_curves); % inf makes sure that undetermined areas evaluate higher than any calculated value
        endpoint_1_distances = inf(number_response_curves, number_true_curves);
        endpoint_2_distances = inf(number_response_curves, number_true_curves);
        for response_curve_index = 1:number_response_curves
            for true_curve_index = 1:number_true_curves
                 
                % Get endpoints and the closest point to on the other curve to the endpoints
                true_endpoints = [true_curve_points{true_curve_index}(1, :); true_curve_points{true_curve_index}(end, :)];
                response_endpoints = [response_curve_points{response_curve_index}(1, :); response_curve_points{response_curve_index}(end, :)];
                [true_closepoints, true_closepoint_indices] = distance_to_bezier(true_curve_points{true_curve_index}, response_endpoints, true);
                [response_closepoints, response_closepoint_indices] = distance_to_bezier(response_curve_points{response_curve_index}, true_endpoints, true);
                
                % Get endpoint to endpoint distances
                r1_t1_dist = sqrt(sum((true_endpoints(1, :) - response_endpoints(1, :)).^2, 2));
                r1_t2_dist = sqrt(sum((true_endpoints(2, :) - response_endpoints(1, :)).^2, 2));
                r2_t2_dist = sqrt(sum((true_endpoints(2, :) - response_endpoints(2, :)).^2, 2));
                r2_t1_dist = sqrt(sum((true_endpoints(1, :) - response_endpoints(2, :)).^2, 2));
                
                % Match endpoints
                if r1_t1_dist < r1_t2_dist
                    endpoint_1_distances(response_curve_index, true_curve_index) = r1_t1_dist;
                    end_group_1 = [response_endpoints(1, :); true_endpoints(1, :); response_closepoints(1, :); true_closepoints(1, :)];
                    r1_index = response_closepoint_indices(1);
                    t1_index = true_closepoint_indices(1);
                else
                    endpoint_1_distances(response_curve_index, true_curve_index) = r1_t2_dist;
                    end_group_1 = [response_endpoints(1, :); true_endpoints(2, :); response_closepoints(2, :); true_closepoints(1, :)];
                    r1_index = response_closepoint_indices(2);
                    t1_index = true_closepoint_indices(1);
                end
                if r2_t2_dist < r2_t1_dist
                    endpoint_2_distances(response_curve_index, true_curve_index) = r2_t2_dist;
                    end_group_2 = [response_endpoints(2, :); true_endpoints(2, :); response_closepoints(2, :); true_closepoints(2, :)];
                    r2_index = response_closepoint_indices(2);
                    t2_index = true_closepoint_indices(2);
                else
                    endpoint_2_distances(response_curve_index, true_curve_index) = r2_t1_dist;
                    end_group_2 = [response_endpoints(2, :); true_endpoints(1, :); response_closepoints(1, :); true_closepoints(2, :)];
                    r1_index = response_closepoint_indices(1);
                    t1_index = true_closepoint_indices(2);
                end
                
                % If the endpoints are not close enough, skip to the next pair of curves
                if endpoint_1_distances(response_curve_index, true_curve_index) > matching_endpoint_cutoff ||...
                    endpoint_2_distances(response_curve_index, true_curve_index) > matching_endpoint_cutoff;
                    continue
                end
                
                % Calculate the index vales for each section of the curves between the closepoints
                response_center_indices = [r1_index:r2_index]'; % We know r1 is closer to t=0 on the curve
                if t1_index < t2_index % t1 may be larger or smaller than t2
                    true_center_indices = [t1_index:t2_index]'; % increasing order
                else
                    true_center_indices = [t1_index:-1:t2_index]'; % decreasing order
                end
                
                % Copy curve segments for convenience
                r_curve_center = response_curve_points{response_curve_index}(response_center_indices, :);
                t_curve_center = true_curve_points{true_curve_index}(true_center_indices, :);
                
                % Get the area between the closepoints as an average of the Riemann sum calculated from each curve. Assumes curves are close to parallel.
                response_values = distance_to_bezier(t_curve_center, r_curve_center, true);
                response_center_area = sqrt(sum((r_curve_center(2:end, :) - r_curve_center(1:end - 1, :)).^2, 2)) .* ((response_values(1:end-1) + response_values(2:end)) / 2);
                true_values = distance_to_bezier(r_curve_center, t_curve_center, true);
                true_center_area = sqrt(sum((t_curve_center(2:end, :) - t_curve_center(1:end - 1, :)).^2, 2)) .* ((true_values(1:end-1) + true_values(2:end)) / 2);
                center_area = (response_center_area + true_center_area) / 2;
                
                % Estimate endpoint areas as triangles
                % Most endpoints will have one endpoint and closepoint as the same point
                end_close_dist_1 = sqrt(sum((end_group_1(3:4, :) - end_group_1(1:2, :)).^2, 2));
                if end_close_dist_1(1) == 0 % Triangle between closepoints and second endpoint
                    end_area_1 = triangle_area(end_group_1(2, :), end_group_1(3, :), end_group_1(4, :));
                elseif end_close_dist_1(2) == 0 % Triangle between closepoints and first endpoint
                    end_area_1 = triangle_area(end_group_1(1, :), end_group_1(3, :), end_group_1(4, :));
                else % Rare, but two triangles between closepoints and first endpoint and endpoints and second closepoint
                    end_area_1 = triangle_area(end_group_1(1, :), end_group_1(3, :), end_group_1(4, :)) +...
                                 triangle_area(end_group_1(1, :), end_group_1(2, :), end_group_1(4, :));
                end
                end_close_dist_2 = sqrt(sum((end_group_2(3:4, :) - end_group_2(1:2, :)).^2, 2));
                if end_close_dist_2(1) == 0 % Triangle between closepoints and second endpoint
                    end_area_2 = triangle_area(end_group_2(2, :), end_group_2(3, :), end_group_2(4, :));
                elseif end_close_dist_1(2) == 0 % Triangle between closepoints and first endpoint
                    end_area_2 = triangle_area(end_group_2(1, :), end_group_2(3, :), end_group_2(4, :));
                else % Rare, but two triangles between closepoints and first endpoint and endpoints and second closepoint
                    end_area_2 = triangle_area(end_group_2(1, :), end_group_2(3, :), end_group_2(4, :)) +...
                                 triangle_area(end_group_2(1, :), end_group_2(2, :), end_group_2(4, :));
                end
                
                % Sum the final area
                area_estimates(response_curve_index, true_curve_index) = center_area + end_area_1 + end_area_2;
            end
        end % end true-response curve pair loops
        
        % Record matches and extra curves
        for curve_index = 1:number_response_curves
            
            % Get best match for each curve
            [min_area, min_index] = min(area_estimates(curve_index, :));
            
            % Record a match if it is within tolerence
            if min_area <= matching_area_cutoff * true_arc_length(min_index)
                found_lines.contrast_ratios = [found_lines.contrast_ratios; info.contrast_ratio];
                found_lines.ECI = [found_lines.ECI; ts.ECI{image_index}];
                found_lines.TCI = [found_lines.TCI; ts.TCI{image_index}];
                found_lines.area = [found_lines.area; min_area];
                found_lines.true_length = [found_lines.true_length; true_arc_length(min_index)];
                found_lines.type = [found_lines.type; info.line_type];
                found_lines.size = [found_lines.size; info.line_size];
                found_true_curve_indices = [found_true_curve_indices; min_index]; % Cleared for each new image above
            
            % Otherwise, the response can't be matched and it's an extra.
            else
                extra_lines.contrast_ratios = [extra_lines.contrast_ratios; info.contrast_ratio];
                extra_lines.ECI = [extra_lines.ECI; ts.ECI{image_index}];
                extra_lines.TCI = [extra_lines.TCI; ts.TCI{image_index}];
                extra_lines.true_length = [extra_lines.true_length; true_arc_length(min_index)];
                extra_lines.type = [extra_lines.type; info.line_type];
                extra_lines.size = [extra_lines.size; info.line_size];
            end
        end 
    end % end if statement for no responses
    
    % Record missing curves
    for curve_index = 1:number_true_curves 
        if ~ismember(curve_index, found_true_curve_indices)
            missing_lines.contrast_ratios = [missing_lines.contrast_ratios; info.contrast_ratio];
            missing_lines.ECI = [missing_lines.ECI; ts.ECI{image_index}];
            missing_lines.TCI = [missing_lines.TCI; ts.TCI{image_index}];
            missing_lines.true_length = [missing_lines.true_length; true_arc_length(min_index)];
            missing_lines.type = [missing_lines.type; info.line_type];
            missing_lines.size = [missing_lines.size; info.line_size];
        end
    end
end % end actin image loop

% Copy into summary structure
ts.actin.found_points = found_points;
ts.actin.extra_points = extra_points;
ts.actin.missing_points = missing_points;

% -------Border analysis-----------

% Initalize results variables
borders = template_struct;
borders.net_area = [];
borders.abs_area = [];
borders.rmsd = [];

for image_index = border_indices;
    info = test_info.image_info{image_index};
    ground_truth = test_info.ground_truth_coords;
    
    % Measure the fractal dimension of the border
    fractal_dim = calc_fractal_dimension(ground_truth);
    
    % Get interpolated response coordinates
    interp_response = struct();
    interp_response.x = interp1(response.y{image_index}, response.x{image_index}, ground_truth(:, 2));
    
    % Get distances between true and response curves
    dist_x = interp_response.x - ground_truth(:, 1);
    delta_x = dist_x(1:end-1) - dist_x(2:end);
    delta_y = ground_truth(1:end-1, 2) - ground_truth(2:end, 2);
    
    % Calculate areas and rmsd
    net_area = sum(delta_x .* delta_y, 1);
    abs_area = sum(abs(delta_x) .* delta_y, 1);
    rmsd = sqrt(mean(dist_x.^2));
    
    % Copy information
    borders.net_area = [borders.net_area; net_area];
    borders.abs_area = [borders.abs_area; abs_area];
    borders.rmsd = [borders.rmsd; rmsd];
    borders.contrast_ratios = [borders.contrast_ratios; info.contrast_ratio];
    borders.ECI = [borders.ECI; ts.ECI{image_index}];
    borders.TCI = [borders.TCI; ts.TCI{image_index}];
end

% Copy into summary structure
ts.borders = borders;

% ----- Final cleanup ------

% Rename and save results
test_summary = ts;
if ~strcmp(output_directory, '');
    filepath = [output_directory, '/test_summary.mat'];
else
    filepath = 'test_summary.mat';
end
save(filepath, 'test_summary');
end
