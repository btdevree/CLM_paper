% Checking accurcy of cubic bezier approximation algorithms

% Generate random beziers on a unit 1 square grid
number_beziers = 1000;
cp_x = rand(number_beziers, 4);
cp_y = rand(number_beziers, 4);

% Get the approximated and actual lengths for each curve
approx_arc_length = zeros(number_beziers, 1);
true_arc_length = zeros(number_beziers, 1);
for curve_index = 1:number_beziers
    
    % Construct the control_points matrix
    control_points = [cp_x(curve_index, :)', cp_y(curve_index, :)'];
    
    % Approximate the arc length
    chord_length = sqrt(sum((control_points(end, :) - control_points(1, :)).^2 , 2));
    control_net_length = sum(sqrt(sum((control_points(2:end, :) - control_points(1:end-1, :)).^2 , 2)), 1);
    approx_arc_length(curve_index) = (chord_length + control_net_length) / 2; 

    % Get points along the bezier curve 
    curve_points = calc_bezier_line(control_points, 100);
    true_arc_length(curve_index) = sum(sqrt(sum((curve_points(2:end, :) - curve_points(1:end-1, :)).^2 , 2)), 1);
end

% Get the approximated arc lengths for each curve
number_points = [2, 5, 10, 20, 50, 100, 200, 500, 1000, 2000, 5000];
arc_lengths = zeros(number_beziers, length(number_points));
for curve_index = 1:number_beziers
    
    % Construct the control_points matrix
    control_points = [cp_x(curve_index, :)', cp_y(curve_index, :)'];
    
    % Loop through each point number
    for number_index = 1:length(number_points)

        % Get points along the bezier curve 
        curve_points = calc_bezier_line(control_points, number_points(number_index));
        arc_lengths(curve_index, number_index) = sum(sqrt(sum((curve_points(2:end, :) - curve_points(1:end-1, :)).^2 , 2)), 1);
    end
end
relative_arc_lengths = arc_lengths - repmat(arc_lengths(:, end), 1, size(arc_lengths, 2));

% Plot
figure
scatter(1:1000, (approx_arc_length - true_arc_length)./true_arc_length)

figure
loglog(number_points, relative_arc_lengths)

% Notes:
% First method - not awesome, 0 to 100% overestimate of error, typical
% around 30%.
% Second method - clear log-log dependence, worst cases at 1% error around
% 20 points, 0.1% around 50 points.