function Parametric_image_tester(GUI_info)
%PARAMETRIC_IMAGE_TESTER Records human parametric interpretation of STORM
% images prepared with function make_test. 

% Set default graphs to black
colordef black

% Create figure, setting up properties
hfig = figure('Name','Parametric STORM Image Interpretation Tester', 'Toolbar', 'none', 'Menubar', 'none',...
              'Position', [25, 50, 1200, 950], 'NumberTitle', 'off', 'Visible', 'off', 'Renderer', 'zbuffer');
set(hfig, 'WindowButtonUpFcn', @clickup_callback);
          
%  --------------Initialize parameters and set options-----------

% Determine coordinate bounds
params = GUI_info.params;
min_x_bound = params.bounds(1);
min_y_bound = params.bounds(2);
max_x_bound = params.bounds(3);
max_y_bound = params.bounds(4);

% Get required constant values
pixel_size = params.STORM_pixel_size;
total_number_images = length(GUI_info.STORM_images);

% Initalize variables
current_image_type = '';
current_answer_x = [];
current_answer_y = [];
current_image_index = 0;
current_image_size = size(GUI_info.STORM_images{1});

%  ---------------Create graphs-------------------

% Create point data variables
X_graph = [];
Y_graph = [];
C_graph = [];

% Create main axes objects and handles
haxes = axes('Units', 'pixels', 'Position', [300, 50, 850, 850]);
hscatter = scatter(haxes, X_graph, Y_graph, 35, C_graph);
set(hscatter,  'PickableParts', 'none'); % Points won't cover up clicking on the axes
set(haxes, 'DataAspectRatio', [1,1,1], 'Xlim', [min_x_bound, max_x_bound], 'Ylim', [min_y_bound, max_y_bound], 'Color', 'none');

% Create the background image axes with a default all black background
hbackaxes = axes('Units', 'pixels', 'Position', get(haxes, 'Position'));
hbackimage = imagesc(zeros(current_image_size, 'uint8'));
uistack(hbackaxes,'bottom');% Move the background axes to the bottom
set(hbackaxes, 'XTickLabel', [], 'YTickLabel', [], 'XTick', [], 'YTick', [], 'Ydir', 'rev');
colormap(hbackaxes, gray);

% ----------Create buttons-----------------

% Create image heading
image_number_text = uicontrol('Parent', hfig, 'Style', 'text',...
    'Position', [25, 890, 225, 25], 'String', 'Image Number 1', 'FontSize', 16);
image_type_text = uicontrol('Parent', hfig, 'Style', 'text',...
    'Position', [25, 850, 225, 25], 'String', 'Image Type', 'FontSize', 14);

% Top of the button stack
top = 800;

% Create button for getting help
help_button = uicontrol('Parent', hfig, 'Style', 'pushbutton',...
    'Position', [25, top-25, 225, 25], 'String', 'Display Instructions',...
    'Callback', @help_button_callback);

% Create button for save and continue.
next_button = uicontrol('Parent', hfig, 'Style', 'pushbutton',...
    'Position', [25, top-115, 225, 25], 'String', 'Start Test',...
    'Callback', @next_button_callback);

% Create undo button
undo_button = uicontrol('Parent', hfig, 'Style', 'pushbutton',...
    'Position', [25, top-155, 225, 25], 'String', 'Undo Last Changes',...
    'Callback', @undo_button_callback);

% Create start over button
start_over_button = uicontrol('Parent', hfig, 'Style', 'pushbutton',...
    'Position', [25, top-195, 225, 25], 'String', 'Start Over (this image only)',...
    'Callback', @start_over_button_callback);

% Create pan button
pan_button = uicontrol('Parent', hfig, 'Style', 'pushbutton',...
    'Position', [25, top-285, 225, 25], 'String', 'Pan View',...
    'Callback', @pan_callback);

% Create zoom button
zoom_button = uicontrol('Parent', hfig, 'Style', 'pushbutton',...
    'Position', [25, top-325, 225, 25], 'String', 'Zoom View',...
    'Callback', @zoom_callback);

% Create view reset button
reset_button = uicontrol('Parent', hfig, 'Style', 'pushbutton',...
    'Position', [25, top-365, 225, 25], 'String', 'Reset View',...
    'Callback', @reset_callback);

% Create overlay toggle
overlay_toggle_button = uicontrol('Parent', hfig, 'Style', 'togglebutton',...
    'Position', [25, top-405, 225, 25], 'String', 'Toggle Response Overlay',...
    'Value', 0, 'Callback', @overlay_toggle_callback);

% -----------Start the GUI---------------

% Display GUI
set(hfig, 'Visible', 'on');

% Update graph
graph_update

%--------- Callback functions -----------

% -- Buttons --

function zoom_callback(source, callbackdata)
    zoom
end

function pan_callback(source, callbackdata)
    pan
end

function reset_callback(source, callbackdata)
    zoom out
end

function overlay_toggle_callback(source, callbackdata)
    % Get the toggle button value
    toggle_val = get(source, 'Value');
    
    % Change the stacking of the plots
    if toggle_val == 1
        uistack(hbackaxes,'bottom');% Move the background axes to the bottom
        graph_update
    elseif toggle_val == 0
        uistack(hbackaxes,'top');% Move the background axes to the top
        graph_update
    end
end

function next_button_callback(source, callbackdata)
    
    % If it's the first time clicked, change button text
    if current_image_index == 0
        set(next_button, 'String', 'Save Response and Load New Image');
    
    % If it's not the first time, we'll have a response to save
    else
        save_response
    end
    
    % Add one to the current_image_index
    current_image_index = current_image_index + 1;
    
    % If this is the last image, call the cleanup routine
    if current_image_index > total_number_images
        final_cleanup
    end
    
    % Clear the response 
    current_answer_x = [];
    current_answer_y = [];
    
    % Get the image type
    info = GUI_info.image_info{current_image_index};
    current_image_type = info.image_type;
    current_image_size = size(GUI_info.STORM_images{current_image_index});
    
    % Put new image into background image data variable
    set(hbackimage, 'CData', GUI_info.STORM_images{current_image_index});
    
    % Reassign the clickdown callback functions depending on the image type
    if strcmp(current_image_type, 'region')
        set(haxes, 'ButtonDownFcn', @region_dots_clickdown_callback);
    elseif strcmp(current_image_type, 'dots')
        set(haxes, 'ButtonDownFcn', @region_dots_clickdown_callback);
    elseif strcmp(current_image_type, 'actin')
        set(haxes, 'ButtonDownFcn', @actin_clickdown_callback);
    elseif strcmp(current_image_type, 'dots')
        set(haxes, 'ButtonDownFcn', @border_clickdown_callback);
    end
    
    % Reset zoom
    zoom out
    
    % Update graphs
    graph_update 
end

% -- Mouse stuff --

function move_point_callback(source, callbackdata, index)
    
    % Moves the current_answer_x/y point at the given index to the current mouse position
    mouse_point = get(haxes, 'CurrentPoint');
    current_answer_x(index) = mouse_point(1, 1);
    current_answer_y(index) = mouse_point(1, 2);
    graph_update
end

function click_bezier_callback(source, callbackdata, index)
    
    % See if there is an active point nearby
    clickdown_coords = callbackdata.IntersectionPoint(1:2);
    existing_point_index = find_existing_point(clickdown_coords, 0.01, true);
    
    % Pass the click along if it is near an active point
    if ~isempty(existing_point_index)
        actin_clickdown_callback(source, callbackdata)
    
    % Otherwise, work with the entire line
        else

        % Determine which button was clicked
        click_type = callbackdata.Button;

        if click_type == 1 % Left click

            % Moves the control points for the line at the indicated index to the last position
            moved_x = current_answer_x(index, :);
            moved_y = current_answer_y(index, :);
            current_answer_x(index, :) = [];
            current_answer_y(index, :) = [];
            current_answer_x = vertcat(current_answer_x, moved_x);
            current_answer_y = vertcat(current_answer_y, moved_y);

        elseif click_type == 3 % Right click

            % Delete the row for the given line
            current_answer_x(index, :) = [];
            current_answer_y(index, :) = [];
        end
    end
    graph_update
end

function create_bezier_callback(source, callbackdata, index)
    
    % Draws new straight line between clickdown point and the current cursor postion
    mouse_point = get(haxes, 'CurrentPoint');
    clickdown_x = current_answer_x(index, 1);
    clickdown_y = current_answer_y(index, 1);
    if size(current_answer_x, 2) == 2;
        current_answer_x(index, 2) = mouse_point(1, 1);
        current_answer_y(index, 2) = mouse_point(1, 2);
    elseif size(current_answer_x, 2) == 3;
        current_answer_x(index, 2) = (mouse_point(1, 1) + clickdown_x) / 2;
        current_answer_y(index, 2) = (mouse_point(1, 2) + clickdown_y) / 2;
        current_answer_x(index, 3) = mouse_point(1, 1);
        current_answer_y(index, 3) = mouse_point(1, 2);
    elseif size(current_answer_x, 2) == 3;
        current_answer_x(index, 2) = (mouse_point(1, 1) + 2 * clickdown_x) / 3;
        current_answer_y(index, 2) = (mouse_point(1, 2) + 2 * clickdown_y) / 3;
        current_answer_x(index, 3) = (2 * mouse_point(1, 1) + clickdown_x) / 3;
        current_answer_y(index, 3) = (2 * mouse_point(1, 2) + clickdown_y) / 3;
        current_answer_x(index, 4) = mouse_point(1, 1);
        current_answer_y(index, 4) = mouse_point(1, 2);
    end
    graph_update
end

function clickup_callback(source, callbackdata)

    % Stop using the motion function
    set(hfig,'WindowButtonMotionFcn', '');
    
    % Set a new undo savedpoint
    set_undo_savepoint;
end

function actin_clickdown_callback(source, callbackdata)
    
    % Get coordinates and type of click
    clickdown_coords = callbackdata.IntersectionPoint(1:2);
    click_type = callbackdata.Button;
    
    % Is the click next to an existing point?
    existing_point_index = find_existing_point(clickdown_coords, 0.01, true);
           
    % If a point was found and it's a left click, we move the point with the pointer  
    if ~isempty(existing_point_index) && click_type == 1;
        
        % Set the window motion function to move the point around with the pointer
        set(hfig,'WindowButtonMotionFcn', {@move_point_callback, existing_point_index});
        
    % If no point is found and it's a left click, add a new line to the end of the answer
    elseif isempty(existing_point_index) && click_type == 1
        
        % We need to determine what type of line we're adding
        image_info = GUI_info.image_info{current_image_index};
        if strcmp(image_info.line_type, 'line_segment')
            number_cp = 2;
        elseif strcmp(image_info.line_type, 'quadratic')
            number_cp = 3;
        elseif strcmp(image_info.line_type, 'cubic')
            number_cp = 4;
        end
        
        % Add a new row at bottom of matrix
        new_index = size(current_answer_x, 1) + 1;
        current_answer_x(new_index, 1:number_cp) = clickdown_coords(1);
        current_answer_y(new_index, 1:number_cp) = clickdown_coords(2);
        
        % Set the window motion function to move the last point around with the pointer
        set(hfig,'WindowButtonMotionFcn', {@create_bezier_callback, new_index});
    end
    
    % Update the graph
    graph_update
end

function region_dots_clickdown_callback(source, callbackdata)
    
    % Get coordinates and type of click
    clickdown_coords = callbackdata.IntersectionPoint(1:2);
    click_type = callbackdata.Button;
    
    % Is the click next to an existing point?
    existing_point_index = find_existing_point(clickdown_coords, 0.01, false);
    
    % If a point was found and it's a right click, we delete the point  
    if ~isempty(existing_point_index) && click_type == 3;     
        current_answer_x(existing_point_index) = [];
        current_answer_y(existing_point_index) = [];
        
    % If a point was found and it's a left click, we move the point with the pointer  
    elseif ~isempty(existing_point_index) && click_type == 1;
        
        % Set the window motion function to move the point around with the pointer
        set(hfig,'WindowButtonMotionFcn', {@move_point_callback, existing_point_index});
        
    % If no point is found and it's a left click, add a new point to the end of the answer
    elseif isempty(existing_point_index) && click_type == 1
        current_answer_x(end+1, 1) = clickdown_coords(1);
        current_answer_y(end+1, 1) = clickdown_coords(2);
    end
    
    % Update the graph
    graph_update
end

% ----------Other Functions--------------

function graph_update 
    
    % Get x, y, and color data for graphing    
    [Xp, Yp, Cp] = get_point_data; % Vectors of x, y, area, and and nx3 matrix of RGB color values
    [Xl, Yl, Cl, Pl] = get_line_data; % Columns are vectors of x and y points defining a line, multiple lines in multiple columns. Cl is a nx3 matrix of RGB color values. Pl is the index value that a line pick points to. 
    
    % redraw scatterplot
    set(hscatter, 'Xdata', Xp)
    set(hscatter, 'Ydata', Yp)
    set(hscatter, 'CData', Cp)
    
    % Delete all existing lines
    delete(findobj(haxes, 'type', 'line'));
    
    % Add a line plot for each particle track onto the graph
    for line_index = 1:size(Xl, 1)
        pick_index = Pl(line_index);
        
        % Plot each line with specified color and picking properties
        if pick_index == 0
            line(Xl{line_index}, Yl{line_index}, 'Parent', haxes, 'Color', Cl(line_index, :), 'PickableParts', 'none'); % Not pickable
        else
            line(Xl{line_index}, Yl{line_index}, 'Parent', haxes, 'Color', Cl(line_index, :), 'ButtonDownFcn', {@click_bezier_callback, pick_index}); % A pick sends the corrosponding index to click_bezier_callback
        end
    end
        
    % Update background zoom/pan
    main_xlim = get(haxes, 'Xlim');
    main_ylim = get(haxes, 'Ylim');
    set(hbackaxes, 'Xlim', main_xlim ./ pixel_size, 'Ylim', current_image_size(1) - ([main_ylim(2), main_ylim(1)] ./ pixel_size));
end

function [X, Y, C] = get_point_data
    % Returns the points and corrosponding colors of the points
    
    % Drawing vertices of the region polygon
    if strcmp(current_image_type, 'region')
        
        % Orange for the first point, Yellow for the last, Red for all others
        X = current_answer_x(:);
        Y = current_answer_y(:);
        C = repmat([1, 0, 0], size(X, 1), 1);
        if size(X, 1) > 0
            C(1, :) = [1, .5, 0];
            C(end, :) = [1, 1, 0];
        end
            
    % Marking the centers of the dots
    elseif strcmp(current_image_type, 'dots')
            
        % Yellow for the last, Red for all others
        X = current_answer_x(:);
        Y = current_answer_y(:);
        C = repmat([1, 0, 0], size(X, 1), 1);
        if size(X, 1) > 0
            C(end, :) = [1, 1, 0];        
        end
        
    % Marking the control points    
    elseif strcmp(current_image_type, 'actin')
        
        % Orange for the first point, Yellow for the last, Red for all others
        % All but the last set of control points are greyish
        [num_lines, num_cp] = size(current_answer_x);
        X = current_answer_x(:);
        Y = current_answer_y(:);
        C = repmat([.5, 0, 0], size(X, 1), 1);
        if size(X, 1) < 0
            C(1:num_lines:end, :) = [.5, .25, 0];
            C(num_cp:num_lines:end, :) = [.5, .5, 0];
            C(end-num_cp+1, :) = [1, 0, 0];
            C(end-num_cp+2:end, :) = [1, .5, 0];
            C(end, :) = [1, 1, 0];
        end
        
    % No points graphed for borders    
    elseif strcmp(current_image_type, 'border')
        
        % Return empty matrices
        X = [];
        Y = [];
        C = zeros(0, 3);
    
    % Return empty matrices for any other image type
    else
        X = [];
        Y = [];
        C = zeros(0, 3);
    end
    
    
end

function [X, Y, C, P] = get_line_data
    % Returns cell arrays of points that define lines on the graph and the corrosponding colors of the lines
    % A line is a column vector of x of y values in the cell, colors are a matrix with columns of RGB triple vectors
    % Matrix P is the row index that a pickable line will return, an index of zero makes the line not pickable.  
    
    % Drawing sides of the region polygon
    if strcmp(current_image_type, 'region')
        
        % No lines for 0 or 1 points
        if size(current_answer_x, 1) < 2
            X = {};
            Y = {};
            C = zeros(0, 3);
            P = zeros(0, 3);
        else % Orange for the line between the first and last point, Red for all others
            X_mat = vertcat(current_answer_x', [current_answer_x(2:end)', current_answer_x(1)]);
            Y_mat = vertcat(current_answer_y', [current_answer_y(2:end)', current_answer_y(1)]);
            X = cell(size(X_mat, 2), 1);
            Y = cell(size(Y_mat, 2), 1);
            for index = 1:size(X_mat, 2)
                X{index} = X_mat(:, index);
                Y{index} = Y_mat(:, index);
            end
            C = repmat([1, 0, 0], size(X_mat, 2), 1);
            C(end, :) = [1, .5, 0];
        end
        P = zeros(size(X));
            
    % No lines for the dots images
    elseif strcmp(current_image_type, 'dots')
            
       % Return empty cell arrays and matrices
        X = {};
        Y = {};
        C = zeros(0, 3);
        P = zeros(0, 3);
        
    % Drawing control net and bezier lines 
    elseif strcmp(current_image_type, 'actin')
        
        % Red for all control net lines
        % All but the last set of control points are greyish
        if size(current_answer_x, 2) == 2
            X_mat = current_answer_x';
            Y_mat = current_answer_y';
        elseif size(current_answer_x, 2) == 3
            X_mat = horzcat(current_answer_x(:, 1:2)', current_answer_x(:, 2:3)');
            Y_mat = horzcat(current_answer_y(:, 1:2)', current_answer_y(:, 2:3)');
        elseif size(current_answer_x, 2) == 4
            X_mat = horzcat(current_answer_x(:, 1:2)', current_answer_x(:, 2:3)', current_answer_x(:, 3:4)');
            Y_mat = horzcat(current_answer_y(:, 1:2)', current_answer_y(:, 2:3)', current_answer_y(:, 3:4)');
        else % no answer yet
            X_mat = [];
            Y_mat = [];
        end  
        X_cn = cell(size(X_mat, 2), 1);
        Y_cn = cell(size(Y_mat, 2), 1);
        for index = 1:size(X_mat, 2)
            X_cn{index} = X_mat(:, index);
            Y_cn{index} = Y_mat(:, index);
        end
        C_cn = repmat([.5, 0, 0], size(X_mat, 2), 1);
        if size(C_cn, 1) > 0;
            C_cn(end, :) = [1, 0, 0];
        end
        P_cn = zeros(size(X_cn));
        
        % Calculate and graph the bezier lines
        % Beziers are orange; all but the last one is greyish        
        X_bz = cell(size(current_answer_x, 1), 1);
        Y_bz = cell(size(current_answer_y, 1), 1);
        P_bz = zeros(size(current_answer_x, 1), 1);
        for index = 1:size(X_bz, 1)
            line_points = calc_bezier_line([current_answer_x(index, :)', current_answer_y(index, :)'], 1000);
            X_bz{index} = line_points(:, 1);
            Y_bz{index} = line_points(:, 2);
            P_bz(index) = index;
        end
        C_bz = repmat([.5, .25, 0], size(current_answer_x, 1), 1);
        if size(C_bz, 1) > 0;
            C_bz(end, :) = [1, .5, 0];
        end
        
        % Add control network and bezier cells together
        X = vertcat(X_cn, X_bz);
        Y = vertcat(Y_cn, Y_bz);
        C = vertcat(C_cn, C_bz);
        P = vertcat(P_cn, P_bz);
        
    % Graph the border line   
    elseif strcmp(current_image_type, 'border')
        
        % Simply return the current answer
        X = {current_answer_x};
        Y = {current_answer_y};
        C = [1, 0, 0];
        P = [0];
    
    % Return empty cells for any other image type
    else
        X = {};
        Y = {};
        C = zeros(0, 3);
        P = zeros(0, 3);
    end
end

function final_cleanup
    
    % Cleanup stuff - placeholder
   
end

function save_response
    
    % Save stuff - placeholder
   
end

function set_undo_savepoint
    
    % Save undo - placeholder
   
end

function [existing_point_index] = find_existing_point(click_coords, relative_pick_radius, last_row_only)
%   Function to determine if a click point is next to and existing answer point, and if so, which one.
%   Inputs: 
%       click_coords: click point given as [x, y] coordinates
%       relative_pick_radius: the radius of the valid picking region around
%           each point, given as a fraction of the axes width. Can be a 
%           scalar or a matrix [x_radius, y_radius].
%       last_row_only: boolean, when true we only compare the click the the
%           last row of current_answer_x/y points.
%   Output:
%       existing_point_index: index value for the nearest point within the
%           picking radius in the current_answer_x/y matrices. If no points
%           are within the picking radius, an empty matrix is returned.
%           Returned value is a linear index.
    
    % Don't bother looking if there are no points yet
    if isempty(current_answer_x) && isempty(current_answer_y)
        existing_point_index = [];
    
    else
        % Deal with both types of pick radii
        rel_radius_x = relative_pick_radius(1);
        if length(relative_pick_radius) == 2;
            rel_radius_y = relative_pick_radius(2);
        else
            rel_radius_y = relative_pick_radius(1);
        end

        % Figure out the point picking radius
        current_xlim = get(haxes, 'Xlim');
        current_ylim = get(haxes, 'Ylim');
        x_range = current_xlim(2) - current_xlim(1);
        y_range = current_ylim(2) - current_ylim(1);
        pick_radius_x = rel_radius_x * x_range; % 1% of the image width, no matter what zoom level we are at
        pick_radius_y = rel_radius_y * y_range;

        % Get the click to point distances in units of picking radii
        if last_row_only
            xdist = (click_coords(1) - current_answer_x(end, :)) / pick_radius_x;
            ydist = (click_coords(2) - current_answer_y(end, :)) / pick_radius_y;
        else
            xdist = (click_coords(1) - current_answer_x(:)) / pick_radius_x;
            ydist = (click_coords(2) - current_answer_y(:)) / pick_radius_y;
        end
        dist = sqrt(xdist.^2 + ydist.^2);

        % If any point is close enough, make sure to pick the one closest to the click 
        if any(dist <= 1)
            [~, existing_point_index] = min(dist);
        else
            existing_point_index = [];
        end
    end
end
    
end