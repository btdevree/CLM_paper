function Parametric_image_tester(GUI_info)
%PARAMETRIC_IMAGE_TESTER Records human parametric interpretation of STORM
% images prepared with function make_test. 

% Set default graphs to black
colordef black

% Create figure, setting up properties
hfig = figure('Name','Parametric STORM Image Interpretation Tester', 'Toolbar', 'none', 'Menubar', 'none',...
              'Position', [25, 50, 1200, 950], 'NumberTitle', 'off', 'Visible', 'off', 'Renderer', 'zbuffer');

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

% Create main axes object and handle
haxes = axes('Units', 'pixels', 'Position', [300, 50, 850, 850]);

% Display scatterplot in the axes and get a handle to the graph
hscatter = scatter(haxes, X_graph, Y_graph, 5, C_graph);
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

    % Reset zoom
    zoom out
    
    % Update graphs
    graph_update 
end


% ----------Other Functions--------------

function graph_update 
    
    % Get x, y, and color data for graphing    
    [Xp, Yp, Cp] = get_point_data; % Vectors of x, y, area, and and nx3 matrix of RGB color values
    [Xl, Yl, Cl] = get_line_data; % Columns are vectors of x and y points defining a line, multiple lines in multiple columns. Cl is a nx3 matrix of RGB color values.  
    
    % redraw scatterplot
    set(hscatter, 'Xdata', Xp)
    set(hscatter, 'Ydata', Yp)
    set(hscatter, 'CData', Cp)
    
    % Delete all existing lines
    delete(findobj(haxes, 'type', 'line'));
    
    % Add a line plot for each particle track onto the graph
    for line_index = 1:size(Xl, 1)
       
        % Plot each line with specified color
        line(Xl{line_index}, Yl{line_index}, 'Parent', haxes, 'Color', Cl(line_index, :));
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

function [X, Y, C] = get_line_data
    % Returns cell arrays of points that define lines on the graph and the corrosponding colors of the lines
    % A line is a column vector, multiple lines are in multiple cells, colors are columns of RGB triple vectors
    
    % Drawing sides of the region polygon
    if strcmp(current_image_type, 'region')
        
        % No lines for 0 or 1 points
        if size(current_answer_x, 1) < 2
            X = {};
            Y = {};
            C = zeros(0, 3);
        else % Orange for the line between the first and last point, Red for all others
            X_mat = vercat(current_answer_x', [current_answer_x(2:end)', current_answer_x(1)]);
            Y_mat = vercat(current_answer_y', [current_answer_y(2:end)', current_answer_y(1)]);
            X = cell(size(X_mat, 2), 1);
            Y = cell(size(Y_mat, 2), 1);
            for index = 1:size(X_mat, 2)
                X{index} = X_mat(:, index);
                Y{index} = Y_mat(:, index);
            end
            C = repmat([1, 0, 0], size(X_mat, 2), 1);
            C(end, :) = [1, .5, 0];
        end
            
    % No lines for the dots images
    elseif strcmp(current_image_type, 'dots')
            
       % Return empty cell arrays and matrices
        X = {};
        Y = {};
        C = zeros(0, 3);
        
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
        
        % Calculate and graph the bezier lines
        % Beziers are orange; all but the last one is greyish        
        X_bz = cell(size(current_answer_x, 1), 1);
        Y_bz = cell(size(current_answer_y, 1), 1);
        for index = 1:size(X_bz, 1)
            line_points = calc_bezier_line([current_answer_x(index, :)', current_answer_y(index, :)'], 1000);
            X_bz{index} = line_points(:, 1);
            Y_bz{index} = line_points(:, 2);
        end
        C_bz = repmat([.5, .25, 0], size(current_answer_x, 1), 1);
        if size(C_bz, 1) > 0;
            C_bz(end, :) = [1, .5, 0];
        end
        
        % Add control network and bezier cells together
        X = vertcat(X_cn, X_bz);
        Y = vertcat(Y_cn, Y_bz);
        C = vertcat(C_cn, C_bz);        
        
    % Graph the border line   
    elseif strcmp(current_image_type, 'border')
        
        % Simply return the current answer
        X = {current_answer_x};
        Y = {current_answer_y};
        C = [1, 0, 0];
    
    % Return empty cells for any other image type
    else
        X = {};
        Y = {};
        C = zeros(0, 3);
    end
end

function final_cleanup
    
    % Cleanup stuff - placeholder
   
end

function save_response
    
    % Save stuff - placeholder
   
end

end