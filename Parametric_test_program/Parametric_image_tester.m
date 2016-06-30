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

%  ---------------Create graphs-------------------

% Create point data variables
X_graph = [];
Y_graph = [];
A_graph = [];
C_graph = [];

% Create main axes object and handle
haxes = axes('Units', 'pixels', 'Position', [300, 50, 850, 850]);

% Display scatterplot in the axes and get a handle to the graph
hscatter = scatter(haxes, X_graph, Y_graph, A_graph, C_graph);
set(haxes, 'DataAspectRatio', [1,1,1], 'Xlim', [min_x_bound, max_x_bound], 'Ylim', [min_y_bound, max_y_bound], 'Color', 'none');

% Create the background image axes with a default all black background
hbackaxes = axes('Units', 'pixels', 'Position', get(haxes, 'Position'));
hbackimage = imagesc(zeros(size(GUI_info.STORM_images{1}), 'uint8'), [0, 1]);
uistack(hbackaxes,'bottom');% Move the background axes to the bottom
set(hbackaxes, 'XTickLabel', [], 'YTickLabel', [], 'XTick', [], 'YTick', [], 'Ydir', 'rev');
colormap(hbackaxes, gray);

% ----------Create buttons-----------------

% Create image heading
image_number_text = uicontrol('Parent', hfig, 'Style', 'text',...
    'Position', [25, 890, 225, 25], 'String', 'Practice Image 1', 'FontSize', 16);
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
    'Position', [25, top-115, 225, 25], 'String', 'Save and Continue to Next Image',...
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

% Load first image




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
    
    % Change the data and color limits to the appropreate values
    if toggle_val == 1
        image_data = image_filename_list{current_movie_index};
        set(hbackimage, 'CData', image_data);
        set(hbackaxes, 'CLim', [min(min(image_data)), max(max(image_data))]); 
    elseif toggle_val == 0
        set(hbackimage, 'CData', zeros(size(image_filename_list{current_movie_index})));
        set(hbackaxes, 'CLim', [0, 1]); 
    end
    
    % Update the zoom levels 
    main_xlim = get(haxes, 'Xlim');
    main_ylim = get(haxes, 'Ylim');
    set(hbackaxes, 'Xlim', main_xlim, 'Ylim', size(image_filename_list{current_movie_index}, 1) - [main_ylim(2), main_ylim(1)]);
end

% ----------Other Functions--------------

function graph_update
    
    % Get x, y, and color data for graphing    
    [Xp,Yp,Ap,Cp] = get_point_data; % Vectors of x, y, area, and and nx3 matrix of RGB color values
    [Xl,Yl,Cl] = get_line_data; % Columns are vectors of x and y points defining a line, multiple lines in multiple columns. Cl is a nx3 matrix of RGB color values.  
    
    % redraw scatterplot
    set(hscatter, 'Xdata', Xp)
    set(hscatter, 'Ydata', Yp)
    set(hscatter, 'CData', Cp)
    set(hscatter, 'SizeData', Ap)
    
    % Delete all existing lines and rectangles
    delete(findobj(haxes, 'type', 'line'));
    delete(findobj(haxes, 'type', 'rectangle'));
    
    % Add a line plot for each particle track onto the graph
    for line_index = 1:size(Xl, 2)
       
        % Determine what color the line should be
        if particle_list{current_movie_index}(track_index) == current_particle_ID
            plot_color = [0, 1, 0];
        else
            plot_color = [.25, .35, 1];
        end
        
        % Collect xy points for the plot and plot them as a line
        plot_XY_data = particle_xy_list{current_movie_index}{particle_list{current_movie_index}(track_index)};
        line(plot_XY_data(:,1), plot_XY_data(:,2), 'Parent', haxes, 'Color', plot_color);
    
        % Plot a search circle at the current frame
        rect_x = plot_XY_data(current_frame_index,1)-search_radius;
        rect_y = plot_XY_data(current_frame_index,2)-search_radius;
        rect_diameter = search_radius*2;
        rectangle('Position', [rect_x, rect_y, rect_diameter, rect_diameter],...
            'Curvature', [1,1], 'EdgeColor', plot_color, 'Parent', haxes)
        
        % Plot a particle circle at the current frame
        rect_x = plot_XY_data(current_frame_index,1)-particle_radius;
        rect_y = plot_XY_data(current_frame_index,2)-particle_radius;
        rect_diameter = particle_radius*2;
        rectangle('Position', [rect_x, rect_y, rect_diameter, rect_diameter],...
            'Curvature', [1,1], 'EdgeColor', plot_color, 'Parent', haxes)
        
        % Update background zoom/pan
        main_xlim = get(haxes, 'Xlim');
        main_ylim = get(haxes, 'Ylim');
        set(hbackaxes, 'Xlim', main_xlim, 'Ylim', size(image_filename_list, 1) - [main_ylim(2), main_ylim(1)]);
    end
end
end