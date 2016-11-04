% Script to create example images for parametric test figures

% Get files
untar('test_archive.tar.gz');
load('test_archive.mat'); % Loads test_info

% Prepare image names
image_names = {'region_cr1-1.5_example.png'; 'region_cr1-2_example.png'; 'region_cr1-5_example.png';...
               'dots_50nm_example.png'; 'dots_100nm_exmaple.png'; 'dots_200nm_example.png';...
               'lines_9nm_example.png'; 'lines_26nm_example.png';...
               'border_fd1.05_example.png'; 'border_fd1.2_example.png'; 'border_fd1.5_example.png'};

% Valid image indices will be used as given, anything with a NaN will be randomly chosen fromthe available appropriate choices 
% image_indices = [NaN; NaN; NaN;...
%                  NaN; NaN; NaN;...
%                  NaN; NaN;...
%                  NaN; NaN; NaN];           
image_indices = [NaN; NaN; NaN;...
                 NaN; NaN; NaN;...
                 NaN; NaN;...
                 NaN; NaN; NaN];

% Find appropreate region images;
for index = 1:3;
   if isnan(image_indices(index))
       
       % Set the variables for matching
       target_image_type = 'region';
       if index == 1
           target_image_param = 0.5;
       elseif index == 2
           target_image_param = 1;
       elseif index == 3
           target_image_param = 4;
       end
       
       % Search through images
       possible_indices = [];
       for target_index = 1:size(test_info.image_info, 1)
           if strcmp(test_info.image_info{target_index}.image_type, target_image_type) && test_info.image_info{target_index}.contrast_ratio == target_image_param
                possible_indices(end + 1) = target_index;
           end
       end
       
       % Randomly choose an image from the possible ones
       number_choices = length(possible_indices(:));
       if number_choices > 0
           image_indices(index) = possible_indices(randi(number_choices));
       else
           warning([' No appropriate images for creating ', image_names{index}, ', please look in another test file.'])
           image_indices(index) = NaN; % will be skipped by the image making section
       end
   end    
end

% Find appropreate dots images;
for index = 4:6;
   if isnan(image_indices(index))
       
       % Set the variables for matching
       target_image_type = 'dots';
       if index == 4
           target_image_param = 50;
       elseif index == 5
           target_image_param = 100;
       elseif index == 6
           target_image_param = 200;
       end
       
       % Search through images
       possible_indices = [];
       for target_index = 1:size(test_info.image_info, 1)
           if strcmp(test_info.image_info{target_index}.image_type, target_image_type) && test_info.image_info{target_index}.dot_sizes == target_image_param
               possible_indices(end + 1) = target_index;
           end
       end
       
      % Randomly choose an image from the possible ones
       number_choices = length(possible_indices(:));
       if number_choices > 0
           image_indices(index) = possible_indices(randi(number_choices));
       else
           warning([' No appropriate images for creating ', image_names{index}, ', please look in another test file.'])
           image_indices(index) = NaN; % will be skipped by the image making section
       end
   end    
end

% Find appropreate lines images;
for index = 7:8;
   if isnan(image_indices(index))
       
       % Set the variables for matching
       target_image_type = 'actin';
       if index == 7
           target_image_param = 9;
       elseif index == 8
           target_image_param = 26;
       end
       
       % Search through images
       possible_indices = [];
       for target_index = 1:size(test_info.image_info, 1)
           if strcmp(test_info.image_info{target_index}.image_type, target_image_type) && test_info.image_info{target_index}.line_width == target_image_param
               possible_indices(end + 1) = target_index;
           end
       end
       
        % Randomly choose an image from the possible ones
       number_choices = length(possible_indices(:));
       if number_choices > 0
           image_indices(index) = possible_indices(randi(number_choices));
       else
           warning([' No appropriate images for creating ', image_names{index}, ', please look in another test file.'])
           image_indices(index) = NaN; % will be skipped by the image making section
       end
   end    
end

% Find appropreate border images;
for index = 9:11;
   if isnan(image_indices(index))
       
       % Set the variables for matching
       target_image_type = 'border';
       if index == 9
           target_image_param = .45;
       elseif index == 10
           target_image_param = .6;
       elseif index == 11
           target_image_param = .75;
       end
       
       % Search through images
       possible_indices = [];
       for target_index = 1:size(test_info.image_info, 1)
           if strcmp(test_info.image_info{target_index}.image_type, target_image_type) && test_info.image_info{target_index}.roughness == target_image_param
               possible_indices(end + 1) = target_index;
           end
       end
       
      % Randomly choose an image from the possible ones
       number_choices = length(possible_indices(:));
       if number_choices > 0
           image_indices(index) = possible_indices(randi(number_choices));
       else
           warning([' No appropriate images for creating ', image_names{index}, ', please look in another test file.'])
           image_indices(index) = NaN; % will be skipped by the image making section
       end
   end    
end

% Now image_indices has a selection of appropriate images for each example
% Graph and print each image
for example_index = 1:size(image_indices, 1)
    summary_index = image_indices(example_index);
    
    % Make the image 
    if ~isnan(summary_index)
        make_nonparameteric_plot(image_names{example_index}, test_info.STORM_images{summary_index}, 21, true, true);
    end
end

% Save the index numbers into a text file
fileID = fopen('image_indices.txt','w');
for index = 1:size(image_names, 1)
    fprintf(fileID,'Image named %s made with image number %d\n', image_names{index}, image_indices(index));
end
fclose(fileID);

% Delete the un-zipped archive file
delete('test_archive.mat');