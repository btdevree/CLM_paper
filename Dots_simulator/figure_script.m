% Set directory information
figure_path_parts = strsplit(pwd, 'CLM_paper');
figure_path = [figure_path_parts{1}, 'CLM_figures_and_data/'];
data_name_no = 'stdev_vs_num_dots_no_over.mat';
data_name_low = 'stdev_vs_num_dots_low_over.mat';
data_name_high = 'stdev_vs_num_dots_high_over.mat';
figure_prefix = 'stdev_vs_num_dots';

% Note - hard-coded for 7 nm pixel distance

% Load data from no overcounting
load([figure_path, data_name_no])

% Create stdev vs number dots curve of the radial averages
stdev_no_over_cells = cell(size(mean_vector_stack_cells));
stdev_no_over_0nm = zeros(size(mean_vector_stack_cells), 1);
stdev_no_over_105nm = zeros(size(mean_vector_stack_cells), 1);
stdev_no_over_210nm = zeros(size(mean_vector_stack_cells), 1);
for num_dots_index = 1:length(mean_vector_stack_cells)
    stdev_vs_dist_vector = std(mean_vector_stack_cells, 0, 2);
    stdev_no_over_cells{num_dots_index} = stdev_vs_dist_vector;
    stdev_no_over_0nm(num_dots_index, 1) = stdev_vs_dist_vector(1);
    stdev_no_over_105nm(num_dots_index, 1) = stdev_vs_dist_vector(16);
    stdev_no_over_210nm(num_dots_index, 1) = stdev_vs_dist_vector(31);
end

% Load data from no overcounting
load([figure_path, data_name_low])

% Create stdev vs number dots curve of the radial averages
stdev_low_over_cells = cell(size(mean_vector_stack_cells));
stdev_low_over_0nm = zeros(size(mean_vector_stack_cells), 1);
stdev_low_over_105nm = zeros(size(mean_vector_stack_cells), 1);
stdev_low_over_210nm = zeros(size(mean_vector_stack_cells), 1);
for num_dots_index = 1:length(mean_vector_stack_cells)
    stdev_vs_dist_vector = std(mean_vector_stack_cells, 0, 2);
    stdev_low_over_cells{num_dots_index} = stdev_vs_dist_vector;
    stdev_low_over_0nm(num_dots_index, 1) = stdev_vs_dist_vector(1);
    stdev_low_over_105nm(num_dots_index, 1) = stdev_vs_dist_vector(16);
    stdev_low_over_210nm(num_dots_index, 1) = stdev_vs_dist_vector(31);
end

% Load data from no overcounting
load([figure_path, data_name_low])

% Create stdev vs number dots curve of the radial averages
stdev_high_over_cells = cell(size(mean_vector_stack_cells));
stdev_high_over_0nm = zeros(size(mean_vector_stack_cells), 1);
stdev_high_over_105nm = zeros(size(mean_vector_stack_cells), 1);
stdev_high_over_210nm = zeros(size(mean_vector_stack_cells), 1);
for num_dots_index = 1:length(mean_vector_stack_cells)
    stdev_vs_dist_vector = std(mean_vector_stack_cells, 0, 2);
    stdev_high_over_cells{num_dots_index} = stdev_vs_dist_vector;
    stdev_high_over_0nm(num_dots_index, 1) = stdev_vs_dist_vector(1);
    stdev_high_over_105nm(num_dots_index, 1) = stdev_vs_dist_vector(16);
    stdev_high_over_210nm(num_dots_index, 1) = stdev_vs_dist_vector(31);
end

% Create figure 0nm
figure
hold on
plot(number_dots_vector, stdev_no_over_0nm, k);
plot(number_dots_vector, stdev_low_over_0nm, b);
plot(number_dots_vector, stdev_high_over_0nm, g);
xlim([10,500]);
ylim([0, 5]);
xlabel('Number of dots per radial average');
ylabel('Standard deviation of radial average value');
title('0 nm radial average');
hold off
figure_name = [figure_path, figure_prefix, '0nm.png'];
print(gcf, figure_name, '-dpng');
delete(gcf);

% Create figure 105nm
figure
hold on
plot(number_dots_vector, stdev_no_over_105nm, k);
plot(number_dots_vector, stdev_low_over_105nm, b);
plot(number_dots_vector, stdev_high_over_105nm, g);
xlim([10,500]);
ylim([0, 5]);
xlabel('Number of dots per radial average');
ylabel('Standard deviation of radial average value');
title('105 nm radial average');
hold off
figure_name = [figure_path, figure_prefix, '105nm.png'];
print(gcf, figure_name, '-dpng');

% Create figure 205nm
figure
hold on
plot(number_dots_vector, stdev_no_over_205nm, k);
plot(number_dots_vector, stdev_low_over_205nm, b);
plot(number_dots_vector, stdev_high_over_205nm, g);
xlim([10,500]);
ylim([0, 5]);
xlabel('Number of dots per radial average');
ylabel('Standard deviation of radial average value');
title('205 nm radial average');
hold off
figure_name = [figure_path, figure_prefix, '205nm.png'];
print(gcf, figure_name, '-dpng');


