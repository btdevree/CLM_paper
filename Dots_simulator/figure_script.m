% Set directory information
figure_path_parts = strsplit(pwd, 'CLM_paper');
figure_path = [figure_path_parts{1}, 'CLM_figures_and_data/'];
data_name_no = 'new_stdev_vs_num_dots_no_over.mat';
data_name_low = 'new_stdev_vs_num_dots_low_over.mat';
data_name_high = 'new_stdev_vs_num_dots_high_over.mat';
figure_prefix = 'new_stdev_vs_num_dots_';
figure_prefix_2 = 'new_dist_vs_corr_';

% Note - hard-coded for 7 nm pixel distance

% Load data from no overcounting
load([figure_path, data_name_no])

% Create stdev vs number dots curve of the radial averages
mean_no_over_cells = cell(size(mean_vector_stack_cells));
stdev_no_over_cells = cell(size(mean_vector_stack_cells));
stdev_no_over_0nm = zeros(size(mean_vector_stack_cells, 2), 1);
stdev_no_over_105nm = zeros(size(mean_vector_stack_cells, 2), 1);
stdev_no_over_210nm = zeros(size(mean_vector_stack_cells, 2), 1);
mean_no_over_boot_cells = cell(size(mean_bootstrap_vector_stack_cells));
stdev_no_over_boot_cells = cell(size(mean_bootstrap_vector_stack_cells));
stdev_no_over_boot_0nm = zeros(size(mean_bootstrap_vector_stack_cells, 2), size(mean_bootstrap_vector_stack_cells, 1));
stdev_no_over_boot_105nm = zeros(size(mean_vector_stack_cells, 2), size(mean_bootstrap_vector_stack_cells, 1));
stdev_no_over_boot_210nm = zeros(size(mean_vector_stack_cells, 2), size(mean_bootstrap_vector_stack_cells, 1));
for num_dots_index = 1:size(mean_vector_stack_cells, 2)
    mean_vs_dist_vector = mean(mean_vector_stack_cells{num_dots_index}, 2);
    mean_no_over_cells{num_dots_index} = mean_vs_dist_vector;
    stdev_vs_dist_vector = std(mean_vector_stack_cells{num_dots_index}, 0, 2);
    stdev_no_over_cells{num_dots_index} = stdev_vs_dist_vector;
    stdev_no_over_0nm(num_dots_index, 1) = stdev_vs_dist_vector(1);
    stdev_no_over_105nm(num_dots_index, 1) = stdev_vs_dist_vector(16);
    stdev_no_over_210nm(num_dots_index, 1) = stdev_vs_dist_vector(31);
    for replicate_index = 1:size(mean_bootstrap_vector_stack_cells, 1);
        mean_vs_dist_boot_vector = mean(mean_bootstrap_vector_stack_cells{replicate_index, num_dots_index}, 2);
        mean_no_over_boot_cells{num_dots_index, replicate_index} = mean_vs_dist_boot_vector;
        stdev_vs_dist_boot_vector = std(mean_bootstrap_vector_stack_cells{replicate_index, num_dots_index}, 0, 2);
        stdev_no_over_boot_cells{num_dots_index, replicate_index} = stdev_vs_dist_boot_vector;
        stdev_no_over_boot_0nm(num_dots_index, replicate_index) = stdev_vs_dist_boot_vector(1);
        stdev_no_over_boot_105nm(num_dots_index, replicate_index) = stdev_vs_dist_boot_vector(16);
        stdev_no_over_boot_210nm(num_dots_index, replicate_index) = stdev_vs_dist_boot_vector(31);
    end
end

% Load data from low overcounting
load([figure_path, data_name_low])

% Create stdev vs number dots curve of the radial averages
mean_low_over_cells = cell(size(mean_vector_stack_cells));
stdev_low_over_cells = cell(size(mean_vector_stack_cells));
stdev_low_over_0nm = zeros(size(mean_vector_stack_cells, 2), 1);
stdev_low_over_105nm = zeros(size(mean_vector_stack_cells, 2), 1);
stdev_low_over_210nm = zeros(size(mean_vector_stack_cells, 2), 1);
mean_low_over_boot_cells = cell(size(mean_bootstrap_vector_stack_cells));
stdev_low_over_boot_cells = cell(size(mean_bootstrap_vector_stack_cells));
stdev_low_over_boot_0nm = zeros(size(mean_bootstrap_vector_stack_cells, 2), size(mean_bootstrap_vector_stack_cells, 1));
stdev_low_over_boot_105nm = zeros(size(mean_vector_stack_cells, 2), size(mean_bootstrap_vector_stack_cells, 1));
stdev_low_over_boot_210nm = zeros(size(mean_vector_stack_cells, 2), size(mean_bootstrap_vector_stack_cells, 1));
for num_dots_index = 1:size(mean_vector_stack_cells, 2)
    mean_vs_dist_vector = mean(mean_vector_stack_cells{num_dots_index}, 2);
    mean_low_over_cells{num_dots_index} = mean_vs_dist_vector;
    stdev_vs_dist_vector = std(mean_vector_stack_cells{num_dots_index}, 0, 2);
    stdev_low_over_cells{num_dots_index} = stdev_vs_dist_vector;
    stdev_low_over_0nm(num_dots_index, 1) = stdev_vs_dist_vector(1);
    stdev_low_over_105nm(num_dots_index, 1) = stdev_vs_dist_vector(16);
    stdev_low_over_210nm(num_dots_index, 1) = stdev_vs_dist_vector(31);
    for replicate_index = 1:size(mean_bootstrap_vector_stack_cells, 1);
        mean_vs_dist_boot_vector = mean(mean_bootstrap_vector_stack_cells{replicate_index, num_dots_index}, 2);
        mean_low_over_boot_cells{num_dots_index, replicate_index} = mean_vs_dist_boot_vector;
        stdev_vs_dist_boot_vector = std(mean_bootstrap_vector_stack_cells{replicate_index, num_dots_index}, 0, 2);
        stdev_low_over_boot_cells{num_dots_index, replicate_index} = stdev_vs_dist_boot_vector;
        stdev_low_over_boot_0nm(num_dots_index, replicate_index) = stdev_vs_dist_boot_vector(1);
        stdev_low_over_boot_105nm(num_dots_index, replicate_index) = stdev_vs_dist_boot_vector(16);
        stdev_low_over_boot_210nm(num_dots_index, replicate_index) = stdev_vs_dist_boot_vector(31);
    end
end

% Load data from high overcounting
load([figure_path, data_name_high])

% Create stdev vs number dots curve of the radial averages
mean_high_over_cells = cell(size(mean_vector_stack_cells));
stdev_high_over_cells = cell(size(mean_vector_stack_cells));
stdev_high_over_0nm = zeros(size(mean_vector_stack_cells, 2), 1);
stdev_high_over_105nm = zeros(size(mean_vector_stack_cells, 2), 1);
stdev_high_over_210nm = zeros(size(mean_vector_stack_cells, 2), 1);
mean_high_over_boot_cells = cell(size(mean_bootstrap_vector_stack_cells));
stdev_high_over_boot_cells = cell(size(mean_bootstrap_vector_stack_cells));
stdev_high_over_boot_0nm = zeros(size(mean_bootstrap_vector_stack_cells, 2), size(mean_bootstrap_vector_stack_cells, 1));
stdev_high_over_boot_105nm = zeros(size(mean_vector_stack_cells, 2), size(mean_bootstrap_vector_stack_cells, 1));
stdev_high_over_boot_210nm = zeros(size(mean_vector_stack_cells, 2), size(mean_bootstrap_vector_stack_cells, 1));
for num_dots_index = 1:size(mean_vector_stack_cells, 2)
    mean_vs_dist_vector = mean(mean_vector_stack_cells{num_dots_index}, 2);
    mean_high_over_cells{num_dots_index} = mean_vs_dist_vector;
    stdev_vs_dist_vector = std(mean_vector_stack_cells{num_dots_index}, 0, 2);
    stdev_high_over_cells{num_dots_index} = stdev_vs_dist_vector;
    stdev_high_over_0nm(num_dots_index, 1) = stdev_vs_dist_vector(1);
    stdev_high_over_105nm(num_dots_index, 1) = stdev_vs_dist_vector(16);
    stdev_high_over_210nm(num_dots_index, 1) = stdev_vs_dist_vector(31);
    for replicate_index = 1:size(mean_bootstrap_vector_stack_cells, 1);
        mean_vs_dist_boot_vector = mean(mean_bootstrap_vector_stack_cells{replicate_index, num_dots_index}, 2);
        mean_high_over_boot_cells{num_dots_index, replicate_index} = mean_vs_dist_boot_vector;
        stdev_vs_dist_boot_vector = std(mean_bootstrap_vector_stack_cells{replicate_index, num_dots_index}, 0, 2);
        stdev_high_over_boot_cells{num_dots_index, replicate_index} = stdev_vs_dist_boot_vector;
        stdev_high_over_boot_0nm(num_dots_index, replicate_index) = stdev_vs_dist_boot_vector(1);
        stdev_high_over_boot_105nm(num_dots_index, replicate_index) = stdev_vs_dist_boot_vector(16);
        stdev_high_over_boot_210nm(num_dots_index, replicate_index) = stdev_vs_dist_boot_vector(31);
    end
end

% Create figure 0nm
figure
hold on
plot(number_dots_vector, stdev_no_over_0nm, 'k');
plot(number_dots_vector, stdev_low_over_0nm, 'b');
plot(number_dots_vector, stdev_high_over_0nm, 'g');
scatter(repmat(number_dots_vector.', size(stdev_no_over_boot_0nm, 2), 1), stdev_no_over_boot_0nm(:), 'k');
scatter(repmat(number_dots_vector.', size(stdev_low_over_boot_0nm, 2), 1), stdev_low_over_boot_0nm(:), 'b');
scatter(repmat(number_dots_vector.', size(stdev_high_over_boot_0nm, 2), 1), stdev_high_over_boot_0nm(:), 'g');
xlim([10,100]);
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
plot(number_dots_vector, stdev_no_over_105nm, 'k');
plot(number_dots_vector, stdev_low_over_105nm, 'b');
plot(number_dots_vector, stdev_high_over_105nm, 'g');
scatter(repmat(number_dots_vector.', size(stdev_no_over_boot_105nm, 2), 1), stdev_no_over_boot_105nm(:), 'k');
scatter(repmat(number_dots_vector.', size(stdev_low_over_boot_105nm, 2), 1), stdev_low_over_boot_105nm(:), 'b');
scatter(repmat(number_dots_vector.', size(stdev_high_over_boot_105nm, 2), 1), stdev_high_over_boot_105nm(:), 'g');
xlim([10,100]);
xlabel('Number of dots per radial average');
ylabel('Standard deviation of radial average value');
title('105 nm radial average');
hold off
figure_name = [figure_path, figure_prefix, '105nm.png'];
print(gcf, figure_name, '-dpng');
delete(gcf);

% Create figure 205nm
figure
hold on
plot(number_dots_vector, stdev_no_over_210nm, 'k');
plot(number_dots_vector, stdev_low_over_210nm, 'b');
plot(number_dots_vector, stdev_high_over_210nm, 'g');
scatter(repmat(number_dots_vector.', size(stdev_no_over_boot_210nm, 2), 1), stdev_no_over_boot_210nm(:), 'k');
scatter(repmat(number_dots_vector.', size(stdev_low_over_boot_210nm, 2), 1), stdev_low_over_boot_210nm(:), 'b');
scatter(repmat(number_dots_vector.', size(stdev_high_over_boot_210nm, 2), 1), stdev_high_over_boot_210nm(:), 'g');
xlim([10,100]);
xlabel('Number of dots per radial average');
ylabel('Standard deviation of radial average value');
title('205 nm radial average');
hold off
figure_name = [figure_path, figure_prefix, '205nm.png'];
print(gcf, figure_name, '-dpng');
delete(gcf);

% Create figure 10 dots, no overcounting
figure
hold on
for replicate_index = 1:size(mean_bootstrap_vector_stack_cells, 1);
    errorbar(distance_vector * 7, mean_no_over_boot_cells{1, replicate_index}, stdev_no_over_boot_cells{1, replicate_index}, 'Color', [.4, .4, .4]);
end
errorbar(distance_vector * 7, mean_no_over_cells{1}, stdev_no_over_cells{1}, 'Color', [1, 0, 0]);
xlim([0,500]);
xlabel('Radial distance (nm)');
ylabel('Normalized correlation value');
title({'True means vs bootstrap estimates','10 dots averaged per value, no overcounting'});
hold off
figure_name = [figure_path, figure_prefix_2, '10_dots_no_over.png'];
print(gcf, figure_name, '-dpng');
delete(gcf);

% Create figure 30 dots, no overcounting
figure
hold on
for replicate_index = 1:size(mean_bootstrap_vector_stack_cells, 1);
    errorbar(distance_vector * 7, mean_no_over_boot_cells{3, replicate_index}, stdev_no_over_boot_cells{3, replicate_index}, 'Color', [.4, .4, .4]);
end
errorbar(distance_vector * 7, mean_no_over_cells{3}, stdev_no_over_cells{3}, 'Color', [1, 0, 0]);
xlim([0,500]);
xlabel('Radial distance (nm)');
ylabel('Normalized correlation value');
title({'True means vs bootstrap estimates','30 dots averaged per value, no overcounting'});
hold off
figure_name = [figure_path, figure_prefix_2, '30_dots_no_over.png'];
print(gcf, figure_name, '-dpng');
delete(gcf);

% Create figure 100 dots, no overcounting
figure
hold on
for replicate_index = 1:size(mean_bootstrap_vector_stack_cells, 1);
    errorbar(distance_vector * 7, mean_no_over_boot_cells{6, replicate_index}, stdev_no_over_boot_cells{6, replicate_index}, 'Color', [.4, .4, .4]);
end
errorbar(distance_vector * 7, mean_no_over_cells{6}, stdev_no_over_cells{6}, 'Color', [1, 0, 0]);
xlim([0,500]);
xlabel('Radial distance (nm)');
ylabel('Normalized correlation value');
title({'True means vs bootstrap estimates','100 dots averaged per value, no overcounting'});
hold off
figure_name = [figure_path, figure_prefix_2, '100_dots_no_over.png'];
print(gcf, figure_name, '-dpng');
delete(gcf);

% Create figure 10 dots, low overcounting
figure
hold on
for replicate_index = 1:size(mean_bootstrap_vector_stack_cells, 1);
    errorbar(distance_vector * 7, mean_low_over_boot_cells{1, replicate_index}, stdev_low_over_boot_cells{1, replicate_index}, 'Color', [.4, .4, .4]);
end
errorbar(distance_vector * 7, mean_low_over_cells{1}, stdev_low_over_cells{1}, 'Color', [1, 0, 0]);
xlim([0,500]);
xlabel('Radial distance (nm)');
ylabel('Normalized correlation value');
title({'True means vs bootstrap estimates','10 dots averaged per value, low overcounting'});
hold off
figure_name = [figure_path, figure_prefix_2, '10_dots_low_over.png'];
print(gcf, figure_name, '-dpng');
delete(gcf);

% Create figure 30 dots, low overcounting
figure
hold on
for replicate_index = 1:size(mean_bootstrap_vector_stack_cells, 1);
    errorbar(distance_vector * 7, mean_low_over_boot_cells{3, replicate_index}, stdev_low_over_boot_cells{3, replicate_index}, 'Color', [.4, .4, .4]);
end
errorbar(distance_vector * 7, mean_low_over_cells{3}, stdev_low_over_cells{3}, 'Color', [1, 0, 0]);
xlim([0,500]);
xlabel('Radial distance (nm)');
ylabel('Normalized correlation value');
title({'True means vs bootstrap estimates','30 dots averaged per value, low overcounting'});
hold off
figure_name = [figure_path, figure_prefix_2, '30_dots_low_over.png'];
print(gcf, figure_name, '-dpng');
delete(gcf);

% Create figure 100 dots, low overcounting
figure
hold on
for replicate_index = 1:size(mean_bootstrap_vector_stack_cells, 1);
    errorbar(distance_vector * 7, mean_low_over_boot_cells{6, replicate_index}, stdev_low_over_boot_cells{6, replicate_index}, 'Color', [.4, .4, .4]);
end
errorbar(distance_vector * 7, mean_low_over_cells{6}, stdev_low_over_cells{6}, 'Color', [1, 0, 0]);
xlim([0,500]);
xlabel('Radial distance (nm)');
ylabel('Normalized correlation value');
title({'True means vs bootstrap estimates','100 dots averaged per value, low overcounting'});
hold off
figure_name = [figure_path, figure_prefix_2, '100_dots_low_over.png'];
print(gcf, figure_name, '-dpng');
delete(gcf);

% Create figure 10 dots, high overcounting
figure
hold on
for replicate_index = 1:size(mean_bootstrap_vector_stack_cells, 1);
    errorbar(distance_vector * 7, mean_high_over_boot_cells{1, replicate_index}, stdev_high_over_boot_cells{1, replicate_index}, 'Color', [.4, .4, .4]);
end
errorbar(distance_vector * 7, mean_high_over_cells{1}, stdev_high_over_cells{1}, 'Color', [1, 0, 0]);
xlim([0,500]);
xlabel('Radial distance (nm)');
ylabel('Normalized correlation value');
title({'True means vs bootstrap estimates','10 dots averaged per value, high overcounting'});
hold off
figure_name = [figure_path, figure_prefix_2, '10_dots_high_over.png'];
print(gcf, figure_name, '-dpng');
delete(gcf);

% Create figure 30 dots, high overcounting
figure
hold on
for replicate_index = 1:size(mean_bootstrap_vector_stack_cells, 1);
    errorbar(distance_vector * 7, mean_high_over_boot_cells{3, replicate_index}, stdev_high_over_boot_cells{3, replicate_index}, 'Color', [.4, .4, .4]);
end
errorbar(distance_vector * 7, mean_high_over_cells{3}, stdev_high_over_cells{3}, 'Color', [1, 0, 0]);
xlim([0,500]);
xlabel('Radial distance (nm)');
ylabel('Normalized correlation value');
title({'True means vs bootstrap estimates','30 dots averaged per value, high overcounting'});
hold off
figure_name = [figure_path, figure_prefix_2, '30_dots_high_over.png'];
print(gcf, figure_name, '-dpng');
delete(gcf);

% Create figure 100 dots, high overcounting
figure
hold on
for replicate_index = 1:size(mean_bootstrap_vector_stack_cells, 1);
    errorbar(distance_vector * 7, mean_high_over_boot_cells{6, replicate_index}, stdev_high_over_boot_cells{6, replicate_index}, 'Color', [.4, .4, .4]);
end
errorbar(distance_vector * 7, mean_high_over_cells{6}, stdev_high_over_cells{6}, 'Color', [1, 0, 0]);
xlim([0,500]);
xlabel('Radial distance (nm)');
ylabel('Normalized correlation value');
title({'True means vs bootstrap estimates','100 dots averaged per value, high overcounting'});
hold off
figure_name = [figure_path, figure_prefix_2, '100_dots_high_over.png'];
print(gcf, figure_name, '-dpng');
delete(gcf);

