% Make plots of the xcor_data files

clear all

% Set directory information
path_parts = strsplit(pwd, 'CLM_paper');
filepath = [path_parts{1}, 'CLM_figures_and_data/'];
data_name = 'aliasing_part_B_bins.mat';

% Load a parameter structure from current dirctory with deisred settings
load([filepath, data_name]) % loads as 'params'

xcor_mean = cell(size(xcor_curves, 1), 1);
xcor_stdev = cell(size(xcor_curves, 1), 1);
resid_mean = cell(size(residuals, 1), 1);
resid_stdev = cell(size(residuals, 1), 1);
% Calculate means and standard deviations
for pixel_length_index = 1:length(pixel_lengths)
    xcor_stack = zeros(length(dist_vectors{pixel_length_index}), size(xcor_curves, 2));
    residual_stack = zeros(length(dist_vectors{pixel_length_index}), size(residuals, 2));
    for replicate_index = 1:size(xcor_curves, 2)
        xcor_stack(:, replicate_index) = xcor_curves{pixel_length_index, replicate_index};
        residual_stack(:, replicate_index) = residuals{pixel_length_index, replicate_index};
    end
    xcor_mean{pixel_length_index}  = mean(xcor_stack, 2);
    xcor_stdev{pixel_length_index} = std(xcor_stack, 0, 2);
    resid_mean{pixel_length_index}  = mean(residual_stack, 2);
    resid_stdev{pixel_length_index} = std(residual_stack, 0, 2); 
end 

% Create figure
figure
offset = .5;
hold on
for replicate_index = 1:size(xcor_mean, 1);
    errorbar(dist_vectors{replicate_index}, xcor_mean{replicate_index} + offset*(replicate_index-1), xcor_stdev{replicate_index}, 'Color', [.5, .5, .5]);
    plot(dist_vectors{replicate_index}, analytical_curves{replicate_index} + offset*(replicate_index-1), 'r');
end
xlim([0,1000]);
xlabel('Radial distance (nm)');
ylabel('Normalized correlation value');
title({'Alaising effects of pixel size'});
hold off
figure_name = [filepath, 'pixel_size_alising_offset_graph_0nm_bins.png'];
print(gcf, figure_name, '-dpng');
delete(gcf);

% Create figure
figure
offset = .5;
hold on
for replicate_index = 1:size(xcor_mean, 1);
    errorbar(dist_vectors{replicate_index}, resid_mean{replicate_index} + offset*(replicate_index-1), resid_stdev{replicate_index}, 'Color', [.5, .5, .5]);
end
xlim([0,1000]);
xlabel('Radial distance (nm)');
ylabel('Residual correlation value');
title({'Residuals between expected and actual crosscorrelaiton curves'});
hold off
figure_name = [filepath, 'pixel_size_alising_residual_graph_0nm_bins.png'];
print(gcf, figure_name, '-dpng');
delete(gcf);