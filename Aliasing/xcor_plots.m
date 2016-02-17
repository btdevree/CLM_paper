% Make plots of the xcor_data files

% Set directory information
path_parts = strsplit(pwd, 'CLM_paper');
filepath = [path_parts{1}, 'CLM_figures_and_data/'];
data_name = 'aliasing_part_B.mat';

% Load a parameter structure from current dirctory with deisred settings
clear all
load([filepath, data_name]) % loads as 'params'

% Make new figure
figure;

% Plot and save average figure
imagesc(av_xcor);
colorbar;
axis square;
print(gcf, 'average_xcor', '-dpng');
clf

% Plot and save cell figures
for cell_index = 1:length(xcor_data)
    
    subplot(2,3,cell_index);
    imagesc(xcor_data{cell_index}.xcor_image);
    axis square;
end
print(gcf, 'cell_xcor', '-dpng');
clf

%Plot and save the radial average
% Error bars are pretty dense, so we'll only plot 1 out of 5
counts = 0:length(dist_vec);
ind = mod(counts, 5) == 0;
errorbar(dist_vec(ind), mean_vec(ind), stdev_vec(ind), 'ko-', 'markerfacecolor', 'k', 'markersize', 2);
hold on
errorbar(dist_vec(ind), mean_vec(ind), sem_vec(ind), 'bo-', 'markerfacecolor', 'b', 'markersize', 2);
plot(dist_vec, mean_vec, 'r-');
set(gca, 'xlim', [0, max(dist_vec(:))]);
set(gcf, 'Position', [100, 100, 800, 400])
set(gcf, 'PaperPosition', [0.5, 0.5, 8, 4])
hold off
print(gcf, 'radial_average', '-dpng');
