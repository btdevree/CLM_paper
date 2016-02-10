% Make plots of the xcor_data files
% Expects:
%   xcor_data: cell array of XC_data structures
%   av_xcor: average crosscorrelation image
%   dist_vec: radial distance vector
%   mean_vec: radial mean vector
%   stdev_vec: radial standard deviation vector

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
