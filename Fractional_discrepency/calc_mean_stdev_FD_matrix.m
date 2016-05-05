function [ mean_FD_values, stdev_FD_values ] = calc_mean_stdev_FD_matrix( FD_matrix )
%CALC_MEAN_STDEV_FDY_MATRIX Calculates the mean and standard
%deviation of a 4-D fractional discrepency matrix.

% Initialize final matricies
mean_FD_values = zeros(size(FD_matrix, 1), size(FD_matrix, 4));
stdev_FD_values = zeros(size(FD_matrix, 1), size(FD_matrix, 4));

% Loop through each set of fractions and event numbers
for fraction_index = 1:size(FD_matrix, 1);
    for event_num_index = 1:size(FD_matrix, 4);
    
        % Take the mean and standard deviation of all the replicates and measurement repeats
        values = FD_matrix(fraction_index, :, :, event_num_index);
        mean_FD_values(fraction_index, event_num_index) = mean(values(:));
        stdev_FD_values(fraction_index, event_num_index) = std(values(:)); 
    end
end
end

