function [distance_vector, permuted_mean_vectors] = eval_bootstrap_mean(number_permutations, params, number_dots, max_correlation_radius, dots_per_cell_mean,...
    dot_radius, dot_correlation_value, label_density_mean, label_density_stdev, label_SN_ratio, event_overcounting,...
    dot_center_precision, event_precision, STORM_pixel_resolution, method_type)
%EVAL_BOOTSTRAP_MEAN Function for parallel evaluation of a boostrap-mean dots
%simulation

% Get a results stack
[correlation_stack, ~] = bootstrap_dot_correlation(params, number_permutations, number_dots, max_correlation_radius, dots_per_cell_mean,...
    dot_radius, dot_correlation_value, label_density_mean, label_density_stdev, label_SN_ratio, event_overcounting,...
    dot_center_precision, event_precision, STORM_pixel_resolution, method_type);

% Radially average the correlations
permuted_mean_vectors = zeros((size(correlation_stack, 1)-1)/2 + 1, number_permutations);
for permutation_index = 1:number_permutations
   
    % Get a radial average
    [distance_vector, mean_vector, ~, ~] = radial_average_2D_correlation(correlation_stack(:, :, permutation_index));
    
    % Add to the results stack
    permuted_mean_vectors(:, permutation_index) = mean_vector;
end
end

