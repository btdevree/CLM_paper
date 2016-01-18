function [distance_vector, mean_vector] = eval_bootstrap_mean(number_permutations, params, number_dots, max_correlation_radius, dots_per_cell_mean,...
    dot_radius, dot_correlation_value, label_density_mean, label_density_stdev, label_SN_ratio, event_overcounting,...
    dot_center_precision, event_precision, STORM_pixel_resolution, method_type)
%EVAL_BOOTSTRAP_MEAN Function for parallel evaluation of a boostrap-mean dots
%simulation

% Get a results stack
[correlation_stack, ~] = collect_dot_correlations(params, number_dots, max_correlation_radius, dots_per_cell_mean,...
    dot_radius, dot_correlation_value, label_density_mean, label_density_stdev, label_SN_ratio, event_overcounting,...
    dot_center_precision, event_precision, STORM_pixel_resolution, method_type);

% Repeat getting the mean with resampling
permuted_image_stack = zeros(2*max_radius_px + 1, 2*max_radius_px + 1, number_permutations);
for permutation_index = 1:number_permutations
    
    % Permute the stack
    permuted_indices = randi(size(correlation_stack, 3));
    permuted_stack = cat(3, correlation_stack(:, :, permuted_indices));
    permuted_image = mean(permuted_stack, 3);
    
    % Add to the results stack
    permuted_image_stack(:, :, permutation_index) = permuted_image;
end

% Average the stack
bootstrap_image = mean(permuted_image_stack, 3);

% Get a radial average
[distance_vector, mean_vector, ~, ~] = radial_average_2D_correlation(bootstrap_image);
end

