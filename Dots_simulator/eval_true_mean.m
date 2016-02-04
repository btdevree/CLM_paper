function [distance_vector, mean_vector] = eval_true_mean(params, number_dots, max_correlation_radius, dots_per_cell_mean,...
    dot_radius, dot_correlation_value, label_density_mean, label_density_stdev, label_SN_ratio, event_overcounting,...
    dot_center_precision, event_precision, STORM_pixel_resolution, method_type)
%EVAL_TRUE_MEAN Funciton for parallel evaluation of a true-mean dots
%simulation

% Get results
[correlation, ~] = collect_dot_correlations(params, number_dots, max_correlation_radius, dots_per_cell_mean,...
    dot_radius, dot_correlation_value, label_density_mean, label_density_stdev, label_SN_ratio, event_overcounting,...
    dot_center_precision, event_precision, STORM_pixel_resolution, method_type);

% Get a radial average
[distance_vector, mean_vector, ~, ~] = radial_average_2D_correlation(correlation);

fprintf('.');
end

