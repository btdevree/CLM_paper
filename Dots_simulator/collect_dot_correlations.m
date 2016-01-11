function [correlation_stack, number_cells] = collect_dot_correlations(param_struct, number_dots, max_corr_length, dots_per_cell_mean,...
    dots_per_cell_stdev, dot_radius, dot_correlation_value, label_density_mean, label_density_stdev, label_SN_ratio, event_overcounting,...
    dot_center_precision, event_precision, STORM_pixel_resolution, STORM_method)
%COLLECT_DOT_CORRELATIONS Returns a stack of 2D correlations from simulated
%   cells with "dot" type well-defined centers of interest
%
%   This function handles the logic needed to run a dot simulation with the
%   given parameters. 
%   
% Inputs:
%   param_struct: structure that holds the parameters needed for simulating
%       a cell with the create_test_data/create_test_STORM_image
%       functions. Parameters in this structure may be superceded due to 
%       requirements of the other inputs.
%   number_dots: the number of dots to simulate
%   max_corr_length: the maximum radius out to which the correlation is 
%       calculated, given in nanometers. 
%   dots_per_cell_mean: mean number of dots to simulate per each simulated
%       cell.
%   dots_per_cell_stdev: standard deviation of the number of dots in each
%       simulated cell.
%   dot_radius: radius of each dot region, given in nanometers.
%   dot_correlation_value: factor for how more/less frequently events are
%       found in the dot regions
%   label_density_mean: mean density of true events in the main body of the 
%      cell, given in events/micrometer^2
%   label_density_stdev: standard deviation of the density of true events 
%       in the main body of the cell, given in events/micrometer^2
%   label_SN_ratio: signal to noise ratio of labels, given as # true
%       events / # spurious events.
%   event_overcounting: average number of times each event is overcounted. 
%   STORM_pixel_resolution: size of pixels in the STORM images that are 
%       created, given in nanometers.
%   STORM_method: string indicating the method used to create the STORM
%       images. Options are: 'pdf' or 'binning'
% Outputs:
%   correlation_stack: 3D stack of the 2D correlations for each dot region.
%   number_cells: the number of cells required to collect the requested
%       number of correlations.

% Calculate required constant values
max_radius_px = ceil(max_corr_length / STORM_pixel_resolution);

% Initialize the results matrix
correlation_stack = zeros(2 * max_radius_px + 1, 2 * max_radius_px + 1, number_dots);

% Edit param values as needed

% Create and correlate cells until we have collected enough correlations
while number_dots_complete < number_dots
    
end


end

