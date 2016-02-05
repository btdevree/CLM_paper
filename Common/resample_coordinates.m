function [ resampled_coordinate_cells ] = resample_coordinates( coordinate_cells )
%RESAMPLE_COORDINATES Resamples coordinates with replacement.
%
% Inputs:
%   coordinate_cells: Cell array containing matrix with m rows of 
%       coordinates.
% Outputs:
%   resampled_coordinate_cells: Cell array containing resampled 
%       coordinates. Returns the same number of coordinates across 
%       all cells, but each individual cell may have a different number 
%       than before. 

% Rewrite the coordinate list as a single matrix with the cell index recorded
number_coord_dims = size(coordinate_cells{1}, 2);
current_index = 0;
for cell_index = 1:length(coordinate_cells(:))
    start_index = current_index + 1;
    end_index = current_index + size(coordinate_cells{cell_index}, 1);
    coord_list(start_index:end_index, 1:number_coord_dims) = coordinate_cells{cell_index};
    coord_list(start_index:end_index, number_coord_dims + 1) = cell_index;
    current_index = end_index;
end
   
% Resample the list
resampled_indices = randi(size(coord_list, 1), size(coord_list, 1), 1);
resampled_coord_list = cat(1, coord_list(resampled_indices, :));

% Put the lists back into cells
resampled_coordinate_cells = cell(size(coordinate_cells));
for cell_index = 1:length(coordinate_cells(:))
    selection_vector = resampled_coord_list(:, end) == cell_index;
    resampled_coordinate_cells{cell_index} = resampled_coord_list(selection_vector, 1:number_coord_dims);
end 
end

