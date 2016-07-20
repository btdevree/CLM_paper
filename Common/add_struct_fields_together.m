function [new_struct] = add_struct_fields_together(struct_1, struct_2, cat_dimension)
% Add together the fields in two structures and returns a larger structure
%
% Input:
%   struct_1/2: structures with identical fields that contain matricies to
%       be stacked together.
%   cat_dimension: Dimension along which to stack the results. Optional, 
%       default = 1.
% Output:
%   new_struct: structures with the same fields and stacked entries.

% preallocate a structure.
new_struct = struct();

% Get fieldnames and loop through each name
fn = fieldnames(struct_1);
for field_index = 1:size(fn, 1)

    % contatinate field and put in new structure
    current_field = fn{field_index};
    data_1 = struct_1.(current_field);
    data_2 = struct_2.(current_field);
    new_struct.(current_field) = cat(cat_dimension, data_1, data_2);
end
end
