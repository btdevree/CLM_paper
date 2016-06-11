function [pdf_map, line_start_coords, line_end_cords] = straight_lines_pdf_map(parameter_struct, number_of_lines, line_width, line_to_background_ratio, line_min_length, line_max_length)
%STRAIGHT_LINES_PDF_MAP Makes a pdf map of line structures in the field of
% view.
%
% Each line is randomly placed within the image,
%
% Inputs:
%   parameter_struct: parameter structure made with
%       test_movie_parameters_dv.
%   number_of_lines: integer, number of lines to put in the pdf map.
%   line_width: width of the lines, in nanometers.
%   line_to_background_ratio: ratio of the density of events inside the
%       line region to the general background. 
%   line_min_length: minimum length of the lines in the image, given in
%       nanometers or as a string in the format 'xxx%' representing the 
%       percentage of the smallest image dimenstion. Optional, default =
%       '5%'.
%   line_max_length: maximum length of the lines in the image, given in
%       nanometers or as a string in the format 'xxx%' representing the 
%       percentage of the smallest image dimenstion. Optional, default =
%       '95%'.
% Outputs:
%   pdf_map: array of floating-point doubles, normalized sampling of the
%       analytical pdf.
%   line_start_coords/line_end_cords: coordinates of the start and end of 
%       the lines in the pdf map, given as an n by 2 matrix of (x, y)
%       values.

% Rename parameter structure for convenience
params = parameter_struct;

map_resolution = params.ch1_distribution_params{2};
end

