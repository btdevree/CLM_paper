function [ image ] = Gauss_STORM_image(resolution, xy_data, covar_inv, covar_det, calc_cutoff_pixels_x, calc_cutoff_pixels_y, xmesh, ymesh)
%Gauss_STORM_image Creates a STORM image from the list of (x,y) data; 
%   mirrors the C MEX function Gauss_STORM_image_MEX.cpp. 
%
% Creates a STORM image as a matrix of floating-point doubles. Assumes a
% Cartesian coordinate system. 
%
% Inputs
%   resolution: floating-point double, number of nanometers per final image
%       pixel
%   xy_data: n by 2 array of floating-point doubles, center of each
%       gaussian pdf
%   covar_inv:  2 by 2 array of floating-point doubles, the inverse of the
%       2D gaussian covariance matrix
%   covar_inv: floating-point double, the determinant of the 2D gaussian
%       covariance matrix
%   calc_cutoff_pixels_x\y: 32 bit integer, extent in number of pixels on 
%       either side of the center pixel to calculate the pdf distribution
%       out to
%   total_number_pixels_x\y: 32 bit integer, number of pixels in the x and
%       y dimension of the final image.
%   xmesh\ymesh: array of floating-point doubles with the coordinate of
%       the center of each pixel, generated with meshgrid.

% Create the output image
image = zeros(size(xmesh));

% 



end

