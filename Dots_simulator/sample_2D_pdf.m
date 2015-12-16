function [ xy_samples ] = sample_2D_pdf( number_samples, pdf_image, resolution )
%SAMPLE_2D_PDF Samples 2D random variables from an arbitrary 2D pdf image.
%
% Uses linear interpolation to determine sampling rates between the pixel
%   centers. Assumes a Cartesian coordinate system with the origin in the 
%   bottom left corner of the bottom left pixel. Algorithm employs inverse 
%   sampling of a greyscale dilated image to sample points at a slightly 
%   higher rate than needed and then corrects the rate with rejection 
%   sampling. 
%
% Inputs:
%   number_samples: number of samples to return.
%   pdf_image: array of floating-point doubles, a discrete pdf to be 
%       sampled from. Does not have to be normalized.
%   resolution: optional, the number of measurment units per pixel. 
%       Default = 1.
% Outputs:
%   xy_samples: n by 2 array to floating point doubles. Given in the
%       measurement units supplied in the resolution parameter.

% Set defaults
if nargin < 3; resolution = 1; end;

% Create a dilated image
se = strel('square',3);
enclosing_pdf = imdilate(pdf_image, se);

% Average dilated and original pdfs to get the maximum sampling rate that will be needed for each pixel (given linear interpolation)
enclosing_pdf = (enclosing_pdf + pdf_image) / 2;

% Calculate the expected rejection rate in rejected points per created point
rejection_rate = 1 - (sum(pdf_image(:)) / sum(enclosing_pdf(:)));

% Calculate the cumulative sum of the enclosing pdf image unraveled as a 1D array
enclosing_cdf = [0; cumsum(enclosing_pdf(:))];
enclosing_cdf = enclosing_cdf / enclosing_cdf(end); % normalize to range [0,1]


end

