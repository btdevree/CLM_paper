function [ECI_mean, ECI_stdev, TCI_mean, TCI_stdev] = calc_ECI_and_TCI(IIC_results, ideal_discrepency_results, fraction_vector)
%CALC_ECI_AND_TCI Calculate the ECI and TCI values from a set of IIC
%   curves.
%
% Inputs:
%   IIC_results: 2D or 3D array of IIC results calculated by 
%       generate_IIC_curve, arranged with the fraction of events used on 
%       the first dimention, number of events on the second dimension, and 
%       replicates on the third dimension.
%   ideal_discrepency: raw ideal image discrepency values, arranged in the 
%       same order as the IIC_results.
%   fraction_vector: column vector of increacing fractions of included
%       points used to make up the IIC curve.
%
% Outputs:
%   ECI_mean/stdev: mean and standard deviation of ECI value for each
%       number of events, returned as a column vector.
%   TCI_mean/stdev: mean and standard deviation of TCI value for each
%       number of events, returned as a column vector.
 
% Calculate the AOC and ECI of the sum of squares IIC curves
dFrac = fraction_vector(2:end) - fraction_vector(1:end-1);
midpoint_II = (IIC_results(2:end, :, :) + IIC_results(1:end-1, :, :))/2;
AOC = sum(repmat(dFrac, 1, size(midpoint_II, 2), size(midpoint_II, 3)) .* midpoint_II, 1);
ECI_data = (2 * AOC - 1);
ECI_mean = squeeze(mean(ECI_data, 3))';
ECI_stdev = squeeze(std(ECI_data, 0, 3))';

% Calculate the corrosponding TCI 
TCI_data = 1 - (ideal_discrepency_results(end, :, :) ./ ideal_discrepency_results(1, :, :));
TCI_mean = squeeze(mean(TCI_data, 3))';
TCI_stdev = squeeze(std(TCI_data, 0, 3))';    
end

