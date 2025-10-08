% 	Written by Thomas J. Aubry, July 2025.
% 	Department of Earth Sciences, University of Oxford, UK
%   E-mail: thomas.aubry@earth.ox.ac.uk
% 	Please cite the corresponding paper if you use or modify this script



function  emis_sanchecks_num(variable,filled,range)
%this function simply checks that the provided variable is within the
%provided range. If filled is 1, it means all values should be filled in
%which case the function also checks there is no NaN value

if filled==1 & sum(isnan(variable))>0
    error('the variable should not have NaN values')
end

if max(variable)>max(range)
    error('the variable maximum is outside range')
end

if min(variable)<min(range)
    error('the variable minimum is outside range')
end