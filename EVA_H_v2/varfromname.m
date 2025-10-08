% 	Written by Thomas J. Aubry, July 2025.
% 	Department of Earth Sciences, University of Oxford, UK
%   E-mail: thomas.aubry@earth.ox.ac.uk
% 	Please cite the corresponding paper if you use or modify this script

function var=varfromname(filename,varname,fv)
%fv = boolean for whether to replace fill value
%filename and varname = file and variable name (string inputs)

%open netcdf and read var
ncid = netcdf.open(filename);
var = double(netcdf.getVar(ncid,netcdf.inqVarID(ncid,varname)));

%replace fill values by NaN if requested
if fv==1
    fillval=netcdf.getAtt(ncid,netcdf.inqVarID(ncid,varname),'_FillValue');
    var(var(:)==fillval)=NaN;
end

end