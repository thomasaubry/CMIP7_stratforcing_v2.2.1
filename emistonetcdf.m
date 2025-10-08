% 	Written by Thomas J. Aubry, July 2025.
% 	Department of Earth Sciences, University of Oxford, UK
%   E-mail: thomas.aubry@earth.ox.ac.uk
% 	Please cite the corresponding paper if you use or modify this script



%NOTE: All netcdf written in this script are not fully formatted to CMIP standard (e.g. CF-compliance, CMIP
%compliance, dimension orders, unit, etc). An additional python script,
%originally developped by Zebedee Nicholls, was used to reformat the netcdf
%produced by this script into the ones provided to the community. These
%python scripts can be provided upon request.


function emistonetcdf(xlsxname,ncfilename)

%inputs are a clean xlsx file you want to readily convert into a netcdf one
%delete any netcdf file with same name to overwrite
delete ncfilename

%read xlsx; em_all contains all cells, em_double mumerical data as double
[em_double em_text em_all]=xlsread(xlsxname);


%get headers
em_header=em_all(1,:);
%get data; if more than one header row this function won't work and I did
%not embed an error message...
em_data=em_all(2:end,:);

%data length (number of eruption in our case
Neru=length(em_data(:,1));
%define dummy eruption number variable
eruption=(1:Neru)';


ncid = netcdf.create(ncfilename,'NETCDF4');
netcdf.close(ncid);

% =========================================================================
%% Write global attribute
%==========================================================================
%write global attribute


ncwriteatt(ncfilename,'/','activity_id','input4MIPs')
ncwriteatt(ncfilename,'/','activity_id','input4MIPs')
ncwriteatt(ncfilename,'/','contact','CMIP7 Forcing Task Team: volcano dataset providers (thom.aubry@gmail.com,Vaishali.Naik@noaa.gov,durack1@llnl.gov)')
ncwriteatt(ncfilename,'/','creation_date',strcat(datestr(now, 'yyyy-mm-dd'),'T',datestr(now, 'HH:MM:SS'),'Z'))
ncwriteatt(ncfilename,'/','data_specs_version','01.00')
ncwriteatt(ncfilename,'/','dataset_category','emissions')
ncwriteatt(ncfilename,'/','dataset_category','Volcanic sulfur emissions (upper-troposphere and above)')
ncwriteatt(ncfilename,'/','frequency','Ill-defined for this; ask advise to CFTT leadership')
ncwriteatt(ncfilename,'/','further_info_url','Link to documentation to be added')
ncwriteatt(ncfilename,'/','grid','Ill-defined for this; ask advise to CFTT leadership')
ncwriteatt(ncfilename,'/','grid_label','Ill-defined for this; ask advise to CFTT leadership')
ncwriteatt(ncfilename,'/','history','')
ncwriteatt(ncfilename,'/','institution','tbd after further discussion with Paul')
ncwriteatt(ncfilename,'/','institution_id','tbd after further discussion with Paul')
ncwriteatt(ncfilename,'/','mip_era','CMIP6plus')
ncwriteatt(ncfilename,'/','product','observations and proxies')
ncwriteatt(ncfilename,'/','realm','upper troposphere and above')
ncwriteatt(ncfilename,'/','references','to be updated')
ncwriteatt(ncfilename,'/','region','global upper atmosphere')
ncwriteatt(ncfilename,'/','release_year','2024')
ncwriteatt(ncfilename,'/','source','MSVOLSO2L4, eVolv2k, Global Volcanism Programme VOTW, IVESPA, and other sources as specified in attributes and documentation')
ncwriteatt(ncfilename,'/','source_id','tbd after further discussion with Paul')
ncwriteatt(ncfilename,'/','source_type','satellite, ice-core and geological records')
ncwriteatt(ncfilename,'/','target_mip','CMIP6Plus')


%write each dimension/variable
%NOTE: if variables name not ordered as assume in this script it will break
%or produce erroneous file...not very elegant again!

% =========================================================================
%% Eruption dimension
%==========================================================================
varname='eruption';
%create new variable
nccreate(ncfilename,varname,'Dimensions',{'eruption',Neru})
%write new variable
ncwrite(ncfilename,varname,eruption);
%write variable attribute
ncwriteatt(ncfilename,varname,'Unit','Arbitrary chronological eruption number')
ncwriteatt(ncfilename,varname,'long_name','Arbitrary chronological eruption number');
ncwriteatt(ncfilename,varname,'standard_name','eruption number');


% =========================================================================
%% Volcano name
%==========================================================================
i=find(strcmp(em_header,'Volcano name'));
varname='volcano_name';

%create new variable
nccreate(ncfilename,varname,'Dimensions',{'eruption',Neru},'Datatype','string','format','netcdf4')
%write new variable
ncwrite(ncfilename,varname,string(em_data(:,i)));
% %write variable attribute
ncwriteatt(ncfilename,varname,'Unit','NA')
ncwriteatt(ncfilename,varname,'long_name','Volcano name');
ncwriteatt(ncfilename,varname,'standard_name','volcano name');


% =========================================================================
%% Eruption year
%==========================================================================
i=find(strcmp(em_header,'Eruption year'));
varname='eruption_year';

%create new variable
nccreate(ncfilename,varname,'Dimensions',{'eruption',Neru})
%write new variable
ncwrite(ncfilename,varname,cell2mat(em_data(:,i)));
%write variable attribute
ncwriteatt(ncfilename,varname,'Unit','Year CE')
ncwriteatt(ncfilename,varname,'long_name','Eruption year');
ncwriteatt(ncfilename,varname,'standard_name','eruption year');



% =========================================================================
%% Eruption month
%==========================================================================
i=find(strcmp(em_header,'Eruption month (1-12)'));
varname='eruption_month';

%create new variable
nccreate(ncfilename,varname,'Dimensions',{'eruption',Neru})
%write new variable
ncwrite(ncfilename,varname,cell2mat(em_data(:,i)));
%write variable attribute
ncwriteatt(ncfilename,varname,'Unit','Month (1-12)')
ncwriteatt(ncfilename,varname,'long_name','Eruption month');
ncwriteatt(ncfilename,varname,'standard_name','eruption month');
%add flag
i=find(strcmp(em_header,'Eruption month fill flag (0=not filled, 1=fill value)'));
flag=cell2mat(em_data(:,i));
ncwriteatt(ncfilename,varname,'Fill flag (0=not filled, 1=fill value)',flag);




% =========================================================================
%% Eruption day
%==========================================================================
i=find(strcmp(em_header,'Eruption day (1-31)'));
varname='eruption_day';

%create new variable
nccreate(ncfilename,varname,'Dimensions',{'eruption',Neru})
%write new variable
ncwrite(ncfilename,varname,cell2mat(em_data(:,i)));
%write variable attribute
ncwriteatt(ncfilename,varname,'Unit','Day of month (1-31)')
ncwriteatt(ncfilename,varname,'long_name','Eruption day');
ncwriteatt(ncfilename,varname,'standard_name','eruption day');
%add flag
i=find(strcmp(em_header,'Eruption day fill flag (0=not filled, 1=fill value)'));
flag=cell2mat(em_data(:,i));
ncwriteatt(ncfilename,varname,'Fill flag (0=not filled, 1=fill value)',flag);


% =========================================================================
%% Eruption longitude
%==========================================================================
i=find(strcmp(em_header,'Longitude (-180 to 180)'));
varname='eruption_longitude';

%create new variable
nccreate(ncfilename,varname,'Dimensions',{'eruption',Neru})
%write new variable
ncwrite(ncfilename,varname,cell2mat(em_data(:,i)));
%write variable attribute
ncwriteatt(ncfilename,varname,'Unit','Degree East')
ncwriteatt(ncfilename,varname,'long_name','Eruption longitude');
ncwriteatt(ncfilename,varname,'standard_name','eruption longitude');
%add flag
i=find(strcmp(em_header,'Location (lat/lon/vent altitude) fill flag (0=not filled, 1=fill value)'));
flag=cell2mat(em_data(:,i));
ncwriteatt(ncfilename,varname,'Fill flag (0=not filled, 1=fill value)',flag);


% =========================================================================
%% Eruption latitude
%==========================================================================
i=find(strcmp(em_header,'Latitude (-90 to 90)'));
varname='eruption_latitude';

%create new variable
nccreate(ncfilename,varname,'Dimensions',{'eruption',Neru})
%write new variable
ncwrite(ncfilename,varname,cell2mat(em_data(:,i)));
%write variable attribute
ncwriteatt(ncfilename,varname,'Unit','Degree North')
ncwriteatt(ncfilename,varname,'long_name','Eruption latitude');
ncwriteatt(ncfilename,varname,'standard_name','eruption latitude');
%add flag
i=find(strcmp(em_header,'Location (lat/lon/vent altitude) fill flag (0=not filled, 1=fill value)'));
flag=cell2mat(em_data(:,i));
ncwriteatt(ncfilename,varname,'Fill flag (0=not filled, 1=fill value)',flag);


% =========================================================================
%% Eruption latitude uncertainty
%==========================================================================
i=find(strcmp(em_header,'Latitude uncertainty (degree)'));
varname='eruption_latitude_uncertainty';

%create new variable
nccreate(ncfilename,varname,'Dimensions',{'eruption',Neru})
%write new variable
ncwrite(ncfilename,varname,cell2mat(em_data(:,i)));
%write variable attribute
ncwriteatt(ncfilename,varname,'Unit','Degree, 95% confidence level')
ncwriteatt(ncfilename,varname,'long_name','Eruption latitude uncertainty (95% confidence level)');
ncwriteatt(ncfilename,varname,'standard_name','eruption latitude uncertainty');
%add flag
i=find(strcmp(em_header,'Location (lat/lon/vent altitude) fill flag (0=not filled, 1=fill value)'));
flag=cell2mat(em_data(:,i));
ncwriteatt(ncfilename,varname,'Fill flag (0=not filled, 1=fill value)',flag);



% =========================================================================
%% Sulfur mass
%==========================================================================
i=find(strcmp(em_header,'Sulfur mass (Tg SO2)'));
varname='eruption_sulfur_mass';

%create new variable
nccreate(ncfilename,varname,'Dimensions',{'eruption',Neru})
%write new variable
ncwrite(ncfilename,varname,cell2mat(em_data(:,i)));
%write variable attribute
ncwriteatt(ncfilename,varname,'Unit','Tg SO2')
ncwriteatt(ncfilename,varname,'long_name','Injected sulfur mass');
ncwriteatt(ncfilename,varname,'standard_name','sulfur mass');



% =========================================================================
%% Sulfur mass uncertainty
%==========================================================================
i=find(strcmp(em_header,'Sulfur mass uncertainty (Tg SO2)'));
varname='eruption_sulfur_mass_uncertainty';

%create new variable
nccreate(ncfilename,varname,'Dimensions',{'eruption',Neru})
%write new variable
ncwrite(ncfilename,varname,cell2mat(em_data(:,i)));
%write variable attribute
ncwriteatt(ncfilename,varname,'Unit','Tg SO2, 95% confidence level')
ncwriteatt(ncfilename,varname,'long_name','Injected sulfur mass uncertainty (95% confidence level)');
ncwriteatt(ncfilename,varname,'standard_name','sulfur mass uncertainty');
%add flag
i=find(strcmp(em_header,'Sulfur mass uncertainty fill flag (0=not filled, 1=fill value)'));
flag=cell2mat(em_data(:,i));
ncwriteatt(ncfilename,varname,'Fill flag (0=not filled, 1=fill value)',flag);



% =========================================================================
%% Sulfur injection height
%==========================================================================
i=find(strcmp(em_header,'Injection height (km a.s.l.)'));
varname='eruption_sulfur_injection_height';

%create new variable
nccreate(ncfilename,varname,'Dimensions',{'eruption',Neru})
%write new variable
ncwrite(ncfilename,varname,cell2mat(em_data(:,i)));
%write variable attribute
ncwriteatt(ncfilename,varname,'Unit','km a.s.l.')
ncwriteatt(ncfilename,varname,'long_name','Volcanic sulfur injection height');
ncwriteatt(ncfilename,varname,'standard_name','sulfur injection height');
%add flag
i=find(strcmp(em_header,'Injection height fill flag (0=not filled, 1=fill value)'));
flag=cell2mat(em_data(:,i));
ncwriteatt(ncfilename,varname,'Fill flag (0=not filled, 1=fill value)',flag);

% =========================================================================
%% Sulfur injection height uncertainty
%==========================================================================
i=find(strcmp(em_header,'Injection height uncertainty'));
varname='eruption_sulfur_injection_height_uncertainty';

%create new variable
nccreate(ncfilename,varname,'Dimensions',{'eruption',Neru})
%write new variable
ncwrite(ncfilename,varname,cell2mat(em_data(:,i)));
%write variable attribute
ncwriteatt(ncfilename,varname,'Unit','km (95% confidence level)')
ncwriteatt(ncfilename,varname,'long_name','Volcanic sulfur injection height uncertainty (95% confidence level)');
ncwriteatt(ncfilename,varname,'standard_name','sulfur injection height uncertainty');
%add flag
i=find(strcmp(em_header,'Sulfur mass uncertainty fill flag (0=not filled, 1=fill value)'));
flag=ones(size(cell2mat(em_data(:,i))));
ncwriteatt(ncfilename,varname,'Fill flag (0=not filled, 1=fill value)',flag);


% =========================================================================
%% Sulfur injection depth
%==========================================================================
i=find(strcmp(em_header,'Injection depth (km, 1-sigma of recommended Gaussian vertical injection profile)'));
varname='eruption_sulfur_injection_depth';

%create new variable
nccreate(ncfilename,varname,'Dimensions',{'eruption',Neru})
%write new variable
ncwrite(ncfilename,varname,cell2mat(em_data(:,i)));
%write variable attribute
ncwriteatt(ncfilename,varname,'Unit','km')
ncwriteatt(ncfilename,varname,'long_name','Volcanic sulfur injection depth (1-sigma of recommended Gaussian vertical injection profile)');
ncwriteatt(ncfilename,varname,'standard_name','sulfur injection depth');
%add flag
% i=find(strcmp(em_header,'Injection height fill flag (0=not filled, 1=fill value)'));
% flag=cell2mat(em_data(:,i));
% ncwriteatt(ncfilename,varname,'Fill flag (0=not filled, 1=fill value)',flag);


%==========================================================================
%% Source dataset
%==========================================================================
i=find(strcmp(em_header,'Main source dataset'));
varname='source_dataset';

%create new variable
nccreate(ncfilename,varname,'Dimensions',{'eruption',Neru},'Datatype','string','format','netcdf4')
%write new variable
ncwrite(ncfilename,varname,string(em_data(:,i)));
%write variable attribute
ncwriteatt(ncfilename,varname,'Unit','NA')
ncwriteatt(ncfilename,varname,'long_name','Source dataset for main eruption source parameters');
ncwriteatt(ncfilename,varname,'standard_name','source dataset');


%==========================================================================
%% Information on match for ice-core-derived data
%==========================================================================
i=find(strcmp(em_header,'Information on match for ice-core-derived data'));
varname='information_eruption_match';

%create new variable
nccreate(ncfilename,varname,'Dimensions',{'eruption',Neru},'Datatype','string','format','netcdf4')
%write new variable
ncwrite(ncfilename,varname,string(em_data(:,i)));
%write variable attribute
ncwriteatt(ncfilename,varname,'Unit','NA')
ncwriteatt(ncfilename,varname,'long_name','Information on match for ice-core-derived data');
ncwriteatt(ncfilename,varname,'standard_name','information on eruption match');


%==========================================================================
%% Match confidence
%==========================================================================
i=find(strcmp(em_header,'Match confidence (1=low, 2=medium, 3=high, 4=very high)'));
varname='eruption_match_confidence';

%create new variable
nccreate(ncfilename,varname,'Dimensions',{'eruption',Neru})
%write new variable
ncwrite(ncfilename,varname,cell2mat(em_data(:,i)));
%write variable attribute
ncwriteatt(ncfilename,varname,'Unit','1=low, 2=medium, 3=high')
ncwriteatt(ncfilename,varname,'long_name','Eruption match confidence');
ncwriteatt(ncfilename,varname,'standard_name','eruption match confidence');



%==========================================================================
%% Information on eruption source parameters (date and location)
%==========================================================================
i=find(strcmp(em_header,'Information on eruption source parameters (date, location, height)'));
varname='information_eruption_source_parameters';

%create new variable
nccreate(ncfilename,varname,'Dimensions',{'eruption',Neru},'Datatype','string','format','netcdf4')
%write new variable
ncwrite(ncfilename,varname,string(em_data(:,i)));
%write variable attribute
ncwriteatt(ncfilename,varname,'Unit','NA')
ncwriteatt(ncfilename,varname,'long_name','Information on eruption source parameters not informed in source datasets');
ncwriteatt(ncfilename,varname,'standard_name','information on eruption source parameters');

%==========================================================================
%% Volcanic Explosivity Index
%==========================================================================
i=find(strcmp(em_header,'Volcanic Explosivity Index'));
varname='VEI';

%create new variable
nccreate(ncfilename,varname,'Dimensions',{'eruption',Neru})
%write new variable
ncwrite(ncfilename,varname,cell2mat(em_data(:,i)));
%write variable attribute
ncwriteatt(ncfilename,varname,'Unit','VEI scale')
ncwriteatt(ncfilename,varname,'long_name','Volcanic Explosivity Index');
ncwriteatt(ncfilename,varname,'standard_name','volcanic explosivity index');


%==========================================================================
%% GVP eruption number
%==========================================================================
i=find(strcmp(em_header,'GVP eruption number'));
varname='GVP_eruption_number';

%create new variable
nccreate(ncfilename,varname,'Dimensions',{'eruption',Neru})
%write new variable
ncwrite(ncfilename,varname,cell2mat(em_data(:,i)));
%write variable attribute
ncwriteatt(ncfilename,varname,'Unit','NA')
ncwriteatt(ncfilename,varname,'long_name','Eruption number from the Global Volcanism Program Volcanoes Of The World database');
ncwriteatt(ncfilename,varname,'standard_name','GVP VOTW eruption number');


%==========================================================================
%% GVP volcano number
%==========================================================================
i=find(strcmp(em_header,'GVP volcano number'));
varname='GVP_volcano_number';

%create new variable
nccreate(ncfilename,varname,'Dimensions',{'eruption',Neru})
%write new variable
ncwrite(ncfilename,varname,cell2mat(em_data(:,i)));
%write variable attribute
ncwriteatt(ncfilename,varname,'Unit','NA')
ncwriteatt(ncfilename,varname,'long_name','Volcano number from the Global Volcanism Program Volcanoes Of The World database');
ncwriteatt(ncfilename,varname,'standard_name','GVP VOTW volcano number');


%==========================================================================
%% Vent altitude (km a.s.l.)
%==========================================================================
i=find(strcmp(em_header,'Vent altitude (km a.s.l.)'));
varname='volcano_vent_altitude';

%create new variable
nccreate(ncfilename,varname,'Dimensions',{'eruption',Neru})
%write new variable
ncwrite(ncfilename,varname,cell2mat(em_data(:,i)));
%write variable attribute
ncwriteatt(ncfilename,varname,'Unit','km a.s.l.')
ncwriteatt(ncfilename,varname,'long_name','Volcano vent altitude');
ncwriteatt(ncfilename,varname,'standard_name','volcano vent altitude');



%==========================================================================
%% Ice-core sulfur deposition asymmetry factor (from -1=only Antartica to 1=only Greenland)
%==========================================================================
i=find(strcmp(em_header,'Ice-core sulfur deposition asymmetry factor (from -1=only Antartica to 1=only Greenland)'));
varname='icecore_sulfur_asymmetry';

%create new variable
nccreate(ncfilename,varname,'Dimensions',{'eruption',Neru})
%write new variable
ncwrite(ncfilename,varname,cell2mat(em_data(:,i)));
%write variable attribute
ncwriteatt(ncfilename,varname,'Unit','Unitless ratio, -1=only Antartica to 1=only Greenland')
ncwriteatt(ncfilename,varname,'long_name','Ice-core sulfur deposition asymmetry factor');
ncwriteatt(ncfilename,varname,'standard_name','ice-core sulfur deposition asymmetry');


%==========================================================================
%% Ice year
%==========================================================================
i=find(strcmp(em_header,'Year of deposition in ice-core'));
varname='ice_sulfur_deposition_year';

%create new variable
nccreate(ncfilename,varname,'Dimensions',{'eruption',Neru})
%write new variable
ncwrite(ncfilename,varname,cell2mat(em_data(:,i)));
%write variable attribute
ncwriteatt(ncfilename,varname,'Unit','Year CE')
ncwriteatt(ncfilename,varname,'long_name','Year of sulfur deposition in ice-core');
ncwriteatt(ncfilename,varname,'standard_name','year of sulfur deposition in ice-core');



%et voila!
end

