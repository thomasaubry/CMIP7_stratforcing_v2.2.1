% 	Written by Thomas J. Aubry, July 2025.
% 	Department of Earth Sciences, University of Oxford, UK
%   E-mail: thomas.aubry@earth.ox.ac.uk
% 	Please cite the corresponding paper if you use or modify this script


function cmip7_var_climatology=makeclimatology(cmip7_var,cmip7_month,climatologybounds);
%cmip7_var = variable of interest
%cmip7_month = time dimension in month from 1750
%climatology bounds are >= for first bound and < for second one. e.g. [1850
%2022] would be jan 150-dec 2021

%convert time to fractional year
cmip7_month=(cmip7_month-0.5)/12+1750;
sz=size(cmip7_var);

%check time dimension is first
if sz(1)~=length(cmip7_month);error('time dimension should be first');end

%pre-define climatology matrix
N=length(sz);
cmip7_var_climatology=NaN([12 sz(2:end)]);

%if array is 2D, then calculate climatology in following way
if N==2
    %select period within climatology bound
    cmip7_var=cmip7_var(cmip7_month>=climatologybounds(1) & cmip7_month<climatologybounds(2),:);
    %calculate mean for each month (jan, feb, etc) and save in pre-defined array
    for j=1:12
        cmip7_var_climatology(j)=mean(cmip7_var(j:12:end,:),1,'omitnan');
    end
end

%below I repeat the exact same for different variable size. Very unelegant,
%there must be a better way! oh well...
if N==3
    cmip7_var=cmip7_var(cmip7_month>=climatologybounds(1) & cmip7_month<climatologybounds(2),:,:);
    for j=1:12
        cmip7_var_climatology(j,:,:)=mean(cmip7_var(j:12:end,:,:),1,'omitnan');
    end
end

if N==4
    cmip7_var=cmip7_var(cmip7_month>=climatologybounds(1) & cmip7_month<climatologybounds(2),:,:,:);
    for j=1:12
        cmip7_var_climatology(j,:,:,:)=mean(cmip7_var(j:12:end,:,:,:),1,'omitnan');
    end
end

if N==5
    cmip7_var=cmip7_var(cmip7_month>=climatologybounds(1) & cmip7_month<climatologybounds(2),:,:,:,:);
    for j=1:12
        cmip7_var_climatology(j,:,:,:,:)=mean(cmip7_var(j:12:end,:,:,:,:),1,'omitnan');
    end
end

if N==6
    cmip7_var=cmip7_var(cmip7_month>=climatologybounds(1) & cmip7_month<climatologybounds(2),:,:,:,:,:);
    for j=1:12
        cmip7_var_climatology(j,:,:,:,:,:)=mean(cmip7_var(j:12:end,:,:,:,:,:),1,'omitnan');
    end
end

%et voila!

end

