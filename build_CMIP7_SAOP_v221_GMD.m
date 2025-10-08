% 	Written by Thomas J. Aubry, July 2025.
% 	Department of Earth Sciences, University of Oxford, UK
%   E-mail: thomas.aubry@earth.ox.ac.uk
% 	Please cite the corresponding paper if you use or modify this script



clear
close all

%set dataset version
fversion='2.2.1';



%%=========================================================================
%% 1) Load and pre-process the CMIP7 emission dataset
%%=========================================================================

%Read the emission dataset created from the emission script
%build_CMIP7_emissions_v221_GMD.m
T = readtable('prelim_CMIP7_volcano_S_emissions_v2_2_1.xlsx');

eru_year=T.EruptionYear;
eru_month=T.EruptionMonth_1_12_;
eru_day=T.EruptionDay_1_31_;
eru_lat=T.Latitude__90To90_;
eru_height=T.InjectionHeight_kmA_s_l__;
eru_mass=T.SulfurMass_TgSO2_;%in Tg SO2
eru_source=T.MainSourceDataset;
eru_VEI=T.VolcanicExplosivityIndex;
eru_conf=T.MatchConfidence_1_low_2_medium_3_high_;

%For VEI 5 events with a default mass (i.e. source dataset is GVP), attribute a mass that is the mean of
%all other VEI 5 events
VEI5_default_mass=mean(eru_mass(eru_VEI==5 & ~strcmp(eru_source,'GVP') & eru_mass<10 & eru_conf~=1 & eru_conf~=2));
eru_mass(strcmp(eru_source,'GVP') & eru_VEI==5)=VEI5_default_mass;


%Change latitude of Agung 1963 to -27 to ensure SH-like transport (a dirty
%trick...). Check this changes the latitude of exactly two events.
ind_agung=find(eru_year==1963 & eru_VEI==5);
if length(ind_agung)~=2;error('Agung should only have two events');end
lat_agung=eru_lat(ind_agung);
eru_lat(ind_agung)=-27;


%%=========================================================================
%% 2) Load and pre-process GloSSAC v2.22
%%=========================================================================

gl_ext=varfromname('./source_datasets/GloSSAC_V2.22.nc','Glossac_Aerosol_Extinction_Coefficient',1);
gl_wl=varfromname('./source_datasets/GloSSAC_V2.22.nc','wavelengths_glossac',0);
gl_alt=varfromname('./source_datasets/GloSSAC_V2.22.nc','alt',0);
gl_lat=varfromname('./source_datasets/GloSSAC_V2.22.nc','lat',0);
gl_time=varfromname('./source_datasets/GloSSAC_V2.22.nc','time',0);
gl_tph=varfromname('./source_datasets/GloSSAC_V2.22.nc','trp_hgt',0);
gl_tph=repmat(gl_tph,[length(gl_time)/12 1]);

%set negative values to 0
gl_ext(gl_ext(:)<0)=0;

%set tropospheric values to NaN
for i=1:size(gl_tph,1)
    for j=1:size(gl_tph,2)
        mask=gl_alt<gl_tph(i,j);
        gl_ext(:,mask,j,i)=NaN;
    end
end

%calculate SAOD as sum of extinction*0.5 because vertical grid
%is regularly spaced with 0.5 km spacing
gl_saod=0.5*squeeze(sum(gl_ext,2,'omitnan'));
gl_gmsaod=NaN(size(gl_saod,1),size(gl_saod,3));

%then calculate global mean SAOD
for i=1:size(gl_saod,3)
    for j=1:size(gl_saod,1)
        gl_gmsaod(j,i)=trapz(squeeze(gl_saod(j,:,i))'.*cosd(gl_lat))/trapz(cosd(gl_lat));
    end
end

%reformat time from GloSSAC format to fractionnal year
gl_time=floor(gl_time/100)+((gl_time-100*floor(gl_time/100))-0.5)/12;



%Correct values with high 525nm extinction and low 1020:525 extinction
%ratio. See documentation paper for more details, this correction affects a
%really small portion of the dataset, mostly confined to a small fraction
%of extratropical lowermost stratosphere post-Pinatubo.
ext525=squeeze(gl_ext(3,:,:,:));
ext1020=squeeze(gl_ext(4,:,:,:));
timemat=permute(repmat(gl_time,[1  80 32]),[2 3 1]);

mask=timemat>1991.4 & timemat<1995;

for i=1:size(gl_ext,2)
    for j=1:size(gl_ext,3)
        for k=1:size(gl_ext,4)
            %find values with high 525 ext and low 1020:525 ratio
            if gl_ext(gl_wl==525,i,j,k)>0.0028 & gl_ext(gl_wl==1020,i,j,k)/gl_ext(gl_wl==525,i,j,k)<0.25

                %find data points not needing correction at same
                %time/latitude, and identify the closest one (in terms of
                %altitude)
                heightpos=1:size(gl_ext,2);
                heightpos=heightpos(squeeze(gl_ext(gl_wl==1020,:,j,k)./gl_ext(gl_wl==525,:,j,k))>=0.25);
                [aa clst]=min(abs(heightpos-i));
                clsti=heightpos(clst);

                %replace 1020 ext of corrected data point to match the
                %1020:525 ratio of closest point not requiring correction
                gl_ext(gl_wl==1020,i,j,k)=(gl_ext(gl_wl==1020,clsti,j,k)/gl_ext(gl_wl==525,clsti,j,k))*gl_ext(gl_wl==525,i,j,k);
            end

        end
    end
end


%plot the global mean SAOD anomaly
figure(1)
subplot(2,2,[1 2])
hts1=plot(gl_time,gl_gmsaod(gl_wl==525,:)'-min(gl_gmsaod(gl_wl==525,:)'))

%for the post-Pinatubo period defined as 1998 onwards, find the mean SAOD
%and the max VEI and SO2 mass
refyear=1998;
glmeansaodanom=mean(gl_gmsaod(gl_wl==525,gl_time>=refyear))'-min(gl_gmsaod(gl_wl==525,:)');
maxVEI=max(eru_VEI(eru_year>=refyear));
maxmass=max(eru_mass(eru_year>=refyear));


%Calculate extinction climatology as 1998-2001 mean.
gl_ext_clim=squeeze(gl_ext(gl_wl==525,:,:,gl_time>=1998 & gl_time<2002));
gl_ext_clim=0.25*(gl_ext_clim(:,:,1:12)+gl_ext_clim(:,:,13:24)+gl_ext_clim(:,:,25:36)+gl_ext_clim(:,:,37:48));
%Extend to +/-82.5 and 87.5 degree latitude by repeating +/- 77.5 values
gl_ext_clim=cat(2,gl_ext_clim(:,1,:),gl_ext_clim(:,1,:),gl_ext_clim,gl_ext_clim(:,end,:),gl_ext_clim(:,end,:));
%restrict to altitudes covered by EVA_H
gl_ext_clim=gl_ext_clim(gl_alt>=5 & gl_alt<40,:,:);



%%=========================================================================
%% 3) Define an optimal SO2 mass for VEI 4 eruptions for which source dataset is GVP
%%=========================================================================

%Add EVA_H v2 path and load model parameters
addpath("EVA_H_v2\")
parameterfile


%find all small eruptions pre-satellite, i.e. VEI and mass smaller than
%post-1998 max VEI and 1.5 * max SO2 mass
masksmall=eru_year>=1850 & eru_year<=1978 & eru_VEI<=maxVEI & eru_mass<=3;

%define a list of plausible mass for VEI 4 with unknown masses (0.02-0.5 Tg
%SO2
mass4list=linspace(0.02,0.5,10)';

%Then for each of these masses, calculate the mean SAOD anomaly associated
%with small eruptions pre-satellite era
presatmean=NaN(size(mass4list));
VEImask=strcmp(eru_source,'GVP') & eru_VEI==4; %find small eruptions for which we need to use the default mass, i.e. VEI=4 and source dataset is GVP


for i=1:length(mass4list)
    length(mass4list)-i
    newmass=eru_mass;
    newmass(VEImask)=mass4list(i);%set mass of GVP VEI 4 eruption to tested default mass


    %run EVA_H v2 with the corresponding SO2 masses, for small eruptions
    %and pre-satellite era only
    [inmass intime]=so2injection_8boxes_cmip7(h1lim,h2lim,latlim,newmass(masksmall),eru_height(masksmall),eru_lat(masksmall),eru_year(masksmall),eru_month(masksmall),eru_day(masksmall));%format volcanic emissions
    tspan = ([1850 1978]-1979)*12;%integration bounds in month after jan 1st 1979
    ic = ones(8,1)*10^(-4);%Set initial conditions
    ode_opts = odeset('RelTol',1e-6,'AbsTol',1e-10,'NonNegative',[1:8]); %Set tolerance for the differential equation solver.
    tref=(tspan(1)+0.5:1:tspan(end)-0.5)';%requested output times in month after jan 1st 1979
    time_yr=(tref)/12+1979;%Create a vector containing time in years
    wl_req=[0.525 0.55]';%list of wavelengths at which output are requested, in um

    %run EVA_H
    [t_ode,SO4mass] = ode45(@(t,y) eightboxequations_CMIP7(t,y,inmass,intime,eru_height(masksmall),newmass(masksmall),modelpara), tspan, ic, ode_opts);
    SO4mass=interp1(t_ode,SO4mass,tref,'pchip');%Interpolate outputs to get a monthly timeseries
    %Run the postprocessing function
    [gmsaod, saod, reff, ext, ssa, asy, lat, alt, sad, vd, nd, gm_reff]=postproc_v5(SO4mass,zeros(length(SO4mass),length(gl_lat)+4,length(gl_alt)-10),modelpara,mstar,R_reff,wl_req,'eva_lut_for_SAD_bimodal_Bier75_215_extended.nc');

    %same the 1850-2023 mean SAOD associated with small eruptions for the
    %chosen default VEI 4 mass
    presatmean(i)=mean(gmsaod(:,1));
end


%plot obtained mean pre-satellite small-eruption SAOD as a function of
%default SO2 mass for VEI 4 eruptions with no SO2 mass
figure(1)
subplot(2,2,3)
h1=plot(mass4list,presatmean,'ko')
hold on
plot([0 1],[1 1]*glmeansaodanom,'--','LineWidth',2,'Color','#0072BD')
xlabel('Assumed VEI 4 mass for GVP-derived pre-satellite eruptions (Tg SO2)')
ylabel({'1850-1978 mean SAOD anomaly';'for small-magnitude eruptions'})
meanVEI4=mean(eru_mass(eru_year>=1979 & eru_VEI==4));
plot([meanVEI4 meanVEI4],[2.6 4.6]/1000,'k:','LineWidth',1.5)
text(0.32,2.95/1000,'GloSSAC 1998-2023 mean','Color','#0072BD')
vei4mlist=eru_mass(strcmp(eru_source,'MSVOLSO2L4') & eru_VEI==4);
text(0.25,3.65/1000,{'MSVOLSO2L4 mean';strcat('VEI4 mass:',num2str(round(mean(vei4mlist),2)),'Tg SO_2');strcat('(min-max:',num2str(round(min(vei4mlist),4)),'-',num2str(round(max(vei4mlist),2)),')')},'Rotation',90)
ylim([2.8 4.6]/1000)
xlim([0 0.52])

%define the optimal VEI 4 mass as the mass enabling to match the
%pre-satellite small-eruption only SAOD to the observed 1998-2023 mean SAOD
%(only small eruptions for that time period)
optVEI4mass=interp1(presatmean,mass4list,glmeansaodanom,'pchip');
h2=plot(optVEI4mass,glmeansaodanom,'kp','MarkerSize',15,'MarkerFaceColor','k')
legend([h1 h2],'Test EVA\_H runs','Optimal VEI4 mass','Location','East')



%%=========================================================================
%% 4) Replace default VEI masses in emission dataset
%%=========================================================================

%do the actual replacing
eru_mass(strcmp(eru_source,'GVP') & eru_VEI==4)=optVEI4mass;


%run EVA_H for small eruptions only pre-satellite with chosen mass.

[inmass intime]=so2injection_8boxes_cmip7(h1lim,h2lim,latlim,eru_mass(masksmall),eru_height(masksmall),eru_lat(masksmall),eru_year(masksmall),eru_month(masksmall),eru_day(masksmall));
[t_ode,SO4mass] = ode45(@(t,y) eightboxequations_CMIP7(t,y,inmass,intime,eru_height(masksmall),eru_mass(masksmall),modelpara), tspan, ic, ode_opts);
SO4mass=interp1(t_ode,SO4mass,tref,'pchip');
[gmsaod, saod, reff, ext, ssa, asy, lat, alt, sad, vd, nd, gm_reff]=postproc_v5(SO4mass,zeros(length(SO4mass),length(gl_lat)+4,length(gl_alt)-10),modelpara,mstar,R_reff,wl_req,'eva_lut_for_SAD_bimodal_Bier75_215_extended.nc');



%add SAOD time series to previous plot
gmsaod=gmsaod(:,1);
figure(1)
subplot(2,2,[1 2])
hold on
hts2=plot(time_yr,gmsaod)
hts3=plot([1850 2024],glmeansaodanom*[1 1],'--','Color','#0072BD','LineWidth',2)
title('Pre-satellite global mean SAOD anomaly at 525nm, small-magnitude eruptions only')
ylabel('SAOD anomaly')
ylim([0 0.035])
xlim([1850 2024])
ylabel('SAOD anomaly')
plot([refyear refyear],[0 0.035],'k:','LineWidth',2)
text(refyear+0.5,0.017,'For 1998-2023:')
text(refyear+0.5,0.013,{'Max SO_2 mass=';strcat(num2str(round(maxmass,5)),'Tg')})
text(refyear+0.5,0.009,{'Max VEI=';num2str(maxVEI)})
legend([hts1 hts3 hts2],'GloSSAC (1979-2023)','GloSSAC 1998-2023 mean',strcat('1850-1978, <= 3 Tg SO_2 eruptions only, default VEI 4 mass = ',num2str(round(optVEI4mass,2)),' Tg SO_2'),'Location','north')


%add a subplot showing probability distribution of SAOD anomalies for
%pre-satellite era small eruptions and for satellite era post 1998. Note
%that although the mean is equal, the distribution is still biased pre-satellite
%era with not enough very small peak and too many moderate peak

figure(1)
bins=[0:0.001:0.012 0.014:0.002:0.02];
subplot(2,2,4)
h1=histogram(gmsaod,bins,'Normalization','probability')
hold on
h2=histogram(gl_gmsaod(gl_wl==525,gl_time>=1998)'-min(gl_gmsaod(gl_wl==525,:)'),bins,'Normalization','probability')
xlabel('Monthly global mean SAOD anomaly (525nm)')
ylabel('Probability')
legend('1850-1978, <=3 Tg SO_2 eruptions only','GloSSAC 1998-2023')
xlim([min(bins) max(bins)])
ylim([0 0.52])


%%=========================================================================
%% 5) Save final emission dataset
%%=========================================================================

%set file name
emisfilename=strcat('utsvolcsulfur_input4MIPs_emissions_CMIP7_',fversion,'_grid_label_1750-2023');

%change latitude of Agung back to actual latitude; the paper documents
%clearly this deviation from the emission dataset when running EVA_H v2

ind_agung=find(eru_year==1963 & eru_VEI==5);
if length(ind_agung)~=2;error('Agung should only have two events');end
eru_lat(ind_agung)=lat_agung;

%update sulfur mass (we changed mass of VEI 4 eruptions for which source
%dataset is GVP) and latitude (we changed Agung latitude back to real one)
T.SulfurMass_TgSO2_=eru_mass;
T.Latitude__90To90_=eru_lat;

%save as an excel table
C = table2cell(T);
addpath("EVA_H_v2\")
cmip7_all_headers={'Volcano name',...
    'Eruption year',...
    'Eruption month (1-12)',...
    'Eruption month fill flag (0=not filled, 1=fill value)',...
    'Eruption day (1-31)',...
    'Eruption day fill flag (0=not filled, 1=fill value)',...
    'Longitude (-180 to 180)',...
    'Latitude (-90 to 90)',...
    'Latitude uncertainty (degree)',...
    'Location (lat/lon/vent altitude) fill flag (0=not filled, 1=fill value)',...
    'Sulfur mass (Tg SO2)',...
    'Sulfur mass uncertainty (Tg SO2)',...
    'Sulfur mass uncertainty fill flag (0=not filled, 1=fill value)',...
    'Injection height (km a.s.l.)',...
    'Injection height fill flag (0=not filled, 1=fill value)',...
    'Injection height uncertainty',...
    'Injection depth (km, 1-sigma of recommended Gaussian vertical injection profile)',...
    'Main source dataset',...
    'Information on match for ice-core-derived data',...
    'Match confidence (1=low, 2=medium, 3=high, 4=very high)',...
    'Information on eruption source parameters (date, location, height)',...
    'Volcanic Explosivity Index',...
    'GVP eruption number',...
    'GVP volcano number',...
    'Vent altitude (km a.s.l.)',...
    'Ice-core sulfur deposition asymmetry factor (from -1=only Antartica to 1=only Greenland)',...
    'Year of deposition in ice-core'};
cmip7_all_final=[cmip7_all_headers;C];
writecell(cmip7_all_final,strcat(emisfilename,'.xlsx'),'FileType','spreadsheet','WriteMode','replacefile')


%save as a netcdf using the bespoke emistonetcdf.m function
emistonetcdf(strcat(emisfilename,'.xlsx'),strcat(emisfilename,'.nc'));


%change latitude of Agung back to -27 to run EVA_H v2 again later down
ind_agung=find(eru_year==1963 & eru_VEI==5);
if length(ind_agung)~=2;error('Agung should only have two events');end
eru_lat(ind_agung)=-27;

%%=========================================================================
%% 6) Run EVA_H vs with the final CMIP7 inventory with modified Agung latitude
%%=========================================================================

%load EVA_H v2 parameter
parameterfile

%calculate volcanic injections with full CMIP7 inventory
[inmass intime]=so2injection_8boxes_cmip7(h1lim,h2lim,latlim,eru_mass,eru_height,eru_lat,eru_year,eru_month,eru_day);

%Set integration parameters
tspan = ([1750 2024]-1979)*12;%integration bounds in month after jan 1st gl_start_year
ic = ones(8,1)*10^(-4);%initial conditions
ode_opts = odeset('RelTol',1e-6,'AbsTol',1e-10,'NonNegative',[1:8]);%Set tolerance for the differential equation solver.
tref=(tspan(1)+0.5:1:tspan(end)-0.5)';%time to which output will be interpolated
time_yr=(tref)/12+1979;%Create a vector containing time in years

%run EVA_H and interpolate output
[t_ode,SO4mass] = ode45(@(t,y) eightboxequations_CMIP7(t,y,inmass,intime,eru_height,eru_mass,modelpara), tspan, ic, ode_opts);
SO4mass=interp1(t_ode,SO4mass,tref,'pchip');


%repeat/permute GloSSAC derived non-volcanic background aerosol climatology
%to have the same dimension as EVA_H extinction
gl_ext_clim=repmat(gl_ext_clim,[1 1 size(SO4mass,1)/12]);
gl_ext_clim=permute(gl_ext_clim,[3 2 1]);

%add a linear trend to the background climatology between 1850 and 2000,
%assuming proportionality to anthropogenic sulfur emissions. See companion
%paper for more.
anthremis1850=73.6;
anthremis2000=109;
anthremis=anthremis1850+(anthremis2000-anthremis1850)*(time_yr-1850)/(2000-1850);
anthremis(anthremis>anthremis2000)=anthremis2000;
anthremis(anthremis<anthremis1850)=anthremis1850;
gl_ext_clim=gl_ext_clim.*repmat(anthremis/max(anthremis),[1 size(gl_ext_clim,2) size(gl_ext_clim,3)]);

%list of wavelengths at which output are requested, in um
%List follows recommendations by Ben Johnson and Anton Laakso
%Mie file covers 0.15-100; don't request wavelengths outside that range
wl_req=[0.16 0.23 0.3 0.39 0.46 0.525 0.53 0.55 0.61 0.7 0.8 0.9 1.01 1.02 1.27 1.46 1.78 2.05 2.33 2.79 3.418 4.016 4.319 4.618 5.154 6.097 6.8 7.782 8.02 8.849 9.708 11.111 13.157 15.037 17.699 20.0 23.529 35 50 75 100];

%Run the postprocessing function; not the first input is the volcanic
%sulfate aerosol mass time series, and the second input is the non-volcanic
%background extinction. The latter was zeroed in previous EVA_H
%runs, but we now add the climatology with a trend derived above.
[gmsaod, saod, reff, ext, ssa, asy, lat, alt, sad, vd, nd]=postproc_v5(SO4mass,gl_ext_clim,modelpara,mstar,R_reff,wl_req,'eva_lut_for_SAD_bimodal_Bier75_215_extended.nc');

%find indiced of 525 and 1020 nm which are the key wavelengths of the
%dataset
i525=find(wl_req==0.525);
i1020=find(wl_req==1.020);


%%=========================================================================
%% 7) Harmonize GloSSAC and EVA_H v2 product over 1979-1981 and get full properties for GloSSAC
%%=========================================================================

gl_start_year=floor(min(gl_time));%GloSSAC start year

%a bit of pre-processing for GloSSAC:
gl_ext=cat(4,squeeze(gl_ext(gl_wl==525,:,:,:)),squeeze(gl_ext(gl_wl==1020,:,:,:)));%select only 525 and 1020nm
gl_ext=permute(gl_ext,[3 2 1 4]);%reorganize dimension to be same as EVA_H extinction
gl_ext=gl_ext(:,:,gl_alt<=max(alt) & gl_alt>=min(alt),:);%limit altitude range to EVA_H altitude (5-39.5 km)
gl_ext=cat(2,repmat(gl_ext(:,1,:,:),[1 sum(lat<(min(gl_lat))) 1 1]),gl_ext);%repeat 77.5S extinction at 82.5S and 87.5S
gl_ext=cat(2,gl_ext,repmat(gl_ext(:,end,:,:),[1 sum(lat>(max(gl_lat))) 1 1]));%repeat 77.5N extinction at 82.5N and 87.5N

% define length of harmonization period
yrsmt=3;

%initialize a merged extinction matrix mg_ext
mg_ext=gl_ext;

%define weight given to GloSSAC during harmonization: linearly ramping from
%0 t start of harmonization period to 1 3 years later, then constant
wgt=linspace(0,1,sum(gl_time>=gl_start_year & gl_time<=gl_start_year+yrsmt))';
wgt=repmat(wgt,[1 size(mg_ext,2:4)]);

%for harmonization period, merged extinction is:
% weight * (GloSSAC extinction) + (1-weight) * (EVA_H extinction)
mg_ext(gl_time<=gl_start_year+yrsmt & gl_time>=gl_start_year,:,:,:)=wgt.*gl_ext(gl_time<=gl_start_year+yrsmt & gl_time>=gl_start_year,:,:,:) + (1-wgt).*ext(time_yr<=gl_start_year+yrsmt & time_yr>=gl_start_year,:,:,[i525 i1020]);

%the latitude and altitude of the merged product are as EVA_H lat/alt
mg_alt=alt;mg_lat=lat;

%From the merged extinction at 525 and 1020, get full aerosol optical
%properties using EVA_H Mie routine
[mg_gmsaod, mg_saod, mg_reff, mg_ext, mg_ssa, mg_asy, mg_sad, mg_vd, mg_nd]=postproc_glwl_v3(squeeze(mg_ext(:,:,:,1)),0.525,squeeze(mg_ext(:,:,:,2)),1.02,wl_req,mg_lat,'eva_lut_for_SAD_bimodal_Bier75_215_extended.nc');


%%=========================================================================
%% 8) Combine the EVA_H and GloSSAC products
%%=========================================================================

%In section 7, we harmonized for the satellite era. We now simply combine
%the pre-satellite (emission+EVA_H v2 derived) and satellite era products

%Define a flag for data source
mask_evah=time_yr<gl_start_year;
cmip7_time=[time_yr(mask_evah);gl_time];
flag_source=ones(size(cmip7_time));%1=derived from emission+EVA_H
flag_source(cmip7_time>=gl_start_year & cmip7_time<=gl_start_year+yrsmt)=2;%2=harmonization period
flag_source(cmip7_time>gl_start_year+yrsmt)=3;%3 = GloSSAC

%altitude, latitude and wavelength
cmip7_alt=alt;
cmip7_lat=mg_lat;
cmip7_wl=wl_req;

%for all aerosol optical properties, simply concatenate EVA_H and merged
%GloSSAC product along time dimension
cmip7_gmsaod=cat(1,gmsaod(mask_evah,:),mg_gmsaod);
cmip7_saod=cat(1,saod(mask_evah,:,:),mg_saod);
cmip7_reff=cat(1,reff(mask_evah,:,:),mg_reff);
cmip7_ext=cat(1,ext(mask_evah,:,:,:),mg_ext);
cmip7_ssa=cat(1,ssa(mask_evah,:,:,:),mg_ssa);
cmip7_asy=cat(1,asy(mask_evah,:,:,:),mg_asy);
cmip7_sad=cat(1,sad(mask_evah,:,:),mg_sad);
cmip7_vd=cat(1,vd(mask_evah,:,:),mg_vd);
cmip7_nd=cat(1,nd(mask_evah,:,:),mg_nd);

%find indices of 525 and 1020 nm wavelengths, it will be handy
i525=find(wl_req==0.525);
i1020=find(wl_req==1.020);
i550=find(cmip7_wl==0.55);



%%=========================================================================
%% 9) Make some basic plots comparing CMIP6 and CMIP7
%%=========================================================================

%load CMIP6 mid-visible extinction
cmip6_ext=varfromname('./CMIP6/CMIP_1850_2014_extinction_550nm_strat_only_v3.nc','ext550',-999);
cmip6_lat=varfromname('./CMIP6/CMIP_1850_2014_extinction_550nm_strat_only_v3.nc','latitude',-999);
cmip6_time=varfromname('./CMIP6/CMIP_1850_2014_extinction_550nm_strat_only_v3.nc','month',-999);


%calculate SAOD by multiplying by 0.5 km the sum of extinction (as
%regularly spaced vertically on 0.5 km grid)
cmip6_saod=squeeze(sum(cmip6_ext,2))*0.5;
cmip6_time=(cmip6_time-0.5)/12+1850;

%calculate global mean
cmip6_gmsaod=NaN(size(cmip6_saod,1),1);
for i=1:size(cmip6_saod,1)
    cmip6_gmsaod(i)=trapz(squeeze(cmip6_saod(i,:))'.*cosd(cmip6_lat))/trapz(cosd(cmip6_lat));
end


%plot latitudinal SAOD distribution
lims=[1850 1895;1895 1940;1940 1985;1985 2023];
figure(3)
colormap(flipud(pink))
for i=1:size(lims,1)

    subplot(4,2,2*i-1)
    aa=pcolor(cmip7_time,cmip7_lat,log10(squeeze(cmip7_saod(:,:,i550))'))
    set(aa,'EdgeColor','none')
    xlim([lims(i,1) lims(i,2)])
    caxis([-2.7 -0.15])
    ylabel('Latitude')
end

for i=1:size(lims,1)

    subplot(4,2,2*i)
    aa=pcolor(cmip6_time,cmip6_lat,log10(cmip6_saod'))
    set(aa,'EdgeColor','none')
    xlim([lims(i,1) lims(i,2)])
    caxis([-2.7 -0.15])
    ylabel('Latitude')
end

subplot(4,2,1)
title('CMIP7')
subplot(4,2,2)
title('CMIP6')

%plot global means
figure(2)
plot(cmip7_time,cmip7_gmsaod(:,i550),'r-','LineWidth',2)
hold on
plot(cmip6_time,cmip6_gmsaod(:,1),'b-','LineWidth',2)
xlim([1750 2023])
ylabel('Global mean SAOD (550 nm)')
legend('CMIP7','CMIP6')



% =========================================================================
%% 10) Write in netcdf
%==========================================================================

%NOTE: All netcdf written in this script (including the emission netcdf
%above) are not fully formatted to CMIP standard (e.g. CF-compliance, CMIP
%compliance, dimension orders, unit, etc). An additional python script,
%originally developped by Zebedee Nicholls, was used to reformat the netcdf
%produced by this script into the ones provided to the community. These
%python scripts can be provided upon request.


%define file name and delete any previous version
ncfilename=strcat('strat_aer_opt_prop_input4MIPs_type_CMIP7_',fversion,'_grid_label_1750-2023.nc');
delete ncfilename
ncid = netcdf.create(ncfilename,'NETCDF4');
netcdf.close(ncid);

%write global attribute

ncwriteatt(ncfilename,'/','activity_id','input4MIPs')
ncwriteatt(ncfilename,'/','activity_id','input4MIPs')
ncwriteatt(ncfilename,'/','contact','CMIP7 Forcing Task Team: volcano dataset providers (thom.aubry@gmail.com,Vaishali.Naik@noaa.gov,durack1@llnl.gov)')
ncwriteatt(ncfilename,'/','creation_date',strcat(datestr(now, 'yyyy-mm-dd'),'T',datestr(now, 'HH:MM:SS'),'Z'))
ncwriteatt(ncfilename,'/','data_specs_version','01.00')
ncwriteatt(ncfilename,'/','dataset_category','Stratospheric aerosol optical properties')
ncwriteatt(ncfilename,'/','frequency','Monthly')
ncwriteatt(ncfilename,'/','further_info_url','Link to documentation to be added')
ncwriteatt(ncfilename,'/','grid','Ill-defined for this; ask advise to CFTT leadership')
ncwriteatt(ncfilename,'/','grid_label','Ill-defined for this; ask advise to CFTT leadership')
ncwriteatt(ncfilename,'/','history','')
ncwriteatt(ncfilename,'/','institution','tbd after further discussion with Paul')
ncwriteatt(ncfilename,'/','institution_id','tbd after further discussion with Paul')
ncwriteatt(ncfilename,'/','mip_era','CMIP6plus')
ncwriteatt(ncfilename,'/','release_year','2024')
ncwriteatt(ncfilename,'/','source','MSVOLSO2L4, eVolv2k, Global Volcanism Programme VOTW, IVESPA, and other sources as specified in attributes and documentation')
ncwriteatt(ncfilename,'/','source_id','tbd after further discussion with Paul')
ncwriteatt(ncfilename,'/','source_type','satellite, model, ice-core and geological records')
ncwriteatt(ncfilename,'/','target_mip','CMIP6Plus')

%define dimensions lengths
cmip7_month=1:length(cmip7_time);
Nmon=length(cmip7_time);
Nlat=length(cmip7_lat);
Nalt=length(cmip7_alt);
Nwl=length(cmip7_wl);

%below, save each dimension/variable in the file

%==========================================================================
% Time dimension
varname='month';
%create new variable
nccreate(ncfilename,varname,'Dimensions',{'month',Nmon})
%write new variable
ncwrite(ncfilename,varname,cmip7_month);
%write variable attribute
ncwriteatt(ncfilename,varname,'Unit','month_from_jan_1750')
ncwriteatt(ncfilename,varname,'long_name','Month from January 1750');
ncwriteatt(ncfilename,varname,'standard_name','Month from January 1750');

%==========================================================================
% Latitude dimension
varname='latitude';
%create new variable
nccreate(ncfilename,varname,'Dimensions',{'latitude',Nlat})
%write new variable
ncwrite(ncfilename,varname,cmip7_lat);
%write variable attribute
ncwriteatt(ncfilename,varname,'Unit','degrees_north')
ncwriteatt(ncfilename,varname,'long_name','Latitude');
ncwriteatt(ncfilename,varname,'standard_name','latitude');

%==========================================================================
% Altitude dimension
varname='altitude';
%create new variable
nccreate(ncfilename,varname,'Dimensions',{'altitude',Nalt})
%write new variable
ncwrite(ncfilename,varname,cmip7_alt);
%write variable attribute
ncwriteatt(ncfilename,varname,'Unit','km_asl')
ncwriteatt(ncfilename,varname,'long_name','Altitude');
ncwriteatt(ncfilename,varname,'standard_name','altitude');


%==========================================================================
% Wavelength dimension
varname='wavelength';
%create new variable
nccreate(ncfilename,varname,'Dimensions',{'wavelength',Nwl})
%write new variable
ncwrite(ncfilename,varname,cmip7_wl);
%write variable attribute
ncwriteatt(ncfilename,varname,'Unit','micrometer')
ncwriteatt(ncfilename,varname,'long_name','Wavelength');
ncwriteatt(ncfilename,varname,'standard_name','wavelength');


%==========================================================================
% Data source
varname='data_source';

%create new variable
nccreate(ncfilename,varname,'Dimensions',{'month',Nmon})
%write new variable
ncwrite(ncfilename,varname,flag_source);
%write variable attribute
ncwriteatt(ncfilename,varname,'long_name','Data source: 1=from CMIP7 volcanic sulfur emission using aerosol model, 3=GloSSAC, 2=mix of 1 and 3');
ncwriteatt(ncfilename,varname,'standard_name','data source');



%==========================================================================
% extinction
varname='ext';

%create new variable
nccreate(ncfilename,varname,'Dimensions',{'month',Nmon,'latitude',Nlat,'altitude',Nalt,'wavelength',Nwl})
%write new variable
ncwrite(ncfilename,varname,cmip7_ext);
%write variable attribute
ncwriteatt(ncfilename,varname,'Unit','1/km')
ncwriteatt(ncfilename,varname,'long_name','Extinction');
ncwriteatt(ncfilename,varname,'standard_name','Extinction');


%==========================================================================
% scattering asymmtery factor
varname='asy';

%create new variable
nccreate(ncfilename,varname,'Dimensions',{'month',Nmon,'latitude',Nlat,'altitude',Nalt,'wavelength',Nwl})
%write new variable
ncwrite(ncfilename,varname,cmip7_asy);
%write variable attribute
ncwriteatt(ncfilename,varname,'Unit','unitless')
ncwriteatt(ncfilename,varname,'long_name','Scattering Asymmtery Factor');
ncwriteatt(ncfilename,varname,'standard_name','scattering asymmtery factor');


%==========================================================================
% single scattering albedo
varname='ssa';

%create new variable
nccreate(ncfilename,varname,'Dimensions',{'month',Nmon,'latitude',Nlat,'altitude',Nalt,'wavelength',Nwl})
%write new variable
ncwrite(ncfilename,varname,cmip7_ssa);
%write variable attribute
ncwriteatt(ncfilename,varname,'Unit','unitless')
ncwriteatt(ncfilename,varname,'long_name','Single Scattering Albedo');
ncwriteatt(ncfilename,varname,'standard_name','single scattering albedo');


%==========================================================================
% effective radius
varname='reff';

%create new variable
nccreate(ncfilename,varname,'Dimensions',{'month',Nmon,'latitude',Nlat,'altitude',Nalt})
%write new variable
ncwrite(ncfilename,varname,cmip7_reff);
%write variable attribute
ncwriteatt(ncfilename,varname,'Unit','micrometers')
ncwriteatt(ncfilename,varname,'long_name','Effective Radius');
ncwriteatt(ncfilename,varname,'standard_name','effective radius');


%==========================================================================
% surface area density
varname='sad';

%create new variable
nccreate(ncfilename,varname,'Dimensions',{'month',Nmon,'latitude',Nlat,'altitude',Nalt})
%write new variable
ncwrite(ncfilename,varname,cmip7_sad);
%write variable attribute
ncwriteatt(ncfilename,varname,'Unit','um^2/cm^3')
ncwriteatt(ncfilename,varname,'long_name','Surface Area Density');
ncwriteatt(ncfilename,varname,'standard_name','surface area density');


%==========================================================================
% volume density
varname='vd';

%create new variable
nccreate(ncfilename,varname,'Dimensions',{'month',Nmon,'latitude',Nlat,'altitude',Nalt})
%write new variable
ncwrite(ncfilename,varname,cmip7_vd);
%write variable attribute
ncwriteatt(ncfilename,varname,'Unit','um^3/cm^3')
ncwriteatt(ncfilename,varname,'long_name','Volume Density');
ncwriteatt(ncfilename,varname,'standard_name','volume density');

%==========================================================================
% number density
varname='nd';

%create new variable
nccreate(ncfilename,varname,'Dimensions',{'month',Nmon,'latitude',Nlat,'altitude',Nalt})
%write new variable
ncwrite(ncfilename,varname,cmip7_nd);
%write variable attribute
ncwriteatt(ncfilename,varname,'Unit','H2SO4 molecules/cm^3')
ncwriteatt(ncfilename,varname,'long_name','Aerosol Number Density');
ncwriteatt(ncfilename,varname,'standard_name','Number density');



%==========================================================================
% saod
varname='saod';

%create new variable
nccreate(ncfilename,varname,'Dimensions',{'month',Nmon,'latitude',Nlat,'wavelength',Nwl})
%write new variable
ncwrite(ncfilename,varname,cmip7_saod);
%write variable attribute
ncwriteatt(ncfilename,varname,'Unit','unitless')
ncwriteatt(ncfilename,varname,'long_name','Stratospheric Aerosol Optical Depth');
ncwriteatt(ncfilename,varname,'standard_name','stratospheric aerosol optical depth');

%==========================================================================
% global mean saod 550nm
varname='gm_saod550';

%create new variable
nccreate(ncfilename,varname,'Dimensions',{'month',Nmon})
%write new variable
ncwrite(ncfilename,varname,cmip7_gmsaod(:,i550));
%write variable attribute
ncwriteatt(ncfilename,varname,'Unit','unitless')
ncwriteatt(ncfilename,varname,'long_name','Global mean stratospheric Aerosol Optical Depth at 550 nm');
ncwriteatt(ncfilename,varname,'standard_name','global mean stratospheric aerosol optical depth at 550 nm');


% =========================================================================
%% 11) Calculate preindustrial climatologies and save in netcdf too
%==========================================================================

%bound is >= and < so 2022 means end in dec 2021
climatologybounds=[1850 2022];
Nmonclim=12;

%for each variable below, calculate climatology using the bespoke function
% makeclimatology.m, and save in the netcdf


%==========================================================================
% Time dimension climatology
varname='month_climatology';
%create new variable
nccreate(ncfilename,varname,'Dimensions',{'month_climatology',Nmonclim})
%write new variable
ncwrite(ncfilename,varname,cmip7_month(1:12));
%write variable attribute
ncwriteatt(ncfilename,varname,'Unit','month_from_jan_1850')
ncwriteatt(ncfilename,varname,'long_name','Month from January 1850');
ncwriteatt(ncfilename,varname,'standard_name','Month from January 1850');


%==========================================================================
% extinction climatology
varname='ext_climatology';
cmip7_ext_climatology=makeclimatology(cmip7_ext,cmip7_month,climatologybounds);
%create new variable
nccreate(ncfilename,varname,'Dimensions',{'month_climatology',Nmonclim,'latitude',Nlat,'altitude',Nalt,'wavelength',Nwl})
%write new variable
ncwrite(ncfilename,varname,cmip7_ext_climatology);
%write variable attribute
ncwriteatt(ncfilename,varname,'Unit','1/km')
ncwriteatt(ncfilename,varname,'long_name','1850-2021 climatology of Extinction');
ncwriteatt(ncfilename,varname,'standard_name','1850-2021 climatology of Extinction');


%==========================================================================
% scattering asymmtery factor climatology
varname='asy_climatology';
cmip7_asy_climatology=makeclimatology(cmip7_asy,cmip7_month,climatologybounds);
%create new variable
nccreate(ncfilename,varname,'Dimensions',{'month_climatology',Nmonclim,'latitude',Nlat,'altitude',Nalt,'wavelength',Nwl})
%write new variable
ncwrite(ncfilename,varname,cmip7_asy_climatology);
%write variable attribute
ncwriteatt(ncfilename,varname,'Unit','unitless')
ncwriteatt(ncfilename,varname,'long_name','1850-2021 climatology of Scattering Asymmtery Factor');
ncwriteatt(ncfilename,varname,'standard_name','1850-2021 climatology of scattering asymmtery factor');


%==========================================================================
% single scattering albedo climatology
varname='ssa_climatology';
cmip7_ssa_climatology=makeclimatology(cmip7_ssa,cmip7_month,climatologybounds);
%create new variable
nccreate(ncfilename,varname,'Dimensions',{'month_climatology',Nmonclim,'latitude',Nlat,'altitude',Nalt,'wavelength',Nwl})
%write new variable
ncwrite(ncfilename,varname,cmip7_ssa_climatology);
%write variable attribute
ncwriteatt(ncfilename,varname,'Unit','unitless')
ncwriteatt(ncfilename,varname,'long_name','1850-2021 climatology of Single Scattering Albedo');
ncwriteatt(ncfilename,varname,'standard_name','1850-2021 climatology of single scattering albedo');


%==========================================================================
% effective radius climatology
varname='reff_climatology';
cmip7_reff_climatology=makeclimatology(cmip7_reff,cmip7_month,climatologybounds);
%create new variable
nccreate(ncfilename,varname,'Dimensions',{'month_climatology',Nmonclim,'latitude',Nlat,'altitude',Nalt})
%write new variable
ncwrite(ncfilename,varname,cmip7_reff_climatology);
%write variable attribute
ncwriteatt(ncfilename,varname,'Unit','micrometers')
ncwriteatt(ncfilename,varname,'long_name','1850-2021 climatology of Effective Radius');
ncwriteatt(ncfilename,varname,'standard_name','1850-2021 climatology of effective radius');


%==========================================================================
% surface area density climatology
varname='sad_climatology';
cmip7_sad_climatology=makeclimatology(cmip7_sad,cmip7_month,climatologybounds);
%create new variable
nccreate(ncfilename,varname,'Dimensions',{'month_climatology',Nmonclim,'latitude',Nlat,'altitude',Nalt})
%write new variable
ncwrite(ncfilename,varname,cmip7_sad_climatology);
%write variable attribute
ncwriteatt(ncfilename,varname,'Unit','um^2/cm^3')
ncwriteatt(ncfilename,varname,'long_name','1850-2021 climatology of Surface Area Density');
ncwriteatt(ncfilename,varname,'standard_name','1850-2021 climatology of surface area density');


%==========================================================================
% volume density climatology
varname='vd_climatology';
cmip7_vd_climatology=makeclimatology(cmip7_vd,cmip7_month,climatologybounds);
%create new variable
nccreate(ncfilename,varname,'Dimensions',{'month_climatology',Nmonclim,'latitude',Nlat,'altitude',Nalt})
%write new variable
ncwrite(ncfilename,varname,cmip7_vd_climatology);
%write variable attribute
ncwriteatt(ncfilename,varname,'Unit','um^3/cm^3')
ncwriteatt(ncfilename,varname,'long_name','1850-2021 climatology of Volume Density');
ncwriteatt(ncfilename,varname,'standard_name','1850-2021 climatology of volume density');

%==========================================================================
% volume density climatology
varname='nd_climatology';
cmip7_nd_climatology=makeclimatology(cmip7_nd,cmip7_month,climatologybounds);
%create new variable
nccreate(ncfilename,varname,'Dimensions',{'month_climatology',Nmonclim,'latitude',Nlat,'altitude',Nalt})
%write new variable
ncwrite(ncfilename,varname,cmip7_nd_climatology);
%write variable attribute
ncwriteatt(ncfilename,varname,'Unit','H2SO4 molecules/cm^3')
ncwriteatt(ncfilename,varname,'long_name','1850-2021 climatology of Aerosol Number Density');
ncwriteatt(ncfilename,varname,'standard_name','1850-2021 climatology of number density');

%==========================================================================
% saod climatology
varname='saod_climatology';
cmip7_saod_climatology=makeclimatology(cmip7_saod,cmip7_month,climatologybounds);
%create new variable
nccreate(ncfilename,varname,'Dimensions',{'month_climatology',Nmonclim,'latitude',Nlat,'wavelength',Nwl})
%write new variable
ncwrite(ncfilename,varname,cmip7_saod_climatology);
%write variable attribute
ncwriteatt(ncfilename,varname,'Unit','unitless')
ncwriteatt(ncfilename,varname,'long_name','1850-2021 climatology of Stratospheric Aerosol Optical Depth');
ncwriteatt(ncfilename,varname,'standard_name','1850-2021 climatology of stratospheric aerosol optical depth');

%==========================================================================
% global mean saod 550nm climatology
varname='gm_saod550_climatology';

cmip7_gmsaod550_climatology=makeclimatology(cmip7_gmsaod(:,i550),cmip7_month,climatologybounds);
%create new variable
nccreate(ncfilename,varname,'Dimensions',{'month_climatology',Nmonclim})
%write new variable
ncwrite(ncfilename,varname,cmip7_gmsaod550_climatology);
%write variable attribute
ncwriteatt(ncfilename,varname,'Unit','unitless')
ncwriteatt(ncfilename,varname,'long_name','1850-2021 climatology of Global mean stratospheric Aerosol Optical Depth at 550 nm');
ncwriteatt(ncfilename,varname,'standard_name','1850-2021 climatology of global mean stratospheric aerosol optical depth at 550 nm');


%add climatology to global mean SAOD plot and add legend
figure(2)
plot(cmip7_time,repmat(cmip7_gmsaod550_climatology,[length(cmip7_time)/12 1]),'r:')
legend('CMIP7','CMIP6','CMIP7 climatology')


%and voila!


