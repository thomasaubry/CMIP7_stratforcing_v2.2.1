% 	Written by Thomas J. Aubry, July 2025.
% 	Department of Earth Sciences, University of Oxford, UK
%   E-mail: thomas.aubry@earth.ox.ac.uk
% 	Please cite the corresponding paper if you use or modify this script


clear
close all


addpath('./source_datasets') %where emission source datasets are located on your computer


%a bit of a silly self-reminder to check source dataset versions are up to
%date but it could be useful to someone :)
x = input(...
    ['Are you happy with versions of core datasets used in this script version?' ...
    '\nThese include MSVOLSO2L4, eVolv2k, D4i and GVP datasets, as well as homemade spreadsheets'...
    '\nwith geochemical eruption match, historical eruption matches, and eruption source parameters.'...
    '\nType ''y'' for yes and ''n'' for no, then press the return key :)'],"s")
if strcmp(x,'n') | ~strcmp(x,'y')
    error('Then better fix these datasets ;)')
end


%%=========================================================================
%% 1) Import MSVOLSO2L4 (satellite emission inventory)
%%=========================================================================

%First load the .txt dataset
opts = delimitedTextImportOptions("NumVariables", 13);

% Specify range and delimiter
opts.DataLines = [2, Inf];
opts.Delimiter = ["\t", " "];

% Specify column names and types
opts.VariableNames = ["volcano", "lat", "lon", "v_alt", "yyyy", "mm", "dd", "type", "vei", "p_alt_obs", "p_alt_est", "so2kt", "VarName13"];
opts.VariableTypes = ["categorical", "double", "double", "double", "double", "double", "double", "categorical", "double", "double", "double", "double", "string"];

% Specify file level properties
opts.ExtraColumnsRule = "ignore";
opts.EmptyLineRule = "read";
opts.ConsecutiveDelimitersRule = "join";

% Specify variable properties
opts = setvaropts(opts, "VarName13", "WhitespaceRule", "preserve");
opts = setvaropts(opts, ["volcano", "type", "VarName13"], "EmptyFieldRule", "auto");

% Import the data
MSVOLSO2L4 = readtable("./source_datasets\MSVOLSO2L4_v04-00-2024m0129_noheader.txt", opts);
% Clear temporary variables
clear opts

%Then extract all variables from the
MSVOLSO2_volcano=MSVOLSO2L4.volcano;
MSVOLSO2_lat=MSVOLSO2L4.lat;%degree north
MSVOLSO2_lon=MSVOLSO2L4.lon;%degree east -180 to 180
MSVOLSO2_ventalt=MSVOLSO2L4.v_alt;% vent altitude km a.s.l.
MSVOLSO2_year=MSVOLSO2L4.yyyy;
MSVOLSO2_month=MSVOLSO2L4.mm;
MSVOLSO2_day=MSVOLSO2L4.dd;
MSVOLSO2_erutype=MSVOLSO2L4.type;%exp = explosive, eff = effusive.
MSVOLSO2_VEI=MSVOLSO2L4.vei;
MSVOLSO2_p_alt_obs=MSVOLSO2L4.p_alt_obs;%km a.s.l., observed from various sources
MSVOLSO2_p_alt_obs(MSVOLSO2_p_alt_obs==-999)=NaN;
MSVOLSO2_p_alt_est=MSVOLSO2L4.p_alt_est;%km a.s.l., estimated plume altitude (km). Fixed at 10 km above vent altitude for explosive eruptions and 5 km above vent altitude for effusive eruptions;
MSVOLSO2_p_alt_est(MSVOLSO2_p_alt_est==-999)=NaN;
MSVOLSO2_mass=MSVOLSO2L4.so2kt/1000;%kt SO2 converted to Tg SO2
MSVOLSO2_height=MSVOLSO2_p_alt_obs;
MSVOLSO2_height(isnan(MSVOLSO2_height))=MSVOLSO2_p_alt_est(isnan(MSVOLSO2_height));
MSVOLSO2_isuts=(abs(MSVOLSO2_lat)>20 & MSVOLSO2_height>=5) | (abs(MSVOLSO2_lat)<=20 & MSVOLSO2_height>=10);


%%=========================================================================
%% 2) Load eVolv2k (emission inventory based on bipolar ice core array)
%%=========================================================================
fname='eVolv2k_v3_ds_1.nc';

%use the bespoke varfromname function to extract variables from netcdf
eVolv2k_year=varfromname(fname,'year',0);%year ISO (year CE same for historical period)
eVolv2k_month=varfromname(fname,'month',0);
eVolv2k_day=varfromname(fname,'day',0);
eVolv2k_lat=varfromname(fname,'lat',0);
eVolv2k_mass=varfromname(fname,'vssi',0)*64.066/32.065;%in TgS, converted to TgSO2
eVolv2k_dmass=varfromname(fname,'sigma_vssi',0)*64.066/32.065;%in TgS, converted to TgSO2, 1-sigma uncertainty
eVolv2k_hemiratio=varfromname(fname,'hemi',0)*64.066/32.065;%hemispheric asymmetry (NH/SH) of aerosol spread for tropical eruptions based on ratio of Greenland to Antarctic flux

%define latitude band as 1 for SH, 2 for tropics and 3 for NH
eVolv2k_latband=NaN(size(eVolv2k_year));
eVolv2k_latband(eVolv2k_lat<-25)=1;eVolv2k_latband(abs(eVolv2k_lat)<=25)=2;eVolv2k_latband(eVolv2k_lat>25)=3;

%convert the eVolv2k deposition asymmetry factor (NH/SH) to the one we use
%((NH-SH)/(NH+SH))
eVolv2k_asym=varfromname(fname,'hemi',0);
eVolv2k_asym(eVolv2k_asym==-1 & eVolv2k_lat>0)=inf;%assume that -1 values with positive lat are those with SH=0;
eVolv2k_asym(eVolv2k_asym==-1 & eVolv2k_lat<=0)=0;%assume that -1 values with negative lat are those with SH=0;
eVolv2k_asym=(1./(1+1./eVolv2k_asym))-1./(eVolv2k_asym+1);

%%=========================================================================
%% 3) Load Sigl et al. (2015) (emission inventory based on bipolar ice core array)
%%=========================================================================
sigl2015=xlsread('Sigl2015.xlsx',2);


S2015_year=sigl2015(:,1);%year CE
S2015_latband=sigl2015(:,2);%latitudinal band; 1=tropical, 2=NH, 3=SH
%change the latitude band to convention we use here
S2015_latband=S2015_latband+1;S2015_latband(S2015_latband==4)=1;%1=SH, 2=tropical, 3=NH

S2015_mass=sigl2015(:,7+1)*64.066/96.06;%in Tg SO4, converted to Tg SO2
S2015_dmass=S2015_mass.*sigl2015(:,11+1)./sigl2015(:,9+1);%1-sigma uncertainty, assuming that the relative uncertainty on their forcing and total loading is the same which sounds valid as no forcing scaling factor uncertainty propagated in their methods

S2015_NHflux=sigl2015(:,3+1);S2015_SHflux=sigl2015(:,5+1);
S2015_NHflux(isnan(S2015_NHflux))=0;S2015_SHflux(isnan(S2015_SHflux))=0;
S2015_asym=(S2015_NHflux-S2015_SHflux)./(S2015_NHflux+S2015_SHflux);%calculate asymmetry factor


%Change in years for Sigl et al below were done for plotting purposes to
%compare the dataset. They do not affect the making of the dataset as
%changes in years for eruption sourced from ice core are ultimately all
% determined from our match spreadsheet. I leave the code below just in
% case deleting introduces any error, e.g. wrong ice-core year for Agung in
% match spreadsheet.

%For events for which for a year N, there is an eruption in evolv2k but not
%in Sigl et al. 2015, but that Sigl et al 2015 has an eruption at year N+1,
%we assume that the eruption year in Sigl et al 2015 is N instead of N+1.

for i=1:length(eVolv2k_year)
    if isempty(find(S2015_year==eVolv2k_year(i))) & ~isempty(find(S2015_year==eVolv2k_year(i)+1))
        S2015_year(find(S2015_year==eVolv2k_year(i)+1))=eVolv2k_year(i);
    end
end

%Correct the year of Agung from 1964 to 1963
S2015_year(S2015_year==1964)=1963;


%%=========================================================================
%% 4) Load D4i (emission inventory based on single Greenland core, but high-res)
%%=========================================================================

D4i=xlsread('D4i_cleaned_noevolv2k.xlsx');% this spreadsheet was pre-processed, e.g. by removing evolv2k events

maskd4i=D4i(:,1)>=1750;%we only retain eruptions for 1750 or later
D4i_year=D4i(maskd4i,1);

%Latitude band is a priori unknown for D4i, so we defined two potential
%masses based on the different conversion factors and the way mass was
%calculated in the D4i paper with NH vs tropical probability (see Fang et
%al 2023).
%if eruption is NH mass will be
D4i_mass_ifNH=D4i(maskd4i,2)*0.103*64.066/32.065;%in TgS, converted to TgSO2, if eruption is NH
%if eruption is SH, mass will be
D4i_mass_ifT=D4i(maskd4i,2)*0.271*64.066/32.065;%in TgS, converted to TgSO2, if eruption is T
%latitude band cannot be inferred for single Greenland ice-core, we enter 2.5 as 2 is tropical and 3 is NH
D4i_latband=2.5*ones(size(D4i_year));



%%=========================================================================
%% 5)Create base CMIP7 emission variables and combine eVolv2k (1750-1900) and Sigl et al. 2015 (1900-1978) and D4i (1730-1900)
%%=========================================================================

%we only use eVolv2k for 1750-1900
eVolv2k_mask=eVolv2k_year<=1900 & eVolv2k_year>=1750;
%we only use Sigl et al for 1901-1978
S2015_mask=S2015_year>1900 & S2015_year<1979;
%note D4i already filtered above to keep only 1750-1900 events, and
% events thoughts to be in eVolv2k were removed in the spreadsheet

%Combine all three ice-core datasets. For now assume D4i mass tropical but
%adequate mass will be used later
cmip7_mass=[D4i_mass_ifT;S2015_mass(S2015_mask);eVolv2k_mass(eVolv2k_mask)];
cmip7_dmass=[D4i_mass_ifT;S2015_dmass(S2015_mask);eVolv2k_dmass(eVolv2k_mask)];
cmip7_latband=[D4i_latband;S2015_latband(S2015_mask);eVolv2k_latband(eVolv2k_mask)];
cmip7_ice_year=[D4i_year;S2015_year(S2015_mask);eVolv2k_year(eVolv2k_mask)];
cmip7_asym=[NaN(size(D4i_year));S2015_asym(S2015_mask);eVolv2k_asym(eVolv2k_mask)];

%source parameter is a flag saying where key source data came from. 1 for
%eVolv2k, 2 for Sigl et al 2015, 3 for D4i; later down we will have 4 for
%GVP and 5 for MSVOLSO2L4
cmip7_source=[3*ones(size(D4i_year));2*ones(size(S2015_asym(S2015_mask)));ones(size(eVolv2k_asym(eVolv2k_mask)))];

%order the dataset chronologically
[a sortind]=sort(cmip7_ice_year);
cmip7_mass=cmip7_mass(sortind);
cmip7_dmass=cmip7_dmass(sortind);
cmip7_latband=cmip7_latband(sortind);
cmip7_ice_year=cmip7_ice_year(sortind);
cmip7_asym=cmip7_asym(sortind);
cmip7_source=cmip7_source(sortind);

%below save all ice-core emissions in a spreadsheet. This spreadsheet was
%used "offline" to match these ice-core events to known eruptive events

tomatchdata=cell(length(cmip7_mass)+1,5);
tomatchdata(1,:)={'Ice Year','Latitude band (1=SH,2=Trop,3=NH,2.5=Trop or NH)','Mass of SO2 (Tg SO2)','Asymetry factor (-1=SH,1=NH)','Source (1=evolv2k, 2=Sigl2015, 3=D4i'};
tomatchdata(2:end,:)=num2cell([cmip7_ice_year cmip7_latband cmip7_mass cmip7_asym cmip7_source]);
writecell(tomatchdata,'./source_datasets\CMIP7_presat_list.csv')


%%=========================================================================
%% 6)Load Global Volcanism Program (GVP) Volcanoes of the World (VOTW) dataset for confirmed eruptions
%%=========================================================================

%First load the xlsx file
% Set up the Import Options and import the data
opts = spreadsheetImportOptions("NumVariables", 24);

% Specify sheet and range
opts.Sheet = "Eruption List";
opts.DataRange = "A3";

% Specify column names and types
opts.VariableNames = ["VolcanoNumber", "VolcanoName", "EruptionNumber", "Var4", "Var5", "VEI", "Var7", "Var8", "StartYear", "Var10", "StartMonth", "Var12", "StartDay", "Var14", "Var15", "Var16", "Var17", "Var18", "Var19", "Var20", "Var21", "Var22", "Latitude", "Longitude"];
opts.SelectedVariableNames = ["VolcanoNumber", "VolcanoName", "EruptionNumber", "VEI", "StartYear", "StartMonth", "StartDay", "Latitude", "Longitude"];
opts.VariableTypes = ["double", "string", "double", "char", "char", "double", "char", "char", "double", "char", "double", "char", "double", "char", "char", "char", "char", "char", "char", "char", "char", "char", "double", "double"];

% Specify variable properties
opts = setvaropts(opts, ["VolcanoName", "Var4", "Var5", "Var7", "Var8", "Var10", "Var12", "Var14", "Var15", "Var16", "Var17", "Var18", "Var19", "Var20", "Var21", "Var22"], "WhitespaceRule", "preserve");
opts = setvaropts(opts, ["VolcanoName", "Var4", "Var5", "Var7", "Var8", "Var10", "Var12", "Var14", "Var15", "Var16", "Var17", "Var18", "Var19", "Var20", "Var21", "Var22"], "EmptyFieldRule", "auto");

% Import the data
gvp = readtable("./source_datasets\GVP_ConfirmedEruption_Search_Result_202501121044.xlsx", opts, "UseExcel", false);
% Clear temporary variables
clear opts


%Then extract variables of use
gvp_volcname=gvp.VolcanoName;
gvp_volcN=gvp.VolcanoNumber;
gvp_eruN=gvp.EruptionNumber;
gvp_vei=gvp.VEI;
gvp_lat=gvp.Latitude;
gvp_lon=gvp.Longitude;
gvp_startyear=gvp.StartYear;
gvp_startmonth=gvp.StartMonth;
gvp_startday=gvp.StartDay;


%%=========================================================================
%% 7) Assess GVP under-recording
%%Note to self: use GVP+LAMEVE in the future?
%%=========================================================================

%Here we assess under-recording per VEI in the GVP dataset. This helps us
%understand how reliable this dataset is overall (of course this is well
%studied but doing our own plots), and in turn later down to think about
%how confident we are when matching an eruption to an ice-core sulfate
%deposition event


%we are interested in eruptions over last 10,000 years and VEI 3-6
mintimegvp=-10000;
veilist=[3:6];

%fitlim is the duration over which we think a given VEI has no
%underrecording, wrt year 2022. For the smallest VEI, we assume that there
%is no underrecording over satellite era, i.e. last ~42 years

fitlim=42%
for i=1:4


    figure(12)
    subplot(2,2,i)
    %find eruption of VEI of interest during period of interest
    mask=gvp_startyear>=mintimegvp & gvp_startyear<= 2022 & gvp_vei==veilist(i);
    %find list of eruption start year
    xx=gvp_startyear(mask);
    %plot cumulative number. In the absence of underrecording and assuming
    %no change in frequency-magnitude distribution over last 10000 years,
    %this should look linear
    her=plot(gvp_startyear(mask),cumsum(ones(size(gvp_startyear(mask))),'reverse'),'k.')
    fitlim0=fitlim;
    yy=cumsum(ones(size(gvp_startyear(mask))),'reverse');



    %for the duration for which we feel confident there is no
    %underrecording, do linear fit of cumulative sum of event vs time
    fitvei=fit(xx(xx>=2022-fitlim),yy(xx>=2022-fitlim),'poly1');
    %calculate return period as 1/slope of fit
    eruperiod=round(1./fitvei.p1,1);
    if veilist(i)<9

        %look at fit prediction/confidence interval at 95% over the full
        %time period and plot
        int_vei = predint(fitvei,xx,0.95,'observation','on');
        hold on
        hconf=patch([xx;flip(xx)], [int_vei(:,1);flip(int_vei(:,2))], 'r', 'FaceAlpha',0.3, 'EdgeColor','none')
        hfit=plot(xx,fitvei(xx),'r-','LineWidth',1.5)

        %find the year before which the cumulative number of eruption is
        %smaller than permitted within fit prediction interval. This is the
        %year from which we deem there is significant underrecording
        ind=find(yy<int_vei(:,2));

        %update the value of fitlim based on the year above. The assumption
        %is that if there is no underrecording for VEI N from a certain
        %year, there also won't be underrecording for VEI N+1 over the same
        %time period, so we can use this time period to do a linear fit of
        %cumulative freq vs time. For our first VEI, we used satellite era
        %as the period with no underrecording (see above)
        fitlim=(2022-xx(ind(end)));
        hlim=plot([2022-fitlim 2022-fitlim],[0 max(yy)*1.02],'k--')
        fitlim=2*fitlim;
    end

    %just adding to the plot, including return freq and year from which no
    %underrecording
    xlim([max(2022-fitlim,mintimegvp) 2023])
    ylim([0.98*min(yy(xx>=2022-fitlim))-1 max(yy)*1.02+1])
    ylabel({'Cumulative number of';'eruptions since -10,000'})
    title(strcat('VEI',{' '},num2str(veilist(i))))
    aa=strcat('Linear fit',{' '},num2str(2022-fitlim0),'-2022 (return period:',{' '},num2str(eruperiod),{' '},'yr)');
    if i==1
        legend([her hfit hconf hlim],'VOTW eruptions',aa{1},'Fit uncertainty','Time limit under-recording','Location','NorthWest')
    else
        legend([hfit],aa{1},'Location','NorthWest')
    end
end


%%=========================================================================
%% 8)Load the eruption match spreadsheet
%%=========================================================================
%This spreadsheet is prepared externally based on the CMIP7_presat_list.csv
% file saved previously,  filled with proposed eruption matches for sulfate
% deposition peak found in ice-core

% Set up the Import Options
opts = spreadsheetImportOptions("NumVariables", 10);

% Specify sheet and range
opts.Sheet = "CMIP7_presat_list";
opts.DataRange = "A2";
% Specify column names and types
opts.VariableNames = ["IceYear", "LatitudeBand1SH2Trop3NH25TropOrNH", "MassOfSO2TgSO2", "AsymetryFactor1SH1NH", "Source1evolv2k2Sigl20153D4i", "MatchEruptionNumber", "MatchVolcanoNumber", "VolcanoNameAndEruptionYear", "MatchConfidence", "Comment"];
opts.VariableTypes = ["double", "double", "double", "double", "double", "double", "double", "string", "double", "string"];

% Specify variable properties
opts = setvaropts(opts, ["VolcanoNameAndEruptionYear", "Comment"], "WhitespaceRule", "preserve");
opts = setvaropts(opts, ["VolcanoNameAndEruptionYear", "Comment"], "EmptyFieldRule", "auto");
% Import the data
matchlist = readtable("./source_datasets\CMIP7_presat_list_match.xlsx", opts ,"UseExcel", false);
% Clear temporary variables
clear opts

%extract variables of interest
matchlist_iceyear=matchlist.IceYear;
matchlist_mass=matchlist.MassOfSO2TgSO2;
matchlist_volcN=matchlist.MatchVolcanoNumber;
matchlist_eruN=matchlist.MatchEruptionNumber;
matchlist_matchinfo=matchlist.Comment;
matchlist_matchconf=matchlist.MatchConfidence;

%%=========================================================================
%% 9)Load the eruption source parameter (ESP) spreadsheet
%%=========================================================================
%This spreadsheet is
%prepared externally and contains key ESPs needed for climate modelling,
%after careful literature review of all eruptions of interest i.e. all VEI
%4+ eruptions between 1750 and 1978, plus any other eruption (VEI 3) deemed
% as a potential match

% Set up the Import Options and import the data
opts = spreadsheetImportOptions("NumVariables", 15);

% Specify sheet and range
opts.Sheet = "Eruption List";
opts.DataRange = "A2";

% Specify column names and types
opts.VariableNames = ["VolcanoNumber", "VolcanoName", "EruptionNumber", "VEI", "StartYear", "EndYear", "Latitude", "Longitude", "EruptionYear", "EruptionMonth", "EruptionDay", "SO2InjectionHeightkmAsl", "SO2InjectionHeightUncertaintykm", "VentAltitudekmAslFromGVP", "CommentForEruptionSourceParametersdateAndHeight"];
opts.VariableTypes = ["double", "string", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "string"];

% Specify variable properties
opts = setvaropts(opts, ["VolcanoName", "CommentForEruptionSourceParametersdateAndHeight"], "WhitespaceRule", "preserve");
opts = setvaropts(opts, ["VolcanoName", "CommentForEruptionSourceParametersdateAndHeight"], "EmptyFieldRule", "auto");

% Import the data
esplist = readtable("./source_datasets\CMIP7_ESP.xlsx", opts, "UseExcel", false);
% Clear temporary variables
clear opts


%extract parameters of interest;
esplist_year=esplist.EruptionYear;
esplist_month=esplist.EruptionMonth;
esplist_day=esplist.EruptionDay;
esplist_volcN=esplist.VolcanoNumber;
esplist_eruN=esplist.EruptionNumber;
esplist_height=esplist.SO2InjectionHeightkmAsl;
esplist_dheight=esplist.SO2InjectionHeightUncertaintykm;
esplist_valt=esplist.VentAltitudekmAslFromGVP;
esplist_espinfo=esplist.CommentForEruptionSourceParametersdateAndHeight;
esplist_VEI=esplist.VEI;
esplist_lat=esplist.Latitude;
esplist_lon=esplist.Longitude;
esplist_volcname=esplist.VolcanoName;
esplist_startyear=esplist.StartYear;

%despite our best effort to fill eruption source parameters, some eruptions
%in GVP might not even have a well-constrained start year, at least for
%explosive/paroxysmal phase. Here we check the user is happy with that and
%to assume than when we don't have a good start year, it's ok to assume the
%reported GVP start year is adequate. We did so to make the CMIP7 dataset,
%but anyone developping/expanding our dataset might want to be aware/push
%ESP search so we leave this check uncommented.

if sum(isnan(esplist_year))>0
    x = input('Warning! Some ESP searches are incomplete. Do you want to continue with the eruption year equal to GVP start year for incomplete ESP searches? (y/n)',"s")
    if strcmp(x,'y')
        esplist_year(isnan(esplist_year))=esplist_startyear(isnan(esplist_year));
    else
        error('Then go fix your broken ESP spreadsheet ;)')
    end
end


%%=========================================================================
%% 10)Pre-allocate CMIP7 parameters and filter GVP data
%%=========================================================================
% Filter GVP data to only keep VEI4+ and 1750-1978 events. These are
% eruptions that will be added to the pre-satellite era emission inventory
% should they not be one of the matches. Also remove  volcano number 600000
% which is unknown volcano, including mystery 1809 eruption

maskgvp = gvp_vei>=4 & gvp_startyear<=1978 & gvp_startyear>=1750 & gvp_volcN~=600000;
gvp4emis_volcname=gvp_volcname(maskgvp);
gvp4emis_volcN=gvp_volcN(maskgvp);
gvp4emis_eruN=gvp_eruN(maskgvp);
gvp4emis_vei=gvp_vei(maskgvp);
gvp4emis_lat=gvp_lat(maskgvp);
gvp4emis_lon=gvp_lon(maskgvp);
gvp4emis_startyear=gvp_startyear(maskgvp);
gvp4emis_startmonth=gvp_startmonth(maskgvp);
gvp4emis_startday=gvp_startday(maskgvp);
gvp4emis = table(gvp4emis_volcname,gvp4emis_volcN,gvp4emis_eruN,gvp4emis_vei,gvp4emis_lat,gvp4emis_lon,gvp4emis_startyear,gvp4emis_startmonth,gvp4emis_startday);


% Then we pre-allocate all parameters for pre-satellite era eruption

cmip7_year=NaN(size(cmip7_mass));
cmip7_month=NaN(size(cmip7_mass));
cmip7_day=NaN(size(cmip7_mass));
cmip7_height=NaN(size(cmip7_mass));
cmip7_dheight=NaN(size(cmip7_mass));
cmip7_ventalt=NaN(size(cmip7_mass));
cmip7_volcN=NaN(size(cmip7_mass));
cmip7_eruN=NaN(size(cmip7_mass));
cmip7_VEI=NaN(size(cmip7_mass));
cmip7_lat=NaN(size(cmip7_mass));
cmip7_lon=NaN(size(cmip7_mass));
cmip7_volcname=cell(size(cmip7_mass));
cmip7_espinfo=cell(size(cmip7_mass));

%for eruptions derived from ice-core matched to an eruption, a cell will
%contain all comments explaining the match, and a confidence flag will
%indicate how confident we are that the match is correct. Flag 1=plausible
% match, 2=dedium confidence match, 3=high confidence match, 4=very high confidence match
cmip7_matchconf=NaN(size(cmip7_mass));
cmip7_matchinfo=cell(size(cmip7_mass));

%%=========================================================================
%% 11) Match pre-satellite ice-core data to GVP database
%%=========================================================================

%for each event in the match list table, do the following...
for imatch=1:length(matchlist_iceyear)

    %find the cmip7 event with matching ice information
    icmip7=find(abs(matchlist_iceyear(imatch)-cmip7_ice_year)<0.01 & abs(matchlist_mass(imatch)-cmip7_mass)<0.001);
    if length(icmip7)~=1
        error(strcat('Error: No CMIP7 event or multiple ones have ice year',num2str(matchlist_iceyear(imatch)),' and SO2 mass',num2str(matchlist_mass(imatch)),' Tg SO2! Please check the match spreadsheet.'))
    end
    cmip7_matchinfo{icmip7}=matchlist_matchinfo(imatch);

    %Then if there is a match i.e. the volcano number is not empty,
    if ~isnan(matchlist_volcN(imatch))

        %find the gvp event in the match ESP spreadsheet with matching eruption number and volcano number
        iesp=find(matchlist_eruN(imatch)==esplist_eruN & matchlist_volcN(imatch)==esplist_volcN);

        %the 1831 has no eruption number yet so it's possible to have volcano
        %number but no eruption number
        if isnan(matchlist_eruN(imatch))
            iesp=find(isnan(esplist_eruN) & matchlist_volcN(imatch)==esplist_volcN);
        end

        if isempty(iesp)
            error(strcat('Error: No information in the ESP spreadsheet for the requested volcano number',num2str(matchlist_volcN(imatch)),' and eruption number',num2str(matchlist_eruN(imatch)),'! Please check the ESP spreadsheet.'))
        end

        %attribute parameters of eruption match including match confidence from
        %spreadsheet

        %if there is a single set of eruption source parameter (i.e. a single event)
        % for the eruption, things are easy. The CMIP7 parameters for the event are
        % simply the match parameters
        if length(iesp)==1
            cmip7_matchconf(icmip7)=matchlist_matchconf(imatch);

            cmip7_year(icmip7)=esplist_year(iesp);
            cmip7_month(icmip7)=esplist_month(iesp);
            cmip7_day(icmip7)=esplist_day(iesp);
            cmip7_height(icmip7)=esplist_height(iesp);
            cmip7_dheight(icmip7)=esplist_dheight(iesp);
            cmip7_ventalt(icmip7)=esplist_valt(iesp);
            cmip7_volcN(icmip7)=esplist_volcN(iesp);
            cmip7_eruN(icmip7)=esplist_eruN(iesp);
            cmip7_VEI(icmip7)=esplist_VEI(iesp);
            cmip7_lat(icmip7)=esplist_lat(iesp);
            cmip7_lon(icmip7)=esplist_lon(iesp);
            cmip7_volcname{icmip7}=esplist_volcname(iesp);
            cmip7_espinfo{icmip7}=esplist_espinfo(iesp);
        end

        %if there is are multiple sets of eruption source parameter (i.e. multiple
        %events) for the eruption, things are a bit more complicated. This happens
        %if an eruption has multiple phases. At v2.2.1, pre-satellite, this is only
        %the case for Agung and Laki so first quick sanity check that we are
        %looking at one of these

        if length(iesp)>1
            ll=length(iesp);

            checklaki= ll==10 & matchlist_volcN==373010 & matchlist_eruN==12809;
            checkagung= ll==2 & matchlist_volcN==264020 & matchlist_eruN==16210;
            if ~checklaki & ~checkagung
                error(strcat('The event with ice year',num2str(matchlist_iceyear(imatch)),' and SO2 mass',num2str(matchlist_mass(imatch)),' Tg SO2 has more than 1 (',num2str(ll),') eruptive events with associated ESPs. It should only be the case for Agung 1963 (2) and Laki 1783 (10).'))
            end

            %if that's indeed the case, for these events, we must insert new rows in
            %the parameter variables as there was one row per eruption, but in the end
            %we want one row per (eruptive) event. So if N is the number of events for
            %the considered eruption, for each parameter, we insert N-1 row. Then for
            %each parameters, we equate the N values for the eruption considered in the
            %CMIP7 spreadsheet to the corresponding N values in the ESP spreadsheet

            cmip7_matchconf=[cmip7_matchconf(1:icmip7);NaN(ll-2,1);cmip7_matchconf(icmip7:end)] ;
            cmip7_matchconf(icmip7:icmip7+ll-1)=matchlist_matchconf(imatch);

            cmip7_year=[cmip7_year(1:icmip7);NaN(ll-2,1);cmip7_year(icmip7:end)];
            cmip7_year(icmip7:icmip7+ll-1)=esplist_year(iesp);

            cmip7_month=[cmip7_month(1:icmip7);NaN(ll-2,1);cmip7_month(icmip7:end)];
            cmip7_month(icmip7:icmip7+ll-1)=esplist_month(iesp);

            cmip7_day=[cmip7_day(1:icmip7);NaN(ll-2,1);cmip7_day(icmip7:end)];
            cmip7_day(icmip7:icmip7+ll-1)=esplist_day(iesp);

            cmip7_height=[cmip7_height(1:icmip7);NaN(ll-2,1);cmip7_height(icmip7:end)];
            cmip7_height(icmip7:icmip7+ll-1)=esplist_height(iesp);

            cmip7_dheight=[cmip7_dheight(1:icmip7);NaN(ll-2,1);cmip7_dheight(icmip7:end)];
            cmip7_dheight(icmip7:icmip7+ll-1)=esplist_dheight(iesp);

            cmip7_ventalt=[cmip7_ventalt(1:icmip7);NaN(ll-2,1);cmip7_ventalt(icmip7:end)];
            cmip7_ventalt(icmip7:icmip7+ll-1)=esplist_valt(iesp);

            cmip7_volcN=[cmip7_volcN(1:icmip7);NaN(ll-2,1);cmip7_volcN(icmip7:end)];
            cmip7_volcN(icmip7:icmip7+ll-1)=esplist_volcN(iesp);

            cmip7_eruN=[cmip7_eruN(1:icmip7);NaN(ll-2,1);cmip7_eruN(icmip7:end)];
            cmip7_eruN(icmip7:icmip7+ll-1)=esplist_eruN(iesp);

            cmip7_VEI=[cmip7_VEI(1:icmip7);NaN(ll-2,1);cmip7_VEI(icmip7:end)];
            cmip7_VEI(icmip7:icmip7+ll-1)=esplist_VEI(iesp);

            cmip7_lat=[cmip7_lat(1:icmip7);NaN(ll-2,1);cmip7_lat(icmip7:end)];
            cmip7_lat(icmip7:icmip7+ll-1)=esplist_lat(iesp);

            cmip7_lon=[cmip7_lon(1:icmip7);NaN(ll-2,1);cmip7_lon(icmip7:end)];
            cmip7_lon(icmip7:icmip7+ll-1)=esplist_lon(iesp);

            cmip7_volcname=[cmip7_volcname(1:icmip7);cell(ll-2,1);cmip7_volcname(icmip7:end)];
            cmip7_volcname(icmip7:icmip7+ll-1)=cellstr(esplist_volcname(iesp));

            cmip7_espinfo=[cmip7_espinfo(1:icmip7);cell(ll-2,1);cmip7_espinfo(icmip7:end)];
            cmip7_espinfo(icmip7:icmip7+ll-1)=cellstr(esplist_espinfo(iesp));

            cmip7_matchinfo=[cmip7_matchinfo(1:icmip7);cell(ll-2,1);cmip7_matchinfo(icmip7:end)];
            cmip7_matchinfo(icmip7:icmip7+ll-1)=cellstr(cmip7_matchinfo(icmip7));

            %When there are multiple eruptive events from an eruption matched to an ice
            % core sulfate deposition event they all have the  same ice_year, source,
            %asymetry factor and latitude bands
            cmip7_latband=[cmip7_latband(1:icmip7);cmip7_latband(icmip7)*ones(ll-2,1);cmip7_latband(icmip7:end)];
            cmip7_ice_year=[cmip7_ice_year(1:icmip7);cmip7_ice_year(icmip7)*ones(ll-2,1);cmip7_ice_year(icmip7:end)];
            cmip7_asym=[cmip7_asym(1:icmip7);cmip7_asym(icmip7)*ones(ll-2,1);cmip7_asym(icmip7:end)];
            cmip7_source=[cmip7_source(1:icmip7);cmip7_source(icmip7)*ones(ll-2,1);cmip7_source(icmip7:end)];

            %the SO2 mass (and uncertainty) from ice core will need to be partitioned,
            %but this is done just below and for now we assume each event has the total
            %mass recorded in ice-core and corresponding uncertainty
            cmip7_mass=[cmip7_mass(1:icmip7);cmip7_mass(icmip7)*ones(ll-2,1);cmip7_mass(icmip7:end)];
            cmip7_dmass=[cmip7_dmass(1:icmip7);cmip7_dmass(icmip7)*ones(ll-2,1);cmip7_dmass(icmip7:end)];

        end


        %Last, if a GVP event was matched to an ice core event, we remove it from the list
        %of eruptions to be added to the emission database to avoid double-counting
        igvp=find(matchlist_eruN(imatch)==gvp4emis.gvp4emis_eruN & matchlist_volcN(imatch)==gvp4emis.gvp4emis_volcN);
        if ~isempty(igvp)
            gvp4emis(igvp,:)=[];
        end

    end
end


%Last (for real this time), for Laki and Agung, we need to distribute the
%SO2 mass (and uncertainty) among multiple events. This is done below. The
%distribution for Laki is based on lava volume for each phase, whereas for
%Agung it's based on mass eruption rate. Further comments are given in the
%ESP info cell

%Laki
ilaki=find(cmip7_volcN==373010 & cmip7_eruN==12809);
laki_frac=[2.9 4.4 5.9 7.7 13.2 8.9 10.8 18.7 13.5 8.3]';
laki_frac=laki_frac/sum(laki_frac);
dmass_rel_laki=cmip7_dmass(ilaki(1))/cmip7_mass(ilaki(1));
cmip7_mass(ilaki)=cmip7_mass(ilaki).*laki_frac;
cmip7_dmass(ilaki)=cmip7_mass(ilaki)*sqrt(dmass_rel_laki^2+0.3);

%Agung
iagung=find(cmip7_volcN==264020 & cmip7_eruN==16210);
agung_frac=[3.7 3.5]';
agung_frac=agung_frac/sum(agung_frac);
dmass_rel_agung=cmip7_dmass(iagung(1))/cmip7_mass(iagung(1));
cmip7_mass(iagung)=cmip7_mass(iagung).*agung_frac;
cmip7_dmass(iagung)=cmip7_mass(iagung)*sqrt(dmass_rel_agung^2+0.3);



%%=========================================================================
%% 12) Attribute ESPs to remaining GVP events and merge with CMIP7 dataset
%%=========================================================================
%To summarize where we are now: we've merged ice core datasets together, we
%identified eruption matches, and we gave them eruption source parameters
%corresponding to these matches. What remains pre-satellite era are all GVP
%events of VEI4+ that were not matched to an ice core events, i.e. the
%gvp4emis table as we removed any event matched to ice core previously. So
%for all these events, we need to extract parameters from GVP, add
%parameters coming from our ESP search, then merge with the rest of
%pre-satellite events. Let's do it! :)


%First sanity check that all VEI 6 are matched to an ice-core sulfate
%deposition event.
if max(gvp4emis.gvp4emis_vei)>5
    error('At least one VEI6+ eruption was not matched to an ice-core event. That''s suspicious...')
end


%Load all parameters from the ESP spreadsheet
gvp4emis_volcname=gvp4emis.gvp4emis_volcname;
gvp4emis_volcN=gvp4emis.gvp4emis_volcN;
gvp4emis_eruN=gvp4emis.gvp4emis_eruN;
gvp4emis_vei=gvp4emis.gvp4emis_vei;
gvp4emis_lat=gvp4emis.gvp4emis_lat;
gvp4emis_lon=gvp4emis.gvp4emis_lon;
gvp4emis_startyear=gvp4emis.gvp4emis_startyear;
gvp4emis_startmonth=gvp4emis.gvp4emis_startmonth;
gvp4emis_startday=gvp4emis.gvp4emis_startday;
gvp4emis_startmonth(gvp4emis_startmonth==0)=NaN;
gvp4emis_startday(gvp4emis_startday==0)=NaN;

%Additional parameters will come from our ESP search
gvp4emis_year=NaN(size(gvp4emis.gvp4emis_vei));
gvp4emis_month=NaN(size(gvp4emis.gvp4emis_vei));
gvp4emis_day=NaN(size(gvp4emis.gvp4emis_vei));
gvp4emis_height=NaN(size(gvp4emis.gvp4emis_vei));
gvp4emis_dheight=NaN(size(gvp4emis.gvp4emis_vei));
gvp4emis_valt=NaN(size(gvp4emis.gvp4emis_vei));
gvp4emis_espinfo=cell(size(gvp4emis.gvp4emis_vei));

%Additional parameters, specific to ice-core events, are not defined for GVP
gvp4emis_asym=NaN(size(gvp4emis.gvp4emis_vei));
gvp4emis_ice_year=NaN(size(gvp4emis.gvp4emis_vei));
gvp4emis_matchconf=NaN(size(gvp4emis.gvp4emis_vei));
gvp4emis_matchinfo=cell(size(gvp4emis.gvp4emis_vei));
gvp4emis_matchinfo(:)={'NA - Eruption from GVP'};

%The latitude band is defined but simply assessed based on known latitude
%instead of ice-core deposition asymetry
gvp4emis_latband=NaN(size(gvp4emis.gvp4emis_vei));
gvp4emis_latband(gvp4emis_lat<-25)=1;gvp4emis_latband(abs(gvp4emis_lat)<=25)=2;gvp4emis_latband(gvp4emis_lat>25)=3;

%Source flag: 4 means GVP
gvp4emis_source=4*ones(size(gvp4emis.gvp4emis_vei)); %4=from GVP

%We don't know the SO2 mass for GVP derived event. For now arbitrarly
%assume 0.1 Tg SO2 if VEI 4 and 1 Tg if VEI 5, but these are placeholder
%masses. Something slightly more sensible will be done in the script to
%build the SAOD dataset and, at the same time, the final emission dataset
gvp4emis_mass=NaN(size(gvp4emis.gvp4emis_vei));
VEI4mass=0.1;%assumed VEI 4 mass if no constraint from ice-core or satellite
gvp4emis_mass(gvp4emis_vei==4)=VEI4mass;
gvp4emis_mass(gvp4emis_vei==5)=VEI4mass*10;%assume VEI 5 mass is 10 times VEI 4 on average
gvp4emis_dmass=gvp4emis_mass;%assume 100% uncertainty but note the upper bound is a lot more than twice the assumed mass.



for igvp=1:length(gvp4emis_vei)

    %for each event in the GVP list, we look whether we have a match in our
    %ESP spreadsheet
    iesp=find(gvp4emis_eruN(igvp)==esplist_eruN & gvp4emis_volcN(igvp)==esplist_volcN);
    if ~isempty(iesp)
        if length(iesp)==1

            %if we do have a match, this gives us a specific year, month, day, height,
            %height uncertainty, vent altitude (if these parameters were found) and
            % ESP info so we fill this in
            gvp4emis_year(igvp)=esplist_year(iesp);
            if isnan(gvp4emis_year(igvp));gvp4emis_year(igvp)=gvp4emis_startyear(igvp);end
            gvp4emis_month(igvp)=esplist_month(iesp);
            if isnan(gvp4emis_month(igvp));gvp4emis_month(igvp)=gvp4emis_startmonth(igvp);end
            gvp4emis_day(igvp)=esplist_day(iesp);
            if isnan(gvp4emis_day(igvp));gvp4emis_day(igvp)=gvp4emis_startday(igvp);end
            gvp4emis_height(igvp)=esplist_height(iesp);
            gvp4emis_dheight(igvp)=esplist_dheight(iesp);
            gvp4emis_valt(igvp)=esplist_valt(iesp);
            gvp4emis_espinfo{igvp}=esplist_espinfo(iesp);

            %check that there was only one event in ESP spreadsheet for the GVP event.
            %In theory there could be multiple ones (it would be an eruption with
            %multiple phases and parameters available for each of those), but it's not
            %the case in our dataset at the moment so best to check there was no mix-up
        else
            error('Multiple events in the ESP spreadsheet for this GVP event. Go fix your ESP spreadsheet :p')
        end

    end
end

%We now have source parameters for all GVP events. So we can append them to
%the existing list of CMIP7 parameters, which so far contained all ice-core
%derived events
cmip7_matchconf=[cmip7_matchconf;gvp4emis_matchconf];
cmip7_year=[cmip7_year;gvp4emis_year];
cmip7_month =[cmip7_month;gvp4emis_month];
cmip7_day=[cmip7_day;gvp4emis_day];
cmip7_height=[ cmip7_height;gvp4emis_height];
cmip7_dheight=[cmip7_dheight;gvp4emis_dheight];
cmip7_ventalt=[cmip7_ventalt;gvp4emis_valt];
cmip7_volcN=[cmip7_volcN;gvp4emis_volcN];
cmip7_eruN=[cmip7_eruN;gvp4emis_eruN];
cmip7_VEI=[cmip7_VEI;gvp4emis_vei];
cmip7_lat=[cmip7_lat;gvp4emis_lat];
cmip7_lon=[cmip7_lon;gvp4emis_lon];
cmip7_volcname=[cmip7_volcname;cellstr(gvp4emis_volcname)];
cmip7_espinfo=[cmip7_espinfo;cellstr(gvp4emis_espinfo)];
cmip7_matchinfo=[cmip7_matchinfo;cellstr(gvp4emis_matchinfo)];
cmip7_mass=[cmip7_mass;gvp4emis_mass];
cmip7_dmass=[cmip7_dmass;gvp4emis_dmass];
cmip7_latband=[cmip7_latband;gvp4emis_latband];
cmip7_ice_year=[cmip7_ice_year;gvp4emis_ice_year];
cmip7_asym=[cmip7_asym;gvp4emis_asym];
cmip7_source=[cmip7_source;gvp4emis_source];


%%=========================================================================
%% 13) Combine pre-satellite and satellite dataset, and create all CMIP7 emission variables
%%=========================================================================

%Progress! We are done with the pre-satellite era part. So all we need to
%do now is add satellite era i.e. the MSVOLSO2L4 dataset :)

%We filter out events before 1979 or after 2023, events
%with a VEI smaller than 2, effusive events, and events with a plume height
%extremely unlikely to be stratospheric. For the latter constraint, we
%find that the minimum value of the monthly zonal mean lapse-rate tropopause
%is well bounded by the simple  function hmin=8+8*exp(-(lat/35)^6). We then
%filter out any event smaller than 50% of this height to account for
%additional temporal/zonal tropopause height variability, error in plume
%height and secondary transport processes e.g. lofting of volcanic cloud
%owing to radiation absorption by gas and/or aerosol.
MSVOLSO2_mask=MSVOLSO2_mass>0 & MSVOLSO2_year>=1979 & MSVOLSO2_year<=2023 & MSVOLSO2_VEI>=3 &  strcmp('exp',string(MSVOLSO2_erutype)) & MSVOLSO2_height>(5+5*exp(-(MSVOLSO2_lat/35).^6));

MSVOLSO2_volcano=cellstr(MSVOLSO2_volcano(MSVOLSO2_mask));
MSVOLSO2_lat=MSVOLSO2_lat(MSVOLSO2_mask);
MSVOLSO2_lon=MSVOLSO2_lon(MSVOLSO2_mask);
MSVOLSO2_year=MSVOLSO2_year(MSVOLSO2_mask);
MSVOLSO2_month=MSVOLSO2_month(MSVOLSO2_mask);
MSVOLSO2_day=MSVOLSO2_day(MSVOLSO2_mask);
MSVOLSO2_VEI=MSVOLSO2_VEI(MSVOLSO2_mask);
MSVOLSO2_mass=MSVOLSO2_mass(MSVOLSO2_mask);
MSVOLSO2_height=MSVOLSO2_p_alt_obs(MSVOLSO2_mask);%do not use height attributed by default from this dataset, only observed heights
MSVOLSO2_ventalt=MSVOLSO2_ventalt(MSVOLSO2_mask);
MSVOLSO2_latband=NaN(size(MSVOLSO2_lat));
MSVOLSO2_latband(MSVOLSO2_lat<-25)=1;MSVOLSO2_latband(abs(MSVOLSO2_lat)<=25)=2;MSVOLSO2_latband(MSVOLSO2_lat>25)=3;

%Then we append MSVOLSO2L4 parameters to CMIP7 ones...
cmip7_ventalt=[cmip7_ventalt;MSVOLSO2_ventalt];
cmip7_year=[cmip7_year;MSVOLSO2_year];
cmip7_month=[cmip7_month;MSVOLSO2_month];
cmip7_day=[cmip7_day;MSVOLSO2_day];
cmip7_mass=[cmip7_mass;MSVOLSO2_mass];
cmip7_dmass=[cmip7_dmass;NaN(size(MSVOLSO2_year))];
cmip7_height=[cmip7_height;MSVOLSO2_height];
cmip7_dheight=[cmip7_dheight;NaN(size(MSVOLSO2_year))];
cmip7_latband=[cmip7_latband;MSVOLSO2_latband];
cmip7_lat=[cmip7_lat;MSVOLSO2_lat];
cmip7_lon=[cmip7_lon;MSVOLSO2_lon];
cmip7_asym=[cmip7_asym;NaN(size(MSVOLSO2_year))];
cmip7_eruN=[cmip7_eruN;NaN(size(MSVOLSO2_year))];
cmip7_volcN=[cmip7_volcN;NaN(size(MSVOLSO2_year))];
cmip7_volcname=[cmip7_volcname;MSVOLSO2_volcano];
satmatch    = cell(size(MSVOLSO2_year)); satmatch(:) = {'NA, eruption from MSVOLSO2L4.'};
cmip7_matchinfo=[cmip7_matchinfo;satmatch];
cmip7_espinfo=[cmip7_espinfo;cell(size(MSVOLSO2_year))];
cmip7_VEI=[cmip7_VEI;MSVOLSO2_VEI];
cmip7_iceyear=[cmip7_ice_year;NaN(size(MSVOLSO2_year))];
cmip7_matchconf=[cmip7_matchconf;NaN(size(MSVOLSO2_year))];
cmip7_source=[cmip7_source;5*ones(size(MSVOLSO2_lat))];%source flag is 5 for MSVOLSO2L4


%%=========================================================================
%% 14)Combine MSVOLSO2 events belonging to the same eruption
%%=========================================================================

%We've got our dataset complete!! Except it of course had lots of gaps...we
%could not find all eruption source parameters pre-satellite era for all
%eruptions, ice-core event with no eruption match need parameters, and even
%some satellite events have no height. So time to fill these gaps! To do
%so, we will get a bit of extra data (pre-processing)

%first, we trust our satellite era data (at least more than the
%pre-satellite era one). So we will use it to derive some empirical
%relationship for gap filling. However, MSVOLSO2L4 typically has daily
%measurmeent for long eruptions, so some eruptions have 30+ events. Before
%analyzing the dataset, we will thus combine eruption phases together


%define some table that will be filled with combined parameter values from MSVOLSO2L4
comb_MSVOLSO2_mass=[];
comb_MSVOLSO2_height=[];
comb_MSVOLSO2_year=[];
comb_MSVOLSO2_VEI=[];
comb_MSVOLSO2_lat=[];
comb_MSVOLSO2_matchconf=[];

%find list of eruption from MSVOLSO2L4. As it does not report GVP volcano
%and eruption number, we assume that two events with the same latitude,
%longitude and VEI could belong to the same eruption.
erulist=unique([MSVOLSO2_lat MSVOLSO2_lon MSVOLSO2_VEI],'rows');
%reformat date
MSVOLSO2dateday=datenum(MSVOLSO2_year,MSVOLSO2_month,MSVOLSO2_day);
%define index list
indexlist=(1:length(MSVOLSO2_lat))';


%then, we will go through the full list of event.
%the while loop above will keep going until all indices in indexlist have
%been voided, and we void indices everytime they are combines into one
%event
while length(indexlist)>0
    i=indexlist(1);

    %for the first event in our list, we find the indices of all events potentially from
    %the same eruption i.e. same lat/lon/VEI
    idx2=find((1:length(MSVOLSO2_lat))'>=i & MSVOLSO2_lat==MSVOLSO2_lat(i) & MSVOLSO2_lon==MSVOLSO2_lon(i) & MSVOLSO2_VEI==MSVOLSO2_VEI(i) & ~isnan(MSVOLSO2_lat(i)) & ~isnan(MSVOLSO2_lat(i)) & ~isnan(MSVOLSO2_lon(i)));
    %then, among those, we check if the difference between two subsequent
    %event is more than 5 days. If it is not we will arbitrarily assume
    %it's the same eruption.
    differentevent=find(diff(MSVOLSO2dateday(idx2))>5);%boolean; true = more than 5 day difference

    %if no event is apart by more than 5 days, then we will combine all in
    %one eruption
    if isempty(differentevent)
        imax=length(idx2);
        %but if some events are apart by more than 5 days, let's only consider the indices up to
        %the first index for which duration difference is more than 5 days, i.e.
        %we just consider the first eruption
    else;  imax=differentevent(1);end;

    %now we redefine the list of index to only consider the first eruption in
    %idx2, which was the list of events originating from same lat/lon/VEI and
    %less than 5 days apart
    idx2=idx2(1:imax);

    %we're about to deal with this list of events, so we void the indices
    %in the full index list to not go over them again in the why loop
    for j=1:length(idx2);indexlist(indexlist==idx2(j))=[];end;



    %list of height (above vent level) for events belonging to same
    %eruption
    h=MSVOLSO2_height(idx2)-MSVOLSO2_ventalt(idx2);
    %list of SO2 mass for events belonging to same
    %eruption
    m=MSVOLSO2_mass(idx2);
    mask=~isnan(h);

    if sum(m(mask))>0
        %to combine the events belonging to the same eruption, we assume
        %the following

        %combined SO2 mass is sum of masses
        comb_MSVOLSO2_mass=[comb_MSVOLSO2_mass;sum(m(mask))];
        %Combined SO2 height avl is mass-weighted average of heights
        comb_MSVOLSO2_height=[comb_MSVOLSO2_height;sum(m(mask).*h(mask))/sum(m(mask))];
        %year, VEI and latitudes are the same
        comb_MSVOLSO2_year=[comb_MSVOLSO2_year;MSVOLSO2_year(idx2(1))];
        comb_MSVOLSO2_VEI=[comb_MSVOLSO2_VEI;MSVOLSO2_VEI(idx2(1))];
        comb_MSVOLSO2_lat=[comb_MSVOLSO2_lat;MSVOLSO2_lat(idx2(1))];
    end

end
%the list of indices has now been voided, and comb_MSVOLSO2_parametername
%contains list of parameters where mutiple events belonging to one eruption
%have been combined together!

%%=========================================================================
%% 15) Load data for eruptions with high-confidence match before 1750
%%=========================================================================
%We've got our dataset complete!! Except it of course had lots of gaps...we
%could not find all eruption source parameters pre-satellite era for all
%eruptions, ice-core event with no eruption match need parameters, and even
%some satellite events have no height. So time to fill these gaps! To do
%so, we will get a bit of extra data (pre-processing)

%Second, we will use ice-core derived events that have been matched
%geochemically to a known eruption, which we consider to be a high
%confidence match. We already have those since 1750 in our CMIP7 dataset,
%but lots of good data before 1750! We've curated these from the literature
%in the following spreadsheet:

prehist=xlsread('geochemical_match_pre1750.xlsx');
prehist_mass=prehist(:,2);%Tg SO2 already
prehist_lat=prehist(:,8);
prehist_asym=prehist(:,3);
%below convert asymetry factor from NH/SH to (NH-SH)/(NH+SH)
prehist_asym(prehist_asym==-1 & prehist_lat<=0)=0;
prehist_asym(prehist_asym==-1 & prehist_lat>0)=inf;
prehist_asym=(1./(1+1./prehist_asym))-1./(prehist_asym+1);

prehist_matchconf=4*ones(size(prehist_mass));%assign high confidence flag
prehist_height=prehist(:,9);
prehist_valt=prehist(:,7);

%%=========================================================================
%% 16) Asymettry vs latitude analysis to fill latitude gap
%%=========================================================================
% Ok finally time to fill gaps! Let's start with latitude. The only cases
% for which we don't have latitude are ice-core derived events. But in
% general, we know the Greenland (NH) vs Antartica (SH) deposition
% asymetry. So let's find an empirical relationship between deposition
% asymetry and latitude to fill gap.

%The data we will use are: i) pre-satellite events from CMIP7; ii) pre-1750
%events with geochemical eruption match; El CHichon and Pinatubo (first two
%events in Sigl 2015 dataset) for which ice core deposition asymmetry is
%known. So we combine this data below
%in the pre-satellite event dataset, we also keep a single Laki phase to
%not have 10 datapoints for the same eruption. I don't seem to do that for
%Agung, I'm unsure why and it might be a mistake but only two phases so
%it's a small bias in terms of Agung weight on the final fit
masklaki=cmip7_eruN==12809 & cmip7_volcN==373010;
idx=find(masklaki);masklaki(idx(1))=(1<0);
asym_mass=[cmip7_mass(cmip7_year<1979 & ~masklaki);prehist_mass;S2015_mass(1:2)];
asym_asym=[cmip7_asym(cmip7_year<1979 & ~masklaki);prehist_asym;S2015_asym(1:2)];
asym_lat=[cmip7_lat(cmip7_year<1979 & ~masklaki);prehist_lat;15.14;17.36];
asym_matchconf=[cmip7_matchconf(cmip7_year<1979 & ~masklaki);prehist_matchconf;3;3];

%for plotting purpose, define a marker size proportional to SO2 mass
mksize=50+750*(asym_mass-min(asym_mass))/(max(asym_mass)-min(asym_mass));

%last, in the above, we only want to analyze events for which we are very
%confident we have the good match (match flag >=3), and we also remove
%outlier events with pure Antartica deposition but northern latitude
%volcano.
confmask=~isnan(asym_lat) & asym_matchconf>=3 & ~(asym_asym==-1 & asym_lat>0);

%plot asymettry vs latitude for these events for which we think data us
%good
figure
h1=scatter(asym_asym(confmask),asym_lat(confmask),mksize(confmask),'k^','MarkerFaceColor','#003f5c')
hold on
h2=scatter(asym_asym(asym_matchconf<3),asym_lat(asym_matchconf<3),mksize(asym_matchconf<3),'ks','MarkerFaceColor','#ffa600')
h3=scatter(asym_asym(asym_asym==-1 & asym_lat>0),asym_lat(asym_asym==-1 & asym_lat>0),mksize(asym_asym==-1 & asym_lat>0),'kp','MarkerFaceColor','r')

%then we just do a linear fit of that and it works decently.
myfittype = fittype("a*x",dependent="y",independent="x",coefficients=["a"])
x0 = [66];
[fit_asym2lat fit_asym2lat_stats]= fit(asym_asym(confmask),asym_lat(confmask),'poly1')
%NOTE: the asym2lat fit is used later to fill latitude from deposition
%asymetry

%add fit to the plot
h4=plot(-1:0.001:1,fit_asym2lat(-1:0.001:1),'-','Color','#003f5c','LineWidth',2)
fitcoef_asym2lat = coeffvalues(fit_asym2lat);
conf_asym2lat = predint(fit_asym2lat,-1:0.001:1,0.95);
h5=plot(-1:0.001:1,conf_asym2lat(:,1),'--','Color','#003f5c','LineWidth',1)
plot(-1:0.001:1,conf_asym2lat(:,2),'--','Color','#003f5c','LineWidth',1,'LineWidth',1)

%add the deposition asymetry-latitude relationship used in eVolv2k
h6=plot(-1:0.001:1,[-45 zeros(1,length(-1:0.001:1)-2) 45],':','Color','#bc5090','LineWidth',2)


%just for additional context, let's add model data to this plot from
%Marshall et al 2021
marsh=xlsread('marshalletal2021.xlsx');marsh_asym=marsh(:,4);marsh_lat=marsh(:,1);
h7=plot(marsh_asym,marsh_lat,'ko','MarkerSize',3)


%tiday the plot/add legend etc
legend([h1 h2 h3 h4 h5 h6 h7],'High-very high confidence matches (flag 3)','Low-medium confidence matches (flag 1-2)','Outlier high-confidence matches (flag 3)',...
    'Best linear fit (high-confidence matches, no outliers)','Fit prediction interval (95 conf.%)','Relationship used in eVolv2k','Interactive stratospheric aerosol model (Marshall et al., 2021)')

xlabel('Modified asymetry factor \alpha','FontSize',12)
ylabel('Eruption latitude \lambda (^oN)','FontSize',12)
ylim([-90 90])
xlim([-1.03 1.03])
text(-0.98,12,'Santa Maria 1902','Color','red')
text(-0.98,20,'El Chichon 1982','Color','red')


figure
plot(-1:0.001:1,fit_asym2lat(-1:0.001:1)-conf_asym2lat(:,1),'ro')
hold on
plot(-1:0.001:1,-fit_asym2lat(-1:0.001:1)+conf_asym2lat(:,2),'b-')
ylabel('latitude fit prediction uncertainty')
xlabel('Modified asymetry factor \alpha','FontSize',12)
legend('Latitude fit lower bound uncertainty','Latitude fit upper bound uncertainty')


%%=========================================================================
%% 17)SO2 height vs mass analysis to fill height gap
%%=========================================================================

% Onwards with gap filling! SO2 height is another one for which we miss
% lots of data, especially pre-satellite. The main parameter related to
% height would be SO2 mass. If you assume proportionality between SO2 mass
% and ash mass, and between ash mass and mass eruption rate, then there
% should be a power-law relationship between SO2 mass and height. The
% analysis below show that it's very noisy, but it's there.

%to analyze SO2 mass vs height, we first use any pre-satellite event from
%CMIP7 with good match confidence flag except Laki (effusive eruption)
% or any pre-1750 from our list of geochemically matched events with
%good match. So we combine this data below:
masklaki=cmip7_eruN==12809 & cmip7_volcN==373010;
hh_mass=[cmip7_mass(cmip7_year<1979 & ~masklaki);mean(cmip7_mass(masklaki));prehist_mass];
hh_height=[cmip7_height(cmip7_year<1979 & ~masklaki);mean(cmip7_height(masklaki));prehist_height];
hh_valt=[cmip7_ventalt(cmip7_year<1979 & ~masklaki);mean(cmip7_ventalt(masklaki));prehist_valt];
hh_lat=[cmip7_lat(cmip7_year<1979);prehist_lat];
hh_matchconf=[cmip7_matchconf(cmip7_year<1979 & ~masklaki);mean(cmip7_matchconf(masklaki));prehist_matchconf];
%we will also used satellite data (where eruption phases have been combined
%above), bit all processing done already

%we define marker size proportional to SO2 masses
hh_mksize=20+980*(hh_mass-min(comb_MSVOLSO2_mass))/(max(hh_mass)-min(comb_MSVOLSO2_mass));
MSVOLSO2_mksize=20+980*(comb_MSVOLSO2_mass-min(comb_MSVOLSO2_mass))/(max(hh_mass)-min(comb_MSVOLSO2_mass));

%we will plot all events pre-satellite but define mask based on whether
%the match confidence is low or high (>=3)
confmask=~isnan(hh_height) & hh_matchconf==3;
lowconfmask=~isnan(hh_height) & hh_matchconf<3;

figure
%plot satellite data
h1=scatter(comb_MSVOLSO2_mass,comb_MSVOLSO2_height,MSVOLSO2_mksize,'kp','MarkerFaceColor','#bc5090')
hold on
%plot high-confidence ice-core match
h2=scatter(hh_mass(confmask),hh_height(confmask)-hh_valt(confmask),hh_mksize(confmask),'ko','MarkerFaceColor','#003f5c')
%plot low-confidence match
h3=scatter(hh_mass(lowconfmask),hh_height(lowconfmask)-hh_valt(lowconfmask),hh_mksize(lowconfmask),'ks','MarkerFaceColor','#ffa600')
set(gca,'Xscale','log')

%do a power-law fit only using satellite data and high-confidence ice-core
%match
%NOTE: this fit will be used to fill gaps in SO2 injection height later.
%ALSO NOTE: all heights have been converted to above vent level, as that's
%the metric relevant for explosivity (e.g. you would expect zero mass
%eruption rate and zero so2 for a 5km a.s.l. plume coming from a 5km-high
%volcano.
[fit_mass2height fit_mass2height_stats]= fit([comb_MSVOLSO2_mass(~isnan(comb_MSVOLSO2_height));hh_mass(confmask)],[comb_MSVOLSO2_height(~isnan(comb_MSVOLSO2_height));hh_height(confmask)-hh_valt(confmask)],'power1','Weight',[comb_MSVOLSO2_mass(~isnan(comb_MSVOLSO2_height));hh_mass(confmask)])

%plot the fit including confidence interval
h4=plot(logspace(-4,2.5,1000),fit_mass2height(logspace(-4,2.5,1000)),'-','Color','#003f5c','LineWidth',2)

coef_mass2height = coeffvalues(fit_mass2height);
conf_mass2height = predint(fit_mass2height,logspace(-4,2.5,1000));
h5=plot(logspace(-4,2.3,1000),conf_mass2height(:,1),'--','Color','#003f5c','LineWidth',1)
plot(logspace(-4,2.3,1000),conf_mass2height(:,2),'--','Color','#003f5c','LineWidth',1)
h6=plot(logspace(-4,2.3,1000),1.66*fit_mass2height(logspace(-4,2.5,1000)),':','Color','#003f5c','LineWidth',1)
plot(logspace(-4,2.3,1000),0.33*fit_mass2height(logspace(-4,2.5,1000)),':','Color','#003f5c','LineWidth',1)


%add legend and tidy plot
legend([h1 h2 h3 h4 h5 h6],'Satellite data','High-confidence ice-core matches','Low-confidence ice-core matches',...
    'Best power-law fit (SO_2 mass-weighted, high-confidence  match+satellite)','Fit prediction interval (95% conf.)','Default relative uncertainty (66%) on height')

xlabel('Mass of SO_2 (Tg SO_2)','FontSize',12)
ylabel('SO_2 height (km a.v.l.)','FontSize',12)
xlim([10^(-4) 10^(2.3)])
ylim([0 46])


%==========================================================================
%% 18) Create additional variables and fill all gaps
%%=========================================================================

%Now that we have empirical relationships to fill missing latitude and
%heights, let's do the actuall gap filling for everything and create a few
%additional variables...
%For numerous variables, we also create fill flag indicating whether we
%have a value derived from actual information e.g. from source datasets or
%peer reviewed literature, or whether parameter was filled

%==========================================================================
%source data flag: create a source data variable more explicit than 1...5
%==========================================================================
cmip7_sourcedata=cell(size(cmip7_year));
cmip7_sourcedata(cmip7_source==5)={'MSVOLSO2L4'};
cmip7_sourcedata(cmip7_source==4)={'GVP'};
cmip7_sourcedata(cmip7_source==3)={'D4i'};
cmip7_sourcedata(cmip7_source==2)={'Sigl et al. (2015)'};
cmip7_sourcedata(cmip7_source==1)={'eVolv2k'};

%==========================================================================
%Fill eruption source parameter information if suggested match
%==========================================================================
%There might several reason why no ESP information so we fill this
emptyesp = cellfun(@isempty,cmip7_espinfo) | strcmp(cmip7_espinfo,'');
cmip7_espinfo(emptyesp & cmip7_source<=4)={'No suggested match or ESPs were not found.'};
cmip7_espinfo(emptyesp & cmip7_source==5)={'NA: All ESPs known from MSVOLSO2L4'};
cmip7_espinfo(strcmp(string(cmip7_espinfo),'') & cmip7_source<=4)={'No suggested match or ESPs were not found.'};

%==========================================================================
%SO2 mass stuff
%==========================================================================

%Create uncertainty fill flag and fill missing uncertainty
cmip7_dmass_fillflag=zeros(size(cmip7_year));%=0 if from source dataset
cmip7_dmass_fillflag(isnan(cmip7_dmass))=1;%=1 if value filled as per documentation
cmip7_dmass(isnan(cmip7_dmass))=cmip7_mass(isnan(cmip7_dmass))/3;%fill value for uncertainty=30% (Pinatubo uncertainty 15+/-5 Tg)

%Change the D4i masses as needed
%If D4i event tropical, do not change mass
%If D4i event NH, mass is 0.103/0.271 smaller (see companion paper)
cmip7_mass(cmip7_source==3 & cmip7_lat>=25)=cmip7_mass(cmip7_source==3 & cmip7_lat>=25)*0.103/0.271;
cmip7_dmass(cmip7_source==3 & cmip7_lat>=25)=cmip7_dmass(cmip7_source==3 & cmip7_lat>=25);

%For D4i events unidentified, there is a bit more work...
%If D4i event is unknown, using the ration of NH to tropical mass and the
%probability that an eruption depositing in Greenland is NH (0.6445) rather
%than tropical, the total mass is chosen as in Fang et al (2023) but we
%assume a single event at mid-latitude instead of one tropical and one NH event.
cmip7_mass(cmip7_source==3 & isnan(cmip7_lat))=cmip7_mass(cmip7_source==3 & isnan(cmip7_lat))*(1-0.6445*(1-0.103/0.271));
cmip7_dmass(cmip7_source==3 & isnan(cmip7_lat))=cmip7_dmass(cmip7_source==3 & isnan(cmip7_lat));

%sanity check no D4i event has SH latitude
if sum(cmip7_source==3 & cmip7_lat<=-25)>0
    error('Some D4i events have SH eruptions as source. Revise your match spreadsheet')
end
%==========================================================================
%Timing information
%==========================================================================

%cmip7_year: no gap to fill, source informed in cmip7_sourcedata_flag + match info
%cmip7_iceyear: document year in ice-core dataset before match for reference (e.g. future match change)

cmip7_month_fillflag=zeros(size(cmip7_year));%=0 if from source dataset or match
cmip7_month_fillflag(cmip7_source<=3 & isnan(cmip7_month))=1;%=1 if filled

%For evolv2k/sigl et al, we assume any missing year is same year as
%ice-core and any missing month is January
cmip7_year(cmip7_source~=3 & isnan(cmip7_year))=cmip7_iceyear(cmip7_source<=2 & isnan(cmip7_year));
cmip7_month(cmip7_source~=3 & isnan(cmip7_month))=1;

%for D4i, there is a fractional year because highres so we choose any
%missing year and month consistently with the fractional year provided
cmip7_year(cmip7_source==3 & isnan(cmip7_year))=floor(cmip7_iceyear(cmip7_source==3 & isnan(cmip7_year)));
cmip7_month(cmip7_source==3 & isnan(cmip7_month))=max(min(floor(12*(cmip7_iceyear(cmip7_source==3 & isnan(cmip7_month))-floor(cmip7_iceyear(cmip7_source==3 & isnan(cmip7_month))))),1),12);

%if day is missing, we assume 1
cmip7_day_fillflag=zeros(size(cmip7_year));%=0 if from source dataset or match
cmip7_day_fillflag(isnan(cmip7_day))=1;%=1 if filled
cmip7_day(isnan(cmip7_day))=1;


%==========================================================================
%location information
%==========================================================================
cmip7_location_fillflag=zeros(size(cmip7_year));%=0 if from source dataset or match
cmip7_location_fillflag(isnan(cmip7_lat))=1;%=1 if we had to fill latitude (and thus longitude)

%if longiture is missing, we assume 0
cmip7_lon(cmip7_location_fillflag==1)=0;%set to 0 and no suggested uncertainty; note longitude from -180 to 180

%if no latitude and not D4i, use deposition asymetry fit to fill latitude
cmip7_lat(cmip7_location_fillflag==1 & cmip7_source~=3)=fit_asym2lat(cmip7_asym(isnan(cmip7_lat) & cmip7_source~=3));%use best fit with asymetry
%if no latitude and D4i, assume 30N eruption latitude. ~65% of SO2 will go in NH EVA_H box and
%rest in tropical box.
cmip7_lat(cmip7_location_fillflag==1 & cmip7_source==3)=30;

cmip7_dlat=zeros(size(cmip7_lat));%no uncertainty on latitude if location known for certain or if match confidence 3
cmip7_dlat(cmip7_location_fillflag==1 | cmip7_matchconf==1 | cmip7_matchconf==2)=40;%latitude uncertainty = 40 as per our fit if unknown or low-medium confidence match

%if vent altitude is missing, fill value is mean of vent altitude value in
%dataset
defaultventalt=mean(unique(cmip7_ventalt),'omitnan');%default vent altitude = mean of know ones
cmip7_ventalt(cmip7_location_fillflag==1)=defaultventalt;

if sum(isnan(cmip7_ventalt) & cmip7_location_fillflag==0)>0
    x=input('Some volcanoes have unfilled vent altitude information. Continue? (y/n)',"s")
    if strcmp(x,'y')
        cmip7_ventalt(isnan(cmip7_ventalt) & cmip7_location_fillflag==0)=defaultventalt;
    else
        error('Then go fix your ESP spreadsheet :p')
    end
end

%==========================================================================
%Height information
%==========================================================================

cmip7_height_fillflag=zeros(size(cmip7_height));%=0 if height is known
cmip7_height_fillflag(isnan(cmip7_height))=1;%=1 otherwise

%If height is missing and the event is not from GVP, use our empirical SO2
%mass-SO2 height fit to fill.
% note that our fit for height is for km avl so need to add vent altitude
cmip7_height(isnan(cmip7_height) & cmip7_source~=4)=fit_mass2height(cmip7_mass(isnan(cmip7_height) & cmip7_source~=4))+cmip7_ventalt(isnan(cmip7_height) & cmip7_source~=4);

%If event is from GVP and has no height constraint, no SO2 mass->can't use
%the fit. So we make arbitrary assumptions
cmip7_height(isnan(cmip7_height) & cmip7_source==4 & cmip7_VEI==4 & cmip7_latband==2)=16;
cmip7_height(isnan(cmip7_height) & cmip7_source==4 & cmip7_VEI==4 & cmip7_latband~=2)=11;
cmip7_height(isnan(cmip7_height) & cmip7_source==4 & cmip7_VEI==5)=18;


%for known eruption or high-confidence match, if satelitte era,  assume
%uncertainty of 38% on avl height (based on IVESPA); arbitrarly assume 50
%if before satellite era
height_mask=isnan(cmip7_height) | cmip7_matchconf==1 | cmip7_matchconf==2;
cmip7_dheight(~height_mask & cmip7_year>=1979)=(cmip7_height(~height_mask & cmip7_year>=1979)-cmip7_ventalt(~height_mask & cmip7_year>=1979))*0.38;
cmip7_dheight(~height_mask & cmip7_year<1979)=(cmip7_height(~height_mask & cmip7_year<1979)-cmip7_ventalt(~height_mask & cmip7_year<1979))*0.5;
%for other eruptions, assume 66% uncertainty on height avl and add 100% vent uncertainty
cmip7_dheight(height_mask)=(((cmip7_height(height_mask)-cmip7_ventalt(height_mask))*0.66).^2+cmip7_ventalt(height_mask).^2).^0.5;

%==========================================================================
%Asymetry factor
%==========================================================================
%cmip7_asym, nothing to change, NaN for satellite date

%==========================================================================
%Links to GVP database and VEI
%==========================================================================
%cmip7_volcname, cmip7_VEI, cmip7_eruN and CMIP7_volcN
%Nothing to do, and always from proposed match if there is a match
%Note that the eruption and volcano number from GVP is not linked in the MSVOLSO2dataset

%correct empty volc names for consistency
emptyvolcname = cellfun(@isempty,cmip7_volcname) | strcmp(cmip7_volcname,'NaN');
cmip7_volcname(emptyvolcname)={'Unknown eruption'};

%==========================================================================
%Emission depth
%==========================================================================
%create emission depth to specify over which charactteristic heigh to spread
% emission around SO2 height. We use the gaussian width constrained from 3D
% plume model simulations in Aubry et al. GRL 2019
cmip7_depth=0.108*(cmip7_height-cmip7_ventalt);


%==========================================================================
%% 19)Sanity checks, sort chronologically and save CMIP7 emission as xlsx
%%=========================================================================

%Almost almost there! Let's just do some sanity check as it's highly likely
% something fishy will have happened in the previous 1300 lines, then we
% save!
%the emis_sanchecks_num function takes a variable and check it's within the
%provided range. A;; ranges should be obvious unless I comment

emis_sanchecks_num(cmip7_year,1,[1750 2024])
emis_sanchecks_num(cmip7_month,1,[1 12])
emis_sanchecks_num(cmip7_day,1,[1 31])
emis_sanchecks_num(cmip7_mass,1,[0 100])%no mass above 100 Tg SO2 for 1750-present
emis_sanchecks_num(cmip7_dmass,1,[0 100])
emis_sanchecks_num(cmip7_height,1,[0 50])%height above 50km asl is a no for SO2, and yes this includes HTHH 2022
emis_sanchecks_num(cmip7_dheight,1,[0 20])
emis_sanchecks_num(cmip7_lat,1,[-90 90])
emis_sanchecks_num(cmip7_dlat,1,[0 90])
emis_sanchecks_num(cmip7_lon,1,[-180 180])
emis_sanchecks_num(cmip7_VEI,0,[0 7])
emis_sanchecks_num(cmip7_ventalt,1,[-0.5 6])%-0.5 for shallow submarine vent if any
emis_sanchecks_num(cmip7_asym,0,[-1 1])%deposition asymetry is -1 to 1 with our definition (NH-SH)/(NH+SH)
emis_sanchecks_num(cmip7_iceyear,0,[1750 1980])
emis_sanchecks_num(cmip7_matchconf,0,[1 4])
emis_sanchecks_num(cmip7_depth,1,[0 5])
emis_sanchecks_num(cmip7_month_fillflag,1,[0 1])
emis_sanchecks_num(cmip7_day_fillflag,1,[0 1])
emis_sanchecks_num(cmip7_location_fillflag,1,[0 1])
emis_sanchecks_num(cmip7_dmass_fillflag,1,[0 1])

%Then we check that text info on ESO choice or match are always there, as
%well as source dataset
if sum(isempty(cmip7_espinfo))>0;error('Missing ESP info field');end;
if sum(isempty(cmip7_matchinfo))>0;error('Missing match info field');end;
srccheck=~strcmp(cmip7_sourcedata,'MSVOLSO2L4') & ...
    ~strcmp(cmip7_sourcedata,'eVolv2k') & ...
    ~strcmp(cmip7_sourcedata,'D4i') & ...
    ~strcmp(cmip7_sourcedata,'Sigl et al. (2015)') & ...
    ~strcmp(cmip7_sourcedata,'GVP');
if sum(srccheck)>0;error('wrong or missing source dataset name');end;

%then we order events chronologically
eventtime=cmip7_year+cmip7_month/12+cmip7_day/12/31;
[a sortedind]=sort(eventtime);

%put it all in a big table
cmip7_all_final=[cmip7_volcname(sortedind),...
    num2cell(cmip7_year(sortedind)),...
    num2cell(cmip7_month(sortedind)),...
    num2cell(cmip7_month_fillflag(sortedind)),...
    num2cell(cmip7_day(sortedind)),...
    num2cell(cmip7_day_fillflag(sortedind)),...
    num2cell(cmip7_lon(sortedind)),...
    num2cell(cmip7_lat(sortedind)),...
    num2cell(cmip7_dlat(sortedind)),...
    num2cell(cmip7_location_fillflag(sortedind)),...
    num2cell(cmip7_mass(sortedind)),...
    num2cell(cmip7_dmass(sortedind)),...
    num2cell(cmip7_dmass_fillflag(sortedind)),...
    num2cell(cmip7_height(sortedind)),...
    num2cell(cmip7_height_fillflag(sortedind)),...
    num2cell(cmip7_dheight(sortedind)),...
    num2cell(cmip7_depth(sortedind)),...
    cmip7_sourcedata(sortedind),...
    cmip7_matchinfo(sortedind),...
    num2cell(cmip7_matchconf(sortedind)),...
    cmip7_espinfo(sortedind),...
    num2cell(cmip7_VEI(sortedind)),...
    num2cell(cmip7_eruN(sortedind)),...
    num2cell(cmip7_volcN(sortedind)),...
    num2cell(cmip7_ventalt(sortedind)),...
    num2cell(cmip7_asym(sortedind)),...
    num2cell(cmip7_iceyear(sortedind))];

%add headers
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
    'Match confidence (1=low, 2=medium, 3=high)',...
    'Information on eruption source parameters (date and location)',...
    'Volcanic Explosivity Index',...
    'GVP eruption number',...
    'GVP volcano number',...
    'Vent altitude (km a.s.l.)',...
    'Ice-core sulfur deposition asymmetry factor (from -1=only Antartica to 1=only Greenland)',...
    'Year of deposition in ice-core'};

%combine headers and data and save
cmip7_all_final=[cmip7_all_headers;cmip7_all_final];
writecell(cmip7_all_final,'prelim_CMIP7_volcano_S_emissions_v2_2_1.xlsx','FileType','spreadsheet','WriteMode','replacefile')

%And voila!! Only took 1383 lines
%IMPORTANT NOTE though: if you followed this script carefully, we entered
%very arbitrary SO2 masses for GVP-derived events...the dataset saved is
%not quite final, these masses get finalized in the SAOP script as we tried
%to be clever-ish when attributing these masses :)
