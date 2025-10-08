% 	Written by Thomas J. Aubry, July 2025.
% 	Department of Earth Sciences, University of Oxford, UK
%   E-mail: thomas.aubry@earth.ox.ac.uk
% 	Please cite the corresponding papers if you use or modify this script,
%   i.e. both the original EVA_H paper (Aubry et al., JGR 2020) and the CMIP7
%   dataset paper documenting version 2 of EVA_H

function [gmsaod, saod, reff, ext, ssa, asy, sad, vd, nd]=postproc_glwl_v3(ext_wl1,wl1,ext_wl2,wl2,wl_req,lat,miename)
%see CMIP paper for variable output name definitions
%This functions returns key aerosol (optical) properties at wavelength wl_req (in um)
% from extinction fields (ext_wl1,ext_wl2, /km) at two different wavelengths(wl1,wl2
%um), using the Mie look-up table miename. lat is latitude dimension of
%fields. In the CMIP7 script, we use it to derive full aerosol property
%fields from the EVA_H v2 Mie look-up table and the 525 and 1020nm
%extinction from GloSSAC

%==========================================================================
%1) Defining some useful parameters
%==========================================================================
%density of sulfate aerosol in kg/m3 at 215K and 0.75wt%
%(https://doi.org/10.1038/s41598-019-52089-6)
rho_aer=(-0.4845-0.7074*0.75)*215+1186.1+ 621.4*0.75+ 573.54*0.75*0.75;
%in g/um3 this is:
rho_aer=rho_aer*1000/(1000000^3);
%Avogadro's number (in /mol)
n_av=6.02214*10^23;
%sulfuric acid concentration in aerosol (mass fraction)
wf_h2so4=0.75;
%molar mass h2so4 in g/mol
mm_h2so4=98.079;
%the number concentration nd of h2so4 in (molecules h2so4/cm3 air), given the
%volume density of aerosol vd in um3 aerosol/cm3 air is then
%nd=(vd*rho_aer*0.75)*n_av/mm_h2so4;

%==========================================================================
%2) Read Mie file amd calculate effective radius to match ratio of ext at provided
%wavelengths
%==========================================================================

%Read Mie look-up tables:
ncid = netcdf.open(miename);

%a) Read effective radius and wavelength grid
reffgrid_mie = varfromname(miename,'reff',-999);
wlgrid_mie = varfromname(miename,'wl',-999);


%b) Read calculated parameters, which are 2D array, with one dimension for
%effective radius and the other one for wavelength.
extrat_mie = varfromname(miename,'extrat',-999);%ratio of extinction to extinction at 550nm (EXT)
ssa_mie = varfromname(miename,'ssa',-999);%single scattering albedo (SSA)
asy_mie = varfromname(miename,'asy',-999);%scattering asymmetry factor (ASY)

sad_mie = varfromname(miename,'sad_unit',-999);%scattering asymmetry factor (ASY)
vad_mie = varfromname(miename,'vol_unit',-999);%scattering asymmetry factor (ASY)
ecs_mie = varfromname(miename,'ecs',-999);%effective cross section

%Mie files are not always perfect so a bit of pre-processing. Note when
% values outside physical range occur, they are only just outside this
% range and are likely interpolation or rounding errors

%Filter any effective radius with negative values of positively defined
%ecs, vad, sad and reff
maskmie=ecs_mie>0 & vad_mie>0 & sad_mie>0 & reffgrid_mie>=0.06;
reffgrid_mie=reffgrid_mie(maskmie);
extrat_mie=extrat_mie(maskmie,:);
ssa_mie=ssa_mie(maskmie,:);
asy_mie=asy_mie(maskmie,:);
sad_mie=sad_mie(maskmie);
vad_mie=vad_mie(maskmie);
ecs_mie=ecs_mie(maskmie);

%Mie file not always perfect so assume negative values of positively
%defined variables are zero, and ssa above 1 are 1.
extrat_mie(extrat_mie<0)=0;
ssa_mie(ssa_mie<0)=0;
ssa_mie(ssa_mie>1)=1;
asy_mie(asy_mie<0)=0;

%==========================================================================
%3) Calculate effective radius
%==========================================================================

%make mess from reff and wl value
[Xwl,Yreff] = meshgrid(wlgrid_mie,reffgrid_mie);
%calculate extinction ratio at the two inputted wavelengths from Mie
%look-up table
mie_ratio_wl2_to_wl1=interp2(Xwl,Yreff,extrat_mie,wl2,reffgrid_mie)./interp2(Xwl,Yreff,extrat_mie,wl1,reffgrid_mie);
%calculate ratio of inputted extinction
ratio_wl2_to_wl1=ext_wl2./ext_wl1;

%calculate reff values required to match Mie-calculated and inputted
%extunction ratio. Only use Reff values smaller than value for which peak
%extinction ratio is used in the file, otherwise two possible solution.
reff=NaN(size(ext_wl1));
[maxrat ind]=max(mie_ratio_wl2_to_wl1);
reff(:)=interp1(mie_ratio_wl2_to_wl1(1:ind),reffgrid_mie(1:ind),ratio_wl2_to_wl1(:),'spline');
%ensure interpolation does not result in values of reff extrapolated
%outside range
reff(reff(:)<min(reffgrid_mie))=min(reffgrid_mie);
reff(reff(:)>reffgrid_mie(ind))=reffgrid_mie(ind);

%==========================================================================
%4) Calculate all properties at requested wavelengths
%==========================================================================


%preallocate memory for calculating EXT, SSA and ASY, vd, nd and sad
ext=NaN(size(ext_wl1,1),size(ext_wl1,2),size(ext_wl1,3),length(wl_req));
ssa=NaN(size(ext));
asy=NaN(size(ext));
sad=NaN(size(squeeze(ext(:,:,:,1))));
vd=NaN(size(squeeze(ext(:,:,:,1))));
nd=NaN(size(squeeze(ext(:,:,:,1))));

% Loop through latitude, altitude and wavelength to calculate them. All
%calculations are done by linearly interpolating the Mie lookup tables at
%the requested wavelength and the effective radius outputted by the model.
for ilat=1:size(ext_wl1,2)
    for ialt=1:size(ext_wl1,3)
        for iwl=1:length(wl_req)
            %ignore points in time where the extinction at 525nm or
            %effective radius are NaNs
            mask=~isnan(squeeze(ext_wl1(:,ilat,ialt))) & ~isnan(squeeze(reff(:,ilat,ialt)));
            %Calculate the raio of extinction at desired wavelength to
            %extinction at wl1
            ratio=interp2(Xwl,Yreff,extrat_mie,wl_req(iwl),reff(mask,ilat,ialt))./interp2(Xwl,Yreff,extrat_mie,wl1,reff(mask,ilat,ialt));
            %multiple above ratio by extinction at 525nm
            ext(mask,ilat,ialt,iwl)=ext_wl1(mask,ilat,ialt).*ratio;
            %calculate SSA and ASY
            ssa(mask,ilat,ialt,iwl)=interp2(Xwl,Yreff,ssa_mie,wl_req(iwl),reff(mask,ilat,ialt));
            asy(mask,ilat,ialt,iwl)=interp2(Xwl,Yreff,asy_mie,wl_req(iwl),reff(mask,ilat,ialt));


        end
        %cross section in mum^2
        i550=find(wl_req==0.55);
        % ecs=interp1(reffgrid_mie,ecs_mie,reff(mask,ilat,ialt));
        %numbconc=1000.*squeeze(ext(mask,ilat,ialt,i550))./ecs; % Sujan Khanal: number concentration (No, per cc) as: 1000*ext550nm/ecs. Here, ext550nm is in per km, ecs is the effective cross section from the LUT and 1000 is a factor to make the units work
        %ext above is in km so consistent with Sujan's formula
        %sad and vd from mie files need to be multiplied by number
        %concentration
        sad(mask,ilat,ialt)=interp1(reffgrid_mie,sad_mie./ecs_mie,reff(mask,ilat,ialt)).*1000.*squeeze(ext(mask,ilat,ialt,i550));
        vd(mask,ilat,ialt)=interp1(reffgrid_mie,vad_mie./ecs_mie,reff(mask,ilat,ialt)).*1000.*squeeze(ext(mask,ilat,ialt,i550));
        nd(mask,ilat,ialt)=(vd(mask,ilat,ialt)*rho_aer*0.75)*n_av/mm_h2so4;

    end
end


%Calculate stratospheric aerosol optical depth. These are simply the sum of
%extinction along vertical dimension, multiplied by 0.5 because the
%vertical grid is regurlarly spaced by 0.5km. All tropospheric values are
%NaNs in the shape functions and thus in ext.
saod=squeeze(sum(ext,3,'omitmissing')*0.5);

%Calculate weights (cosinus(latitude)) for calculating global mean
latweight=permute(repmat(cosd(lat),[1 length(wl_req) size(ext,1)]),[3 1 2]);
latweight=latweight/squeeze(sum(latweight(1,:,1),2));
%Calculate global mean SAOD
gmsaod=squeeze(sum(saod.*latweight,2,'omitmissing'));

end