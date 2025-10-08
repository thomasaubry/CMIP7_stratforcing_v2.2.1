% 	Written by Thomas J. Aubry, July 2025.
% 	Department of Earth Sciences, University of Oxford, UK
%   E-mail: thomas.aubry@earth.ox.ac.uk
% 	Please cite the corresponding papers if you use or modify this script,
%   i.e. both the original EVA_H paper (Aubry et al., JGR 2020) and the CMIP7
%   dataset paper documenting version 2 of EVA_H



%This functions calculate the time derivative of the mass of sulfate in
%each of the 8 boxes of the model (dydt) given the mass of sulfate y at the
%time t, the volcanic SO2 injections (inmass and intime),and the model
%parameters (coef and backinj)
function dydt = eightboxequations_CMIP7(t,y,inmass,intime,inheight,totmass,coef)
%t is the time (real number) in month after January 1st of a user-chosen
%year
%dydt is the 8x1 vector containing the time derivative of the sulfate mass
%in each box (in Tg S/month)
%y is the 8x1 vector containing the sulfate mass in each box (in Tg S)
%inmass is a 8xNeru matrix containing the SO2 mass injected in the 8 boxes by
%the Neru eruptions (in Tg S)
%intime is a 1xNeru matrix containing the corresponding times of injections
%(same unit as t)
%inheight and totmass are a Neru x1 vectors containing the SO2 injection
%height in km a.s.l. and total (not stratospheric only) injected SO2 mass
%in Tg SO2
%coef is a 53x1 vector containing values for all model parameters (cf
%parameterfile.m for units)


%==========================================================================
%1) "unwrap" vector containing model parameters
%==========================================================================
%Separate model parameters into different, more explicit variables. See
%parameterfile.m for details on each parameters and how they were combined
%into a single vector.

A=coef(1);% SAOD = A x MSO4tot
tauprod=max(coef(2).*(1+coef(3).*inheight').*sum(inmass,1).^coef(4),0.1);% production timescales using the relationship introduced in the CMIP7 paper
tauloss=coef(10:17);%loss timescales
taumixm=coef(18:23);%mixing timescales
amix=coef(24:29);%amplitude of mixing seasonal cycle
smix=coef(30:35);%peak mixing season
tauowmm=coef(36:41);%one-way mixing (OWM) timescales
aowm=coef(42:47);%OWM seasonal cycle amplitude
sowm=coef(48:53);%OWM peak mixing season
%coef=[A;tauprod;tauloss;taumixm;amix;smix;tauowmm;aowm;sowm;backinj];

%Calculate value of mixing timescales given the month of the year and
%parameters defining seasonal cycles
taumix=taumixm.*(1+amix.*cos((t-smix)*pi/6));
tauowm=tauowmm.*(1+aowm.*cos((t-sowm)*pi/6));

%==========================================================================
%2) Calculate the time derivative of sulfate mass
%==========================================================================

%This time derivative is calculated as dydt=C * y+ production + background
%The matrix C and term C*y encompasse all transport terms
%The term production correspond to new sulfate producted from SO2
%injections
%The term background correspond to background sulfate injections


%==========================================================================
%2.a) Build the matrix C
%==========================================================================

%Initiate the matrix with all terms being 0. We will specify non-zero terms
%hereafter.
C=zeros(8);

%Non zero coefficient are specified row by row, i.e. box by box. Detailed
%explainations are provided for select examples

%BOX 1
i=1;%number of the matrix row corresponding to box number

C(i,1)=-(1/tauloss(i))-(1/taumix(1));
%All boxes have a loss term equal to their mass divided by the loss
%timescale of the box, so that C(i,i) have a term -(1/tauloss(i)) as in
%line above

C(i,2)=+(1/tauowm(1))+(1/taumix(1)) ;
%Flux related to mixing is related to mass difference between two boxes. In
%the case of box 1, mixing result in a mass flux (y(2)-y(1))/taumix(1)
%explaining terms -(1/taumix(1)) in C(1,1) and +(1/taumix(1)) in C(1,2).

%BOX 2
i=2;
C(i,1)=+1/taumix(1);
C(i,2)=-(1/tauloss(i))-(1/taumix(1)+1/taumix(2)+1/tauowm(1)+1/tauowm(2)) ;
%In addition to mixing, one-way mixing terms transport mass from tropical
%to extratropical box, with the flux being proportional to the mass of the
%tropical box. Thus, box 2, which is tropical, lose a flux -y(2)/tauowm(1) to
%box 1 and y(2)/tauowm(2) to box 3. This explains terms
%-1/tauowm(1)-1/tauowm(2) for C(2,2) and, for example, the term
%+(1/tauowm(1)) for C(1,2)

C(i,3)=+1/taumix(2) ;

%BOX 3
i=3;
C(i,2)=1/taumix(2)+1/tauowm(2) ;
C(i,3)=-(1/tauloss(i))-1/taumix(2) ;

%BOX 4
i=4;
C(i,1)=1/tauloss(1);
%The mass of aerosol lost by a box by settling (loss flux) correspond to a
%mass gain for the model box immediatly below. Thus, the term -1/tauloss(1)
%for C(1,1) correspond to a term +1/tauloss(1) for C(4,1) as box 4 is below
%box 1.

C(i,4)=-(1/tauloss(i))-1/taumix(3) ;
C(i,5)=1/taumix(3)+1/tauowm(3) ;

%BOX 5
i=5;
C(i,2)=1/tauloss(2) ;
C(i,4)=1/taumix(3) ;
C(i,5)= -(1/tauloss(i))-(sum(1./taumix(3:6))+sum(1./tauowm(3:6)));
%Mixing fluxes are only between boxes that are neighbour and at the same
%altitude, except for box 5 (tropical lower stratosphere) for which we
%allow mixing with boxes 7 and 8 (lowermost extratropical stratosphere.
%This explains, for example, the terms -1/taumix(5) in C(5,5) and
%+1/taumix(5) in C(5,7)
C(i,6)=1/taumix(4) ;
C(i,7)=1/taumix(5) ;
C(i,8)=1/taumix(6) ;

%BOX 6
i=6;
C(i,3)=1/tauloss(3) ;
C(i,5)=1/taumix(4)+1/tauowm(4) ;
C(i,6)= -(1/tauloss(i))-1/taumix(4);

%BOX 7
i=7;

C(i,4)= 1/tauloss(4);
C(i,5)=1/taumix(5)+1/tauowm(5) ;
C(i,7)= -(1/tauloss(i))-1/taumix(5);

%BOX 8
i=8;
C(i,6)= 1/tauloss(6);
C(i,5)= 1/taumix(6)+1/tauowm(6);
C(i,8)=-(1/tauloss(i))-1/taumix(6) ;



%==========================================================================
%2.b) Calculate the mass of SO2 in each box
%==========================================================================

%filter out all eruptions which occured after the current timestep as these
%do not matter to calculate the mass of SO2 at current timestep.

mask=intime<=t;
inmassint=inmass(:,mask);
intimeint=intime(mask);
tauprodint=tauprod(mask);

%if no eruption has occured yet, the SO2 mass in all boxes is 0
if sum(mask)==0
    so2mass=0;
else
    % Else, the SO2 mass is calculated by assuming that the mass of SO2
    % declines exponentially after each eruption, with an e-folding time equal
    % to the production timescale for each box
    so2mass=inmassint.*exp(-repmat(t-intimeint',[8 1])./repmat(tauprodint,[8 1]));
end

%Note that there is no SO2 transport in our model. Some of the transport
%that would occur is accounted for by the latitudinal distribution of the
%mass of SO2 injected by an eruption (cf. so2injection_8boxes.m and
%companion paper).

%==========================================================================
%2.c) Calculate the time derivative of sulfate mass by summing the 3
%components
%==========================================================================


dydt = C*y+sum(so2mass./repmat(tauprodint,[8 1]),2);
% Note: no more background injection in EVA_H v2
% Sulfate production from volcanic SO2 is simply the SO2 mass calculated in
% 2.b divided by sulfate production timescales, so2mass./tauprod
% Transport terms are the only ones dependent on sulfate mass at current
% timestep, and are the product of the matrix C calculated in 2.a with the
% mass of sulfate in the 8 boxes.


