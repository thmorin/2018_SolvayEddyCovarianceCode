%Olentangy Wetlands Research Park post processing code
%Written by Tim Morin
%Last edited Jan. 2014

%cd C:/Users/morin.37/Desktop/Matlab_Tim/ORWRP/PostProcessing/
%cd /home/morin.37/poolA/wetland_data/Matlab_Tim/ORWRP/PostProcessing/
% clear all;clc;
%%
%%%%%%%%%%%%%%%%%%%%%%USER INPUTS%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Input season to run as array values. Input one season. If multiple seasons
%are required, please run Matlab instances in parallel as this will greatly
%improve run times.
Season=1;

disp(['RUNNING: Post processing for season: ' num2str(Season)]);
%Iterations to run on neural network
NNRuns=1000;
%NNRuns=8;

footprintdirectory = './';

%Where 30 minute processed data is stored at
datadirectory = '../';

%Where subroutines are located. Will add subfolders recursively
codedirectory = './../../PostProcessSubs/';

%Geographic location of tower
lat=40.02; 
lon=-83.02;
GMToffset=-5;

%Height of the forest surrounding the tower
CanopyHeight=10;

SeasonDates=[
    2011 001,   2011 118 %Season 1: Winter 2010-2011. Only has from January 2011 forward
    2011 119,   2011 282 %Season 2: Summer 2011
    2011 283,   2012 106 %Season 3: Winter 2011-2012
    2012 107,   2012 281 %Season 4: Summer 2012 (Tower height change enacted on April 23)
    2012 282,   2013 122 %Season 5: Winter 2012-2013 
    2013 123,   2013 284 %Season 6: Summer 2013
    2013 285,   2014 125 %Season 7: Winter 2013-2014 (LI-7700 lowered on November 07 DOY 311)
    2014 126,   2014 277 %Season 8: Summer 2014 (LI-7700 raised on April 09 DOY 99)
    2014 278,   2014 300 %Season 9: Winter 2014-2014 (LI-7700 lowered on October 17 DOY 290)
    ];

%%%%%%%%%%%%%%%%%%%%%%%END USER INPUTS%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%
%%%%%%%%%%%%%%%%%%%%%%%DATA VALIDATION%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if exist(footprintdirectory,'dir')~=7 && exist(footprintdirectory,'file')~=2
    error('footprintdirectory does not exist. Please check this over and retry');
end
if exist(datadirectory,'dir')~=7
    error('datadirectory does not exist. Please check this over and retry');
end
if exist(codedirectory,'dir')~=7
    error('codedirectory does not exist. Please check this over and retry');
end
%%%%%%%%%%%%%%%%%%%%%%%END DATA VALIDATION%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%
%%%%%%%%%%%%%%%%%%%%%%%BEGIN PROCESSING%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
addpath(genpath(codedirectory));

if mod(Season,1)==0
    SeasonStatus='Summer';
else
    SeasonStatus='Winter';
end

start_year=SeasonDates(Season,1);
start_DOY=SeasonDates(Season,2);
end_year=SeasonDates(Season,3);
end_DOY=SeasonDates(Season,4);

%Minimum
MinUstar=0.2;
if strcmp(SeasonStatus,'Winter')
    MinUstar=0.15;
end

%Implemented for ORWRP. Will be ignored if z1 and z2 are input
z1=??;
if exist('z1','var')~=1 && exist('z2','var')~=1
    if Season==4 && Season<7 %Summer heights - no DS2
      z1=15;
      z2=17;
      NumWindMeasures=2;
    elseif Season==7 || Season==9 %Lowered winter heights
      z1=12;
      z2=17;
      NumWindMeasures=3;
    elseif Season==8 %Summer heights - with DS2 (which stopped working right away)
      z1=15;
      z2=17;
      NumWindMeasures=3;
    else %Really early heights
      z1=9.6;
      z2=11;
      CanopyHeight=7.5;
      NumWindMeasures=2;
    end
elseif exist('z2','var')~=1
    NumWindMeasures=1;
else
    NumWindMeasures=2;
end

[Dspk]=PostProcDespikeLimits(SeasonStatus);

%%
[extra_vars,DOY,DECDAY,DECYear,tair,tair_gf,pressure,pressure_gf,Fc,Fc_gf,LE,~,SWinbar,SWinbar_gf,SWoutbar,SWoutbar_gf,net_sw,net_sw_gf,LWinbar,LWinbar_gf,LWoutbar,LWoutbar_gf,net_lw,net_lw_gf,Q,Q_gf,...
Total_PAR,RNET,RNET_gf,rH,rH_gf,vpd,vpd_gf,L,L_gf,ww,ww_gf,SoilTemp,SoilTemp_gf,vv_1,vv_1_gf,H,~,wind_speed_1,wind_speed_1_gf,wind_dir_1,wind_dir_1_gf,ustar_1,ustar_1_gf,...
	   wind_speed_2,wind_speed_2_gf,wind_dir_2,wind_dir_2_gf,ustar_2,ustar_2_gf,vv_2,vv_2_gf,Spectral_correction,C1,C2,C]=...
    ORW_build_arrays(datadirectory,NumWindMeasures,Dspk,start_year,start_DOY,end_year,end_DOY,'FCH4','Spectral_correctionLI7700');

extra_vars.FCH4=extra_vars.FCH4./extra_vars.Spectral_correctionLI7700;
%%
%
if Season==2
%%%%%%Special case run for ORWRP to gapfill soil temperature before installation was complete%%%%%%%
    Jan01=2929;
%     load('/home/morin.37/poolA/wetland_data/Matlab_Liel/ModelingSoilRespirationCode/ModelingSoilRespirationCode/PostProcess/ModelWetFunc/SoilTemp2011_2013.mat','ST_gf');
    load([datadirectory '/PostProcessing/SoilTemp2011_2013.mat'],'ST_gf');
    SoilTemp_gf=ST_gf(Jan01+1+(start_DOY-1)*48:Jan01+end_DOY*48);
end
    
%Calculate marker for day vs night times
if sum(~isnan(Total_PAR))>0.85*length(Total_PAR)
    daynight=Daynightcalc(DECYear,Total_PAR,lat,lon,GMToffset);
else
    disp('Warning: PAR not available. Using SW radiation to generate daynight');
    daynight=Daynightcalc(DECYear,SWinbar_gf,lat,lon,GMToffset);
end

%%Default to RMY
%
GoodData=true(size(ustar_1));
if Season==4
    GoodData(1:(115-107)*48)=false;
elseif Season==7
    GoodData(1:1250)   = false; %Tower change 
    Fc(1250:1296)      = nan;
    GoodData(8622:end) = false;
elseif Season==9
    GoodData(1:2)=false;%FIXME when you run this season
end

%%
%Storage
if Season==8
    C1(1:2028)=nan;C2(1:2028)=nan;
 else
   C1=nan(size(C1));C2=C1;
end
C1=PpmToumolPerM3(C1,tair_gf,pressure_gf);
C2=PpmToumolPerM3(C2,tair_gf,pressure_gf);
C =C*1000;
[dCdt,Stor]=CarbonStorageFlux([C1 C2 C],[2 5 z1],1800);
Fc_true=nansum([Fc dCdt],2);Fc_true(isnan(Fc))=nan;

%%
%Roughness length

disp('Attempting rougness length calculation with data from lower anemometer')
[d,z0,n]=RoughnessLength(z1,CanopyHeight,L(GoodData),ustar_1(GoodData),wind_speed_1(GoodData),daynight(GoodData));
if isnan(d) || d<0.5*CanopyHeight
    disp('Defaulting to d=7 m and z0=2 m')
    d=7; %displacement height (CHANGE VERO) 2/3
    z0=.1; %roughness lenght
end
%%
%u* filter
UstLimits=nan;
if exist('z2','var')==1
    if sum(~isnan(ustar_2))>0.5*length(ustar_2)
        [UstLimits,USTflag]=ustarfiltdataprocess(Fc_true,ustar_2,tair,daynight,MinUstar);
    end
end
if isnan(UstLimits) && exist('z1','var')==1
    if sum(~isnan(ustar_1))>.5*length(ustar_1) && Season~=6 && Season~=2
        disp('Unable to process u* accurately with info from upper anemometer, attemping with lower');
        [UstLimits,USTflag]=ustarfiltdataprocess(Fc_true,ustar_1,tair,daynight,MinUstar);
    end
end
if isnan(UstLimits) && exist('z1','var')==1 && exist('z2','var')==1
    disp('Unable to process u* accurately with info from either anemometer, defaulting to combo of two');
    ust_all=ustar_1;ust_all(isnan(ust_all))=ustar_2(isnan(ust_all));
    [UstLimits,USTflag]=ustarfiltdataprocess(Fc_true,ust_all,tair,daynight,MinUstar);
end
if isnan(UstLimits)
    disp(['Defaulting u* value to minimum (' num2str(MinUstar) ')']);
    UstLimits=MinUstar;
    USTflag=ustar_1<UstLimits;
end
ust_killed=sum(USTflag)/length(USTflag);
disp(['u* filter set at: ' num2str(UstLimits) ' m/s']);
disp(['u* eliminates ' num2str(floor(ust_killed*100)) '% of data']);

try
if matlabpool('size')==0
    matlabpool(8)
end

%%
%Calculates where air packet originates
if exist(footprintdirectory,'dir')~=7
    load(footprintdirectory);
else
    VegYear=floor(DECYear(1));
    load([footprintdirectory num2str(VegYear) 'ORWRP_veg_update.mat'])
end
FX=Fx1(1,:)';
FY=Fy1(:,1);

PlotVec=zeros(size(ustar_1));

%%
footprint_full_temp_1=[];
footprint_full_temp_2=[];
if Season==4
    footprint_full_temp_1=FP_process_ORWRP(ustar_1_gf(~GoodData),vv_1_gf(~GoodData),L_gf(~GoodData),wind_dir_1_gf(~GoodData),1.5,9.6,7.5,FX,FY,WetlandMap,PlotVec(~GoodData));
    footprint_full_temp_2=[];
elseif Season==7 || Season ==9 || Season==11
    early=false(size(GoodData));early(1:3000)=true;
    footprint_full_temp_1=FP_process_ORWRP(ustar_1_gf(~GoodData & early),vv_1_gf(~GoodData & early),L_gf(~GoodData & early),wind_dir_1_gf(~GoodData & early),1.5,15,10,FX,FY,WetlandMap,PlotVec(~GoodData & early));
    footprint_full_temp_2=FP_process_ORWRP(ustar_1_gf(~GoodData &~early),vv_1_gf(~GoodData &~early),L_gf(~GoodData &~early),wind_dir_1_gf(~GoodData &~early),1.5,15,10,FX,FY,WetlandMap,PlotVec(~GoodData &~early));
end

footprint_full=FP_process_ORWRP(ustar_1_gf(GoodData),vv_1_gf(GoodData),L_gf(GoodData),wind_dir_1_gf(GoodData),z0,z1,d,FX,FY,WetlandMap,PlotVec(GoodData));
%PlotPoint=4815;footprint_full=FP_process_ORWRP(ustar_1_gf(PlotPoint),vv_1_gf(PlotPoint),L_gf(PlotPoint),wind_dir_1_gf(PlotPoint),z0,z1,d,FX,FY,WetlandMap,1);
if Season==4 || Season>6
    footprint_full=[footprint_full_temp_1;footprint_full;footprint_full_temp_2];
end

footprint_other_water=footprint_full(:,1+6);
footprint_roads=footprint_full(:,2+6);
footprint_lawns=footprint_full(:,3+6);
footprint_shrubs=footprint_full(:,4+6);
footprint_trees=footprint_full(:,5+6);
footprint_buildings=footprint_full(:,6+6);
footprint_wetland_upland=footprint_full(:,7+6);
footprint_wetland_water=footprint_full(:,8+6);
footprint_wetland_macro=footprint_full(:,9+6)+footprint_full(:,10+6);
footprint_wetland_edge=footprint_full(:,11+6);
footprint_wetland_tower=footprint_full(:,12+6);
InSite=sum(footprint_full(:,7:size(footprint_full,2)),2);

footprint_wetlands_combo=footprint_wetland_water+footprint_wetland_macro+footprint_wetland_edge+...
    footprint_wetland_tower+footprint_wetland_upland;

footprint_kidneys=footprint_wetland_water+footprint_wetland_edge+footprint_wetland_macro;
[~,footprint_kidneys_gf]=gapfill(footprint_kidneys,length(ustar_1)/48,48);
[~,footprint_lawns_gf]=gapfill(footprint_lawns,length(ustar_1)/48,48);


%%
%Adjust the filter to pick the right height cutoff
footprint_cutoff=0.7;
footprint_flag=~(InSite>=footprint_cutoff);
footprint_killed=sum(footprint_flag&~USTflag)/length(USTflag);
disp(['footprint_flag eliminates additional ' num2str(floor(100*footprint_killed)) '% of data'])

DayNightLineVec=nan(size(DECDAY));
TimeOfDay=(DECDAY-DOY);
DayNightLineVec(TimeOfDay<=0.5)=TimeOfDay(TimeOfDay<=0.5)*4-1;
DayNightLineVec(TimeOfDay>0.5)=TimeOfDay(TimeOfDay>0.5)*-4+3;

LENNVars = [net_sw_gf,net_lw_gf,SoilTemp_gf,ustar_1_gf,vpd_gf,DayNightLineVec];
HNNVars  = [net_sw_gf,net_lw_gf,SoilTemp_gf,ustar_1_gf,rH_gf,DayNightLineVec];

[LE,LE_gf,LE_Recon,LE_Error,H,H_gf,H_Recon,H_Error,LE_Full,LE_gf_Full,H_Full,H_gf_Full]=LEH_ORWRP_NN(NNRuns,USTflag,daynight,LE,H,LENNVars,HNNVars);

RespNNVars=[tair_gf,footprint_kidneys_gf,footprint_lawns_gf,wind_speed_1_gf,LWoutbar_gf,rH_gf,DayNightLineVec];
GPPNNVars=[tair_gf,footprint_kidneys_gf,footprint_lawns_gf,wind_speed_1_gf,LWinbar_gf,SWinbar_gf,rH_gf,DayNightLineVec];

[Resp,Resp_gf,Resp_Reconstructed,Resp_Error,...
 GPP,GPP_gf,GPP_Reconstructed,GPP_Error,...
 NEE,NEE_gf,NEE_Reconstructed,NEE_Error,...
 Resp_gf_Full,GPP_gf_Full,NEE_gf_Full]...
    =G_NN_OWRP(NNRuns,USTflag,daynight,Dspk,Fc_true,RespNNVars,GPPNNVars);

delta_P_gf=nan(size(pressure_gf));
delta_P_gf(2:end)=pressure_gf(2:end)-pressure_gf(1:end-1);delta_P_gf(isnan(delta_P_gf))=0;
if strcmp(SeasonStatus,'Summer')
    MethNNVarsDay=[DayNightLineVec,LE_gf,footprint_kidneys_gf,SoilTemp_gf,ustar_1_gf,tair_gf,NEE_gf,H_gf,net_sw_gf,pressure_gf];
    MethNNVarsNight=[DayNightLineVec,NEE_gf,SoilTemp_gf,pressure_gf,Q_gf,tair_gf,wind_speed_1_gf,H_gf];    
%     MethNNVarsDay=[DayNightLineVec,net_sw_gf,footprint_kidneys_gf,SoilTemp_gf,ustar_1_gf,wind_speed_1_gf,Q_gf,rH_gf];
%     MethNNVarsNight=[DayNightLineVec,SoilTemp_gf,footprint_kidneys_gf,pressure_gf];
else
    MethNNVarsDay=[DayNightLineVec,LE_gf,H_gf,NEE_gf,ustar_1_gf,net_sw_gf,rH_gf,footprint_kidneys_gf,SoilTemp_gf,Q_gf];
    MethNNVarsNight=[DayNightLineVec,NEE_gf,SoilTemp_gf,tair_gf,net_lw_gf,LE_gf,delta_P_gf,rH_gf];
%     MethNNVarsDay=[DayNightLineVec,net_sw_gf,wind_speed_1_gf,SoilTemp_gf,footprint_kidneys_gf,tair_gf,delta_P_gf,ustar_1_gf];
%     MethNNVarsNight=[DayNightLineVec,SoilTemp_gf,ustar_1_gf,rH_gf,wind_speed_1_gf,delta_P_gf];
end

%Neural networks for methane
%if strcmp(SeasonStatus,'Summer')
if Season~=5
[meth_gen,meth_gen_gf,meth_gen_Reconstructed,meth_gen_Error,meth_gen_gf_Full]...
   =Methane_NN_OWRP(NNRuns,[USTflag,footprint_flag],extra_vars.FCH4,daynight,Dspk,MethNNVarsNight,MethNNVarsDay);
else
   meth_gen=extra_vars.FCH4;meth_gen(USTflag | footprint_flag )=nan;
   meth_gen=DeSpike(meth_gen,Dspk.Interval,Dspk.STD,Dspk.Methane.min,Dspk.Methane.max,Dspk.trimp);
   R=NNRuns;
   meth_gen_Reconstructed_Big=nan(length(meth_gen),R);
   meth_gen_gf_Big=nan(length(meth_gen),R);
   R2_MDay=nan(1,R);
   parfor j=1:R
       [meth_gen_gf_Big(:,j), ~, meth_gen_Reconstructed_Big(:,j),R2_MDay(1,j)] = wetland_ANN_gf_Final(meth_gen,MethNNVarsDay);
       [p,~,~,~,stats]=regress(meth_gen,[meth_gen_Reconstructed_Big(:,j),ones(size(meth_gen))]);
       disp(['Methane iteration:' num2str(j) ', r2 of ' num2str(stats(1)) ', slope of ' num2str(p(1))])
   end
   BestToKeep=ceil(0.1*R);
   [~,I]=sort(R2_MDay);
   
   meth_gen_gf_Full=meth_gen_gf_Big(:,(I(end-BestToKeep+1:end)));
   meth_gen_Full=meth_gen_Reconstructed_Big(:,(I(end-BestToKeep+1:end)));

   meth_gen_gf=nanmean(meth_gen_gf_Full,2);
   meth_gen_Reconstructed=nanmean(meth_gen_Reconstructed_Big(:,(I(end-BestToKeep+1:end))),2);
   [p,~,~,~,stats]=regress(meth_gen,[meth_gen_Reconstructed,ones(size(meth_gen))]);
   disp(['Overall Methane r2: ' num2str(stats(1)) ' slope of ' num2str(p(1))]);
end

matlabpool close

save(['../VariableStorage/SeasonNum' num2str(Season) ' ' datestr(now,'yyyy-mm-dd') '.mat']);

catch err
    matlabpool close
    rethrow(err)
end
