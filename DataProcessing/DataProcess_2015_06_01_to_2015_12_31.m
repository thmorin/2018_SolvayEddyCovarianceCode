clear;close;clc
% Process data
%% Delete After finished editing 
% NOTE: Delete to everything that has to do with slow, DS, or RMY. Line 337 is where we
% start reading fast data
%%
tic
% clear all
readfromDAT=1;
%% This is the process of data from the Solvay Sediment Basin
% Created by Tim Morin on 2015_06_01. 
% Edited by Veronica Davies on 2019_10_24
%% Equipment on Tower
%CR3000 Datalogger
%Licor 7700 (CH4 sensor)
%Licor 7500 (CO2 sensor)
%CSAT3B (anemometer)
%% Data format
% File name format for 10Hz data is: thmremote_ts_data_2011_08_12_0000.dat
% File name format for 1/60Hz data is: thmremote_met_data_2011_08_12_0000.dat
%
% The signals are being processed in half hour windows. (TQUEST: would this
% be an hour files or does this mean that when we are finished, we would
% process it into 30 min intervals?)

%% Delete after finishing editing (using this portion to track names in code)
% 'Instruments on tower :
% '1.  RMyoung - sonic anemometer (model 81000) at 17 meter
% '2.  CSAT3 - sonic anemometer at 15 meter
% '3.  EC150 IRGA - combined with CSAT at 15 meter
% '4.  LI-7700 - Methane analyzer (Power Supply in the other box 24V) at 15 meter
% '5.  HMP45C - Temperature Probe at 15 meter
% '6.  NR01 - Four component radiation sensor at 15 meter
% '7.  AM1632B multiplexe with 14 107L Soil temp probs
% '8.  BF5 - Direct diffuse PAR sensor
% '9.  RF450 - Radio Transmiter in supplly box
% '10. DS2 - slow read sonic anemometer at 15 meter
%%
year=2019;
% cd('C:\Users\vldavies\Documents\GitHub\2018_SolvayEddyCovarianceCode\DataProcessing\');
cd('C:\Users\thmorin\Documents\Projects\SolvayTower\EC_code\DataProcessing\');
datadirectory =  '../SolvayRawData/';
codedirectory =  '../PreProcessSubs/';
savedir       = ['../Processed' num2str(year) '/'];
loaddir       = ['../Processed' num2str(year) '/'];

nfcol=19; %number of columns in fast file in ts data (target on one, and make a new for a different timewindow)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

addpath(genpath(codedirectory));
%DS2 installed 9 24, previous start date is 08 21

%DOY_range=Date2DOY( 06 , 26 , 2019 ):Date2DOY( 10 , 24 , 2019 ); %dates range for when we camptured data, when the col # changed
DOY_range=datenum(2019,11,05):datenum(0,0,1):datenum(2019,12,24);
fastfile = 'thmremote_ts_data_';

site = 'SLVY';                                                             % Solvay
concatenating_files = 0;                                                   % Binary. '0' not concatenating files , '1' concatenating files.

LI7500 = 1;                                                                % Infra Red Gas Analizer Type : 1='Open', 2='Closed'(AF46m), 3='Closed'(FASET32m) VERY IMPORTANT
LI7700 = 1;                                                                % Infra Red Gas Analizer Type : 1='Open', 2='Closed'(AF46m), 3='Closed'(FASET32m) VERY IMPORTANT

%% this is where we are going to make the most edits (file set up for raw data)
NhrsPfile = 1;                                                            % number of hours per file
freqF = 10;                                                                % frequency of Fast data measurements [Hz]
freqD = 1/10;                           
AVGtime = 1800 ;                                                           % Averaging time for the data input file [sec] TQUEST: So this is 30 min intervals (no change)?

ntsf = NhrsPfile*60*60*freqF;                                              % # of timesteps per file - fast data
% 24 [hr/file] * 3600 [sec/hr] *10
% [1/sec] = [number/file]

nw = floor(NhrsPfile*3600/AVGtime);                                        % number of windows per file [hrs/file] * [sec/hrs] * [1/sec] = [number/file]

NFavg = ntsf*AVGtime/(NhrsPfile*3600);                                     % Number of fast values that are meaned into one value per window
%% dont need to change
% parameters of nanWinSlaughtererT(inVec,min_nans,window,num_wins)
min_nans=600;                                                              % At least a minute of missing data
window=600;                                                                % a minut increment for calculating the STD to find max STD and throw out data for up to that point
num_wins=120;                                                              % number of increments for calculating STD
nan_STDcutoff = 3 ;                                                        % threshold for determining end of nan-edge error by using std in the problem region compared to std up to 10 windows away
%% These are the measured distances from the field (needs to be in meters)
Sdist_LI7500init = 0.23495;                                                    % Separaton distance between LI7500 and Sonic anamometer. (TQUEST: These values are twice the size of what you had previously, issue?)
Sdist_LI7700init = 0.6096;                                                    % Separaton distance between LI7700 and Sonic anamometer.
                                                     
CSAT_angle = 184;                                                          % Angle in [Deg] of CSAT from the north ClockWise (184 South) (TQUEST: so looking @ a compass, would this be 184 or 4?)
MagneticCorrection = -(6+51/60);                                           %correction due to magnetic declination for the site, degrees (leave alone)

Sonic_Height = 2.8702 ;                                                        %m
Sonic_pathLength = 0.10;                                                   %keep the same
CanpH = 1.3462;                                                                % Canopy height 
dispH = 0.67*CanpH;                                                        % Displacement Height (leave alone)
% Constants
% Gas constant dry air
R_d     =   287.0586;                                                      % Unit [Nm kg^{-1} K^{-1}] according to Foken (2005), Springer, i.e. Rg (8.314510 J mol-1 K-1 ./ 0.029002 kg mol-1), derived by Rg/ Ma with Rg= universal gas constant, Ma= molecular mass of dry air
R_v     =   461.495;                                                       % [J/(kg·K)] from wikipedia
% Molar mass of dry air
Ma = 0.029002;                                                             % [kg/mol]
% Heat capacity of the air
Cp = 1004.67;
% Universal Gas constant
R_g=8.314510;                                                              % [J/K/mol]

%%
% File name format for 10Hz data is: thmremote_ts_data_2011_03_12_0000.dat

for CD = 2:length(DOY_range) %changing 2 to 1 
%    [filemonth, fileday] = DOYtoMoDay(DOY_range(CD),year); %juilan day correction
    [fileyear,filemonth,fileday]=datevec(DOY_range(CD));
    disp(['Beginning ' num2str(year) '/' num2str(filemonth,'%02d') '/' num2str(fileday,'%02d')]);
     
     
    %% Converting Binary Files to TOA3

    tsDatOut=[];obs_time=[];
    if CD>1
        [fileyear_2,filemonth_2,fileday_2]=datevec(DOY_range(CD-1));
        
        files=dir([datadirectory 'thmremote_ts_data_' ...
            num2str(fileyear_2) '_' num2str(filemonth_2, '%02d') '_' num2str(fileday_2,'%02d') '_' num2str(23, '%02d') '*']);
        
        [tsData,tsHeaders,Var_names,raw_year,raw_month,raw_day,raw_hour,raw_minut,raw_secon]=...
             readTOB1([files(1).folder '\' files(1).name]);
         
         tsDatOut=[tsDatOut;tsData];
         obs_time=[obs_time;datenum(raw_year,raw_month,raw_day,raw_hour,raw_minut,raw_secon)];
    end
     for i=00:23
         files=dir([datadirectory 'thmremote_ts_data_' ...
            num2str(fileyear) '_' num2str(filemonth, '%02d') '_' num2str(fileday,'%02d') '_' num2str(i, '%02d') '*']);
         [tsData,tsHeaders,Var_names,raw_year,raw_month,raw_day,raw_hour,raw_minut,raw_secon]=...
             readTOB1([files(1).folder '\' files(1).name]);
         tsDatOut=[tsDatOut;tsData];
         obs_time=[obs_time;datenum(raw_year,raw_month,raw_day,raw_hour,raw_minut,raw_secon)];
     end
     
     %% Filtering daily files
     tens_of_sec = 864000; %number of tenths of a seconds in a day
    %% make perfect structure for one day
    date_skel=nan(tens_of_sec,6);
    date_skel(:,1)=fileyear;
    date_skel(:,2)=filemonth;
    date_skel(:,3)=fileday;
    date_skel(:,4)=reshape(repmat(0:23,60*60*10,1),tens_of_sec,1);
    date_skel(:,5)= reshape(repmat(0:59,60*10,24),tens_of_sec,1);
    date_skel(:,6)= reshape(repmat(0:0.1:59.9,1,24*60),tens_of_sec,1);
    date_skel_num=datenum(date_skel);
    
    [yn,whereat]=ismember(obs_time,date_skel_num);
    ts_filt=nan(tens_of_sec,size(tsDatOut,2));
    ts_filt(whereat(yn==1),:)=tsDatOut(yn==1,:);
    %      for b = 1:length(DOY_range)
%          filnamday = day(datetime(DOY_range(b),'ConvertFrom','datenum'),'dayofyear');
%          filnamday = repmat(filnamday,1,nanosec)';
%      end
     
     
     
    
      datday = day(datetime(datenum(tsDatOut(:,1)),'ConvertFrom','datenum'),'dayofyear');
 %    datday = datetime(datday,'ConvertFrom','datenum');
 %    datday = day(datday,'dayofyear');
 
%      for p = 1: length(DOY_range)
%          if (filnamday > datday)
%              greatfilt = filter(filnamday > datday);
%               tsDatOut(p+1)= [tsDatOut(p+1) greatfilt];  
%           end
%           if (filnamday < datday)
%               lessfilt = filter(filnamday < datday);
%               tsDatOut(p-1)= [tsDatOut(p-1) lessfilt];       
%           end
%       end
     
     Fdata = ts_filt;
    %%
    [DespikeW,~,fluxW,~,~,~,~,~,~,~,~,~,~,~,~,~,Header]=ORWPrealocatingVariables(nw,ntsf,NFavg,nfcol);
    currentday = DOY_range(CD);
    currentday = datetime(currentday, 'ConvertFrom','datenum');
    currentday = day(currentday,'dayofyear');
    Dspk=getDspkPrmtrs(year,currentday);
        
   
    if readfromDAT==1
        fast = dir([ datadirectory fastfile num2str(fileyear) '_' num2str(filemonth,'%02d') '_' num2str(fileday,'%02d') '*']);    % List all ts files in directory
    else
        fast = dir([loaddir 'FastData/' site '_' num2str(year) '_' num2str(currentday,'%03.0f') '*']);
    end

    
%    if ~isempty(fast)
%        [numoff,~]= size(fast);
%        if readfromDAT==1
%            [Fdata, Header.F1, Header.F2, Header.F3, Header.F4 ] = T0A2Matrix( [datadirectory fast(1).name] ); 
%            
%            [Fdata] = MakeStandard( Fdata , year , currentday , freqF ,24);
%        else
%            load([loaddir 'FastData/' site '_' num2str(year) '_' num2str(currentday,'%03.0f') '_FastData.mat']);
%       end
     
        %currentdayf = Fdata(1,2);
        currentdayf = currentday;
        %%%%%%%%%%%%%% Header of fast ts file %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %"TOA5","CR3000Olen","CR3000","3257","CR3000.Std.26","CPU:2013_08_20_ORWRP.CR3","36738","ts_data"
        %"TIMESTAMP","RECORD","UxCSAT","UyCSAT","UzCSAT","TsCSAT","diag_CSAT","CO2","H2O","LI7500Diag","LI7700Diag","LI7700_CH4_dens","LI7700_Press","RSSI","checksum_datalogger", "Timer1_m", "nBytes_LI7700"
        %   "TS"       "RN"     "m/s"    "m/s"    "m/s"    "m/s"  "unitless"   "mmol/m^3"   "unitless"   "unitless"    "mmol/m^3"           "kPa"      "unitless"   "unitless"      "unitless"  "unitless"     
        %     1          2        4        5        6        7         8        9     10       11           12            13                14         15            16                17        18   19
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        % Adjust column numbers from data including column numbers of diagnostics
        DespikeW.Time = Fdata(:,1); %changed from (:,1:3) TQUEST: what is this? 
        
        [DespikeW.CSAT_U2, DespikeW.CSAT_UQA2]= DeSpike(Fdata(:, 4),2400,Dspk.VWSTD,Dspk.VWmn,Dspk.VWmx,Dspk.VWtpr,'CSAT3',Fdata(:,8)); % CSAT Ux [m/s]
        fluxW.CSAT_UQA2 = DataQA (DespikeW.CSAT_UQA2, NFavg);
        [DespikeW.CSAT_V2, DespikeW.CSAT_VQA2] = DeSpike(Fdata(:,5),2400,Dspk.VWSTD,Dspk.VWmn,Dspk.VWmx,Dspk.VWtpr,'CSAT3',Fdata(:,8)); % CSAT Uy [m/s]
        fluxW.CSAT_VQA2 = DataQA (DespikeW.CSAT_VQA2, NFavg);
        [DespikeW.CSAT_W2, DespikeW.CSAT_WQA2] = DeSpike(Fdata(:,6),2400,Dspk.HWSTD,Dspk.HWmn,Dspk.HWmx,Dspk.HWtpr,'CSAT3',Fdata(:,8)); % CSAT Uz [m/s]
        fluxW.CSAT_WQA2 = DataQA (DespikeW.CSAT_WQA2, NFavg);
        [DespikeW.CSAT_Tmp2, DespikeW.CSAT_TmpQA2] = DeSpike(Fdata(:,7),3600,Dspk.TSTD,Dspk.Tmn,Dspk.Tmx, Dspk.Ttpr,'CSAT3',Fdata(:,8)); % CSAT Ts [Deg C]
        fluxW.CSAT_TmpQA2 = DataQA (DespikeW.CSAT_TmpQA2, NFavg);
        
        RSSIper=30; %RSSI cutoff (changed from 20)
        
        % LI7500 in mode 6 exporting CO2 & H2O [mmol/m^3], Pressure [kPa]
        % Diagnostic value (for interpetation use errorcodeLI7500.m).
        % (CRB command is: CS7500 (LI7500_CO2,1,7,6) )
	
        [DespikeW.LI7500_C2, DespikeW.LI7500_CQA2] = DeSpike (Fdata(:,9)/44,3600,Dspk.CSTD,Dspk.Cmn,Dspk.Cmx,Dspk.Ctpr,'LI7500',Fdata(:,12));      % LI7500 CO2 [mmol/m^3] { Fdata(:,17) = LI7500 Diag}
        DespikeW.LI7500_C2=nanWinSlaughterer(DespikeW.LI7500_C2,min_nans,window,num_wins,nan_STDcutoff);
        fluxW.LI7500_CQA2 = DataQA (DespikeW.LI7500_CQA2, NFavg);
        [DespikeW.LI7500_Q2, DespikeW.LI7500_QQA2] = DeSpike (Fdata(:,10)*1000/18,3600,Dspk.QSTD,Dspk.Qmn,Dspk.Qmx,Dspk.Qtpr,'LI7500',Fdata(:,12));    % LI7500 H2O [mmol/m^3] TQUEST: we do not have h20signal or CO2 signal
        DespikeW.LI7500_Q2=nanWinSlaughterer(DespikeW.LI7500_Q2,min_nans,window,num_wins,nan_STDcutoff);
        fluxW.LI7500_QQA2 = DataQA (DespikeW.LI7500_QQA2, NFavg);
        [DespikeW.LI7500_P2, DespikeW.LI7500_PQA2] = DeSpike(Fdata(:,11),2400,Dspk.PSTD,Dspk.Pmn,Dspk.Pmx,Dspk.Ptpr,'LI7500',Fdata(:,12));     % LI7500 Pressure [ Pa ]
        fluxW.LI7500_PQA2 = DataQA (DespikeW.LI7500_PQA2, NFavg);
        DespikeW.LI7500_P2=1000*DespikeW.LI7500_P2;

        
        % LI7700 Signal Strength [ max 70 min 0, the higher the better ]
        
        %Z_old = -5.02946; Z_new = -14.2544;
        %S_old = 1.97176e-06; S_new = 0.000106478;
        [DespikeW.LI7700_M, DespikeW.LI7700_MQA] = DeSpike (Fdata(:,13),6000,Dspk.MSTD,Dspk.Mmn,Dspk.Mmx,Dspk.Mtpr,'LI7700',Fdata(:,15),Fdata(:,12),RSSIper);      % LI7700 CH4 [mmol/m^3] 
        
        DespikeW.LI7700_M=nanWinSlaughterer(DespikeW.LI7700_M,min_nans,window,num_wins,nan_STDcutoff);
        fluxW.LI7700_MQA = DataQA (DespikeW.LI7700_MQA, NFavg);
        [DespikeW.LI7700_P, DespikeW.LI7700_PQA] = DeSpike(Fdata(:,14),2400,Dspk.PSTD,Dspk.Pmn,Dspk.Pmx,Dspk.Ptpr,'LI7700',Fdata(:,12));     % LI7700 Pressure [Pa] 
        fluxW.LI7700_PQA = DataQA (DespikeW.LI7700_PQA, NFavg);
        DespikeW.LI7700_P=1000*DespikeW.LI7700_P;
        %Forces methane to nan when P is nan. P will normally only be nan if mirror is spinning
        DespikeW.LI7700_M(isnan(DespikeW.LI7700_P))=nan;
%    end %~isempty(fast)
    
        save([savedir 'DespikeW/' site '_' num2str(year) '_' num2str(currentday,'%03.0f') '_2_DespikeW.mat'],'DespikeW');
        save([savedir 'FastData/' site '_' num2str(year) '_' num2str(currentday,'%03.0f') '_2_FastData.mat'],'Fdata');
        save([savedir 'DataHeader/' site '_' num2str(year) '_' num2str(currentday,'%03.0f') '_2_DataHeader.mat'],'Header');
    
    clear Fdata Header tmpTEMPRMY
    
    for CW = 1:nw   % CW = Current Window ; nw = Number Of Windows per file
        
        fluxW=ORW_fluxW_process(fluxW,CW,NFavg,ntsf,DespikeW,...
            MagneticCorrection, CSAT_angle,...
            LI7500,Sdist_LI7500init,LI7700,Sdist_LI7700init,...
            Sonic_pathLength,Sonic_Height,AVGtime,dispH, freqF,nw, concatenating_files,...
            R_d);
        
    end %for CW (CW = Current Window)
    
    disp([savedir site '_' num2str(year) '_' num2str(currentday,'%03.0f') '_flux.mat; of raw file: ' fast.name])
    
    save([savedir 'flux/' site '_' num2str(year) '_' num2str(currentday,'%03.0f') '_flux.mat'],'fluxW');
    toc
end     % for the files in Data folder
