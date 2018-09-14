% Process data
tic
% clear all
readfromDAT=1;

%%%%% This is the process of data from the Olentangy River Wetland
%%%%% Research Park
% 2013 10 10: Creating this process to run data after the switch to the
% EC150 sonic anemometer. Includes DS2 sonic anemometer at 15 m also
%
% Data Logger is CR3000
%
% File name format for 10Hz data is: ts_data_2011_08_12_0000.dat
% File name format for 1/60Hz data is: met_data_2011_08_12_0000.dat
%
% The signals are being processed in half hour windows.
%
% Header for ts file
%
% Header for met file
%
%
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
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
year=2015;

% cd('/home/morin.37/poolA/wetland_data/Matlab_Tim/ORWRP/DataProcessing/');
datadirectory =  ['../ORWRawData/' num2str(year) '/'];
codedirectory =  '../../PreProcessSubs/';
savedir       = ['../Processed' num2str(year) '/'];
loaddir       = ['../Processed' num2str(year) '/'];

nscol=48; %number of columns in slow file
nfcol=33; %number of columns in fast file

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

addpath(genpath(codedirectory));
%DS2 installed 9 24, previous start date is 08 21
DOY_range=Date2DOY( 07 , 13 , 2015 ):Date2DOY( 07 , 26 , 2015 );
% Special treatment is needed for 05_06_2013

fastfile = 'ts_data_';
slowfile = 'met_data_';
dsfile   = 'ds_data_';

site = 'OWR';                                                              % Olentangy Wetland Research
concatenating_files = 0;                                                   % Binary. '0' not concatenating files , '1' concatenating files.

LI7500 = 1;                                                                % Infra Red Gas Analizer Type : 1='Open', 2='Closed'(AF46m), 3='Closed'(FASET32m) VERY IMPORTANT
LI7700 = 1;                                                                % Infra Red Gas Analizer Type : 1='Open', 2='Closed'(AF46m), 3='Closed'(FASET32m) VERY IMPORTANT

NhrsPfile = 24;                                                            % number of hours per file
freqF = 10;                                                                % frequency of Fast data measurements [Hz]
freqS = 1/60;                                                              % frequency of Slow data measurements [Hz] (one measurment per minute)
freqD = 1/10;                           
AVGtime = 1800 ;                                                           % Aveging time for the data input file [sec]

ntsf = NhrsPfile*60*60*freqF;                                              % # of timesteps per file - fast data
% 24 [hr/file] * 3600 [sec/hr] *10
% [1/sec] = [number/file]
ntss = NhrsPfile*60*60*freqS;                                              % # of timesteps per file - slow data

ntsd = NhrsPfile*60*60*freqD;                                              % # of timesteps per file - ds data

nw = floor(NhrsPfile*3600/AVGtime);                                        % number of windows per file [hrs/file] * [sec/hrs] * [1/sec] = [number/file]

NFavg = ntsf*AVGtime/(NhrsPfile*3600);                                     % Number of fast values that are meaned into one value per window
NSavg = ntss*AVGtime/(NhrsPfile*3600);                                     % Number of slow values that are meaned into one value per window

% parameters of nanWinSlaughtererT(inVec,min_nans,window,num_wins)
min_nans=600;                                                              % At least a minute of missing data
window=600;                                                                % a minut increment for calculating the STD to find max STD and throw out data for up to that point
num_wins=120;                                                              % number of increments for calculating STD
nan_STDcutoff = 3 ;                                                        % threshold for determining end of nan-edge error by using std in the problem region compared to std up to 10 windows away

Sdist_LI7500init = 0.1;                                                    % Separaton distance between LI7500 and Sonic anamometer.
Sdist_LI7700init = 0.3;                                                    % Separaton distance between LI7700 and Sonic anamometer.
%*********************
RMY_angle = 212-180;                                                       % Angle in [Deg] of RMY from the north CW
CSAT_angle = 359;                                                          % Angle in [Deg] of CSAT from the north CW
MagneticCorrection = -(6+51/60);                                           %correction due to magnetic declination for the site, degrees

Sonic_Height = 12 ;
Sonic_pathLength = 0.10;
CanpH = 10;                                                                % Canopy height
dispH = 0.67*CanpH;                                                        % Displacement Height
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

% File name format for 10Hz data is: ts_data_2011_03_12_0000.dat
% File name format for 1/60Hz data is: met_data_2011_03_12_0000.dat


for CD = 2:length(DOY_range)
    
    [filemonth, fileday] = DOYtoMoDay(DOY_range(CD),year);
    
    disp(['Beginning ' num2str(year) '/' num2str(filemonth,'%02d') '/' num2str(fileday,'%02d')]);

    [W1min,DespikeW,~,fluxW,~,~,~,~,~,~,~,~,~,~,~,~,Sdata,Fdata,Header]=ORWPrealocatingVariables(ntss,nw,ntsf,NFavg,nscol,nfcol);
    currentday = DOY_range(CD);
    
    DsData = nan(1440,6);    
    Dspk=getDspkPrmtrs(year,currentday);
        
    if readfromDAT==1
        slow = dir([datadirectory slowfile num2str(year) '_' num2str(filemonth,'%02.0f') '_' num2str(fileday,'%02.0f') '*']);    % List all met files in directory
    else
        slow = dir([loaddir 'SlowData/' site '_' num2str(year) '_' num2str(currentday,'%03.0f') '*']);
    end

    W1minWin=10;
    
    if ~isempty(slow)
        
        [numosf,~]= size(slow);
        if readfromDAT==1
            [Sdata, Header.S1, Header.S2, Header.S3, Header.S4 ] = T0A2Matrix( [datadirectory slow(1).name] );
            [Sdata] = MakeStandard( Sdata , year , currentday , freqS ,24);
        else
            load([loaddir 'SlowData/' site '_' num2str(year) '_' num2str(currentday,'%03.0f') '_SlowData.mat']);
        end
        
        
        
        W1min.year = Sdata(:,1);
        W1min.doy = Sdata(:,2);
        W1min.dectime = Sdata(:,3);
        
        %%%%%%%%%%%%%% Header of slow met file %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %"TOA5","CR3000Olen","CR3000","3257","CR3000.Std.09.02","CPU:2012_11_01_2_ORWRP.CR3","25871","met_data"
        %"TIMESTAMP","RECORD","t_hmp_Avg","rh_hmp_Avg","e_hmp_Avg","batt_volt_Avg","panel_temp_Avg","SR01Up_Avg","SR01Dn_Avg","IR01Up_Avg","IR01Dn_Avg","TC_Avg","Total_PAR_Avg","Diffuse_PAR_Avg","CO2at2m_Avg","CO2at5m_Avg","CO2at10m_Avg","STat8cmOW_W1_Avg","STat25cmOW_W1_Avg","STat8cmOWr_W1_Avg","STat8cmIM_W1_Avg","STat25cmIM_W1_Avg","STat8cmIMr_W1_Avg","STat8cmUL_W1_Avg","STat8cmOW_W2_Avg","STat25cmOW_W2_Avg","STat8cmOWr_W2_Avg","STat8cmIM_W2_Avg","STat25cmIM_W2_Avg","STat8cmIMr_W2_Avg","STat8cmUL_W2_Avg","PropAnemWS_Avg","PropAnemWD_Avg","DS2Gust_Max"
        %  "TS"       "RN"       "Deg C"       "%"        "kPa"           "V"         "Deg C"           "W/m2"      "W/m2"       "W/m2"       "W/m2"     "Deg C"      "W/m2"         "W/m2"           "ppm"           "ppm"         "ppm"         "Deg C"             "Deg C"              "Deg C"           "Deg C"              "Deg C"             "Deg C"            "Deg C"             "Deg C"            "Deg C"           "Deg C"             "Deg C"              "Deg C"             "Deg C"              "Deg C"           "m/s"            "deg"          "m/s"
        %%% 1           2           4           5           6              7             8                9           10          11             12          13         14             15               16              17           18             19                  20                   21                 22                   23                  24                25                  26                 27                28                  29                    30                  31                  32                33               34            35
        %%%%%%%%%%%%%% Header of fast ts file %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %"TOA5","CR3000Olen","CR3000","3257","CR3000.Std.26","CPU:2013_08_20_ORWRP.CR3","36738","ts_data"
        %"TIMESTAMP","RECORD","UxRMY","UyRMY","UzRMY","sosRMY","diag_RMY","UxCSAT","UyCSAT","UzCSAT","TsCSAT","diag_CSAT","CO2","H2O","diag_irga","Ta_irga","press_irga","co2_signal","h2o_signal","LI7700Diag","LI7700_CH4_dens","LI7700_Press","RSSI"
        %   "TS"       "RN"    "m/s"   "m/s"   "m/s"    "m/s"  "unitless"  "m/s"     "m/s"    "m/s"   "DegC"  "unitless"  "mg/m3"g/m3" "unitless"   "DegC"      "kPa"     "fraction"       ""           ""         "mmol/m^3"         "kPa"        "%"
        %     1          2       4       5       6        7         8        9        10        11      12        13        14    15        16        17         18           19           20           21            22                23          24
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
           % 1                    4            5           6           7                 8               9            10           11           12         13          14               15                  16                17                     18              19                  20                  21                   22               23                  24                   25                 26                 27                    28                29                  30                 31                 32                    33                   34                     35               36        37          38              39              40                41               42             43               44                 45                     46                47          48        
        %"TIMESTAMP","RECORD","t_hmp_Avg","rh_hmp_Avg","e_hmp_Avg","batt_volt_Avg","panel_temp_Avg","SR01Up_Avg","SR01Dn_Avg","IR01Up_Avg","IR01Dn_Avg","TC_Avg","Total_PAR_Avg","Diffuse_PAR_Avg","STat8cmOW_W1_Avg","STat25cmOW_W1_Avg","STat8cmOWr_W1_Avg","STat8cmIM_W1_Avg","STat25cmIM_W1_Avg","STat8cmIMr_W1_Avg","STat8cmUL_W1_Avg","STat8cmOW_W2_Avg","STat25cmOW_W2_Avg","STat8cmOWr_W2_Avg","STat8cmIM_W2_Avg","STat25cmIM_W2_Avg","STat8cmIMr_W2_Avg","STat8cmUL_W2_Avg","Wind_Dir_min_Min","Wind_Dir_mean_Avg","Wind_Dir_max_Max","Wind_Speed_min_Min","Wind_Speed_mean_Avg","Wind_Speed_max_Max","Air_temp_Avg","RH_Avg","Pressure_Avg","Rain_accum_Max","Rain_dur_Max","Rain_intens_Avg","Hail_accum_Max","Hail_dur_Max","Hail_intens_Avg","Rain_pk_intens_Max","Hail_pk_intens_Max","GMP_2_m_Avg","GMP_5_m_Avg"

        [W1min.Tm, W1min.TmQA] = DeSpike( Sdata(:,28) ,W1minWin,Dspk.TSTD,Dspk.Tmn,Dspk.Tmx, Dspk.Ttpr);   % HMP45C Temp [ Deg C ]
        fluxW.tairQA = DataQA(W1min.TmQA,NSavg);
        Sdata((Sdata(:,29)>100 & Sdata(:,29)<110),29)=100;
        [W1min.Hum, W1min.HumQA] = DeSpike(Sdata(:,29),W1minWin,Dspk.HSTD,Dspk.Hmn,Dspk.Hmx,Dspk.Htpr );    % HMP45C Relative humidity [ % ]
        fluxW.humbarQA= DataQA(W1min.HumQA,NSavg);
        %[W1min.Pvapor, W1min.PvaporQA] = DeSpike(Sdata(:,6),W1minWin,Dspk.VPSTD,Dspk.VPmn,Dspk.VPmx,Dspk.VPtpr );   % HMP45C Vapor pressure [mPa] = [N/m^2] = [Kg/(m*s^2)]=[J/m^3]
        %fluxW.p_vapor_barQA= DataQA(W1min.PvaporQA,NSavg);
        W1min.Pvapor = W1min.Pvapor*1000;   % [Pa]
        W1min.Q = W1min.Pvapor/R_g./(W1min.Tm+273.16);    % HMP45C H2O Density [mol/m^3]  {Ideal gas law, n/V=Q=P/(R_gT)}
        W1min.Pvaporsat = 611.2*exp(17.67*W1min.Tm./(243.51+W1min.Tm)); % HMP45C Saturated Vapor Pressure [Pa]
        
        [W1min.Tcnr01, W1min.Tcnr01QA] = DeSpike(Sdata(:,13-3),W1minWin,Dspk.TSTD,Dspk.Tmn,Dspk.Tmx, Dspk.Ttpr );  % NR01 TC Temp [ Deg C ]
        % Calculations for Radiation data
        [W1min.NetRs, W1min.NetRl, W1min.Albedo, W1min.UpTot, W1min.DownTot, W1min.NetRad, W1min.IRUp, W1min.IRDown, W1min.SWUp, W1min.SWDown, QAstruct] = ...
            calculate_radiation( Sdata(:,9-3), Sdata(:,10-3), Sdata(:,11-3), Sdata(:,12-3) , W1min.Tcnr01, Dspk);
        W1min.SWUpQA = QAstruct.SR01UpQA ; fluxW.SWinbarQA = DataQA(W1min.SWUpQA,NSavg);
        W1min.SWDownQA = QAstruct.SR01DnQA ; fluxW.SWoutbarQA = DataQA(W1min.SWDownQA,NSavg);
        W1min.IRUpQA = QAstruct.IRUpQA ; fluxW.LWinbarQA = DataQA(W1min.IRUpQA,NSavg);
        W1min.IRDownQA = QAstruct.IRDownQA ; fluxW.LWoutbarQA = DataQA(W1min.IRDownQA,NSavg);
        
        
        [W1min.Total_PAR, W1min.Total_PARQA] = DeSpike( Sdata(:,14-3),W1minWin,Dspk.TPSTD,Dspk.TPmn,Dspk.TPmx,Dspk.TPtpr);  %BF5 Total PAR [W/m^2](1 mV = 0.5 W/m^2)
        fluxW.Total_PARQA = DataQA(W1min.Total_PARQA,NSavg);
        W1min.Total_PAR = 0.5*W1min.Total_PAR;
        [W1min.Diffuse_PAR, W1min.Diffuse_PARQA] = DeSpike( Sdata(:,15-3),W1minWin,Dspk.DPSTD,Dspk.DPmn,Dspk.DPmx,Dspk.DPtpr);  %BF5 Diffuse PAR [W/m^2](1 mV = .5 W/m^2)
        fluxW.Diffuse_PARQA = DataQA(W1min.Diffuse_PARQA,NSavg);
        W1min.Diffuse_PAR = 0.5* W1min.Diffuse_PAR;
    
        %[W1min.C1, W1min.C1QA] = DeSpike( Sdata(:,16),W1minWin,Dspk.CO2STD,Dspk.CO2mn,Dspk.CO2mx,Dspk.CO2tpr);  %ID200 CO2 concentration [ppm]
        %[W1min.C2, W1min.C2QA] = DeSpike( Sdata(:,17),W1minWin,Dspk.CO2STD,Dspk.CO2mn,Dspk.CO2mx,Dspk.CO2tpr);  %ID200 CO2 concentration [ppm]
        %[W1min.C3, W1min.C2QA] = DeSpike( Sdata(:,18),W1minWin,Dspk.CO2STD,Dspk.CO2mn,Dspk.CO2mx,Dspk.CO2tpr);  %ID200 CO2 concentration [ppm]
%         fluxW.C1QA = DataQA (W1min.C1QA, NSavg);
%         fluxW.C2QA = DataQA (W1min.C2QA, NSavg);
%         fluxW.C3QA = DataQA (W1min.C3QA, NSavg); 

        [W1min.STat8cmOW_W1, W1min.STat8cmOW_W1QA] = DeSpike(Sdata(:,16-3),W1minWin,Dspk.STSTD,Dspk.STmn,Dspk.STmx,Dspk.STtpr);  % 107L soil temperature: 8cm into soil, in wetland 1, under open water
        fluxW.STat8cmOW_W1QA = DataQA (W1min.STat8cmOW_W1QA, NSavg);
        [W1min.STat25cmOW_W1, W1min.STat25cmOW_W1QA] = DeSpike(Sdata(:,17-3),W1minWin,Dspk.STSTD,Dspk.STmn,Dspk.STmx,Dspk.STtpr);  % 107L soil temperature: 25cm into soil, in wetland 1, under open water
        fluxW.STat25cmOW_W1QA = DataQA (W1min.STat25cmOW_W1QA, NSavg);
        [W1min.STat8cmOWr_W1, W1min.STat8cmOWr_W1QA] = DeSpike(Sdata(:,18-3),W1minWin,Dspk.STSTD,Dspk.STmn,Dspk.STmx,Dspk.STtpr);  % 107L soil temperature: 8cm into soil, in wetland 1, under open water replica
        fluxW.STat8cmOWr_W1QA = DataQA (W1min.STat8cmOWr_W1QA, NSavg);
        [W1min.STat8cmIM_W1, W1min.STat8cmIM_W1QA] = DeSpike(Sdata(:,19-3),W1minWin,Dspk.STSTD,Dspk.STmn,Dspk.STmx,Dspk.STtpr);    % 107L soil temperature: 8cm into soil, in wetland 1, intermediate
        fluxW.STat8cmIM_W1QA = DataQA (W1min.STat8cmIM_W1QA, NSavg);
        [W1min.STat25cmIM_W1, W1min.STat25cmIM_W1QA] = DeSpike(Sdata(:,20-3),W1minWin,Dspk.STSTD,Dspk.STmn,Dspk.STmx,Dspk.STtpr);   % 107L soil temperature: 25cm into soil, in wetland 1, intermediate
        fluxW.STat25cmIM_W1QA = DataQA (W1min.STat25cmIM_W1QA, NSavg);
        [W1min.STat8cmIMr_W1, W1min.STat8cmIMr_W1QA] = DeSpike(Sdata(:,21-3),W1minWin,Dspk.STSTD,Dspk.STmn,Dspk.STmx,Dspk.STtpr);   % 107L soil temperature: 8cm into soil, in wetland 1, intermediate replica
        fluxW.STat8cmIMr_W1QA = DataQA (W1min.STat8cmIMr_W1QA, NSavg);
        [W1min.STat8cmUL_W1, W1min.STat8cmUL_W1QA] = DeSpike(Sdata(:,22-3),W1minWin,Dspk.STSTD,Dspk.STmn,Dspk.STmx,Dspk.STtpr);     % 107L soil temperature: 8cm into soil, in wetland 1, upland
        fluxW.STat8cmUL_W1QA = DataQA (W1min.STat8cmUL_W1QA, NSavg);
        
        [W1min.WD_min, W1min.WD_min_QA] = DeSpike(Sdata(:,30-10),W1minWin,6,0,360,1);
        fluxW.WD_min_QA = DataQA (W1min.WD_min_QA, NSavg);
        [W1min.WD_bar, W1min.WD_bar_QA] = DeSpike(Sdata(:,31-10),W1minWin,6,0,360,1);
        fluxW.WD_bar_QA = DataQA (W1min.WD_bar_QA, NSavg);
        [W1min.WD_max, W1min.WD_max_QA] = DeSpike(Sdata(:,35-10),W1minWin,6,0,360,1);
        fluxW.WD_max_QA = DataQA (W1min.WD_max_QA, NSavg);
        
        [W1min.WS_min, W1min.WS_min_QA] = DeSpike(Sdata(:,34-10),W1minWin,Dspk.VWSTD,Dspk.VPmn,Dspk.VWmx,Dspk.VWtpr);
        fluxW.WS_min_QA = DataQA (W1min.WS_min_QA, NSavg);
        
        [W1min.WS_bar, W1min.WS_bar_QA] = DeSpike(Sdata(:,33-10),W1minWin,Dspk.VWSTD,Dspk.VPmn,Dspk.VWmx,Dspk.VWtpr);
        fluxW.WS_bar_QA = DataQA (W1min.WS_bar_QA, NSavg);
        
        [W1min.WS_max, W1min.WS_max_QA] = DeSpike(Sdata(:,35-10),W1minWin,Dspk.VWSTD,Dspk.VPmn,Dspk.VWmx,Dspk.VWtpr);
        fluxW.WS_max_QA = DataQA (W1min.WS_max_QA, NSavg);
        
        [W1min.tair_2, W1min.tair_2_QA] = DeSpike(Sdata(:,38-10),W1minWin,Dspk.TSTD,Dspk.Tmn,Dspk.Tmx,Dspk.Ttpr);
        fluxW.tair_2_QA = DataQA (W1min.tair_2_QA, NSavg);
        
        [W1min.rH_2, W1min.rH_2_QA] = DeSpike(Sdata(:,39-10),W1minWin,6,0,100,1);
        fluxW.rH_2_QA = DataQA (W1min.rH_2_QA, NSavg);
        
        [W1min.pres,W1min.pres_QA] = DeSpike(Sdata(:,40-10),W1minWin,6,900,1200,1);
        fluxW.pres_QA = DataQA (W1min.pres_QA, NSavg);
        
       
        W1min.rain_accum    = Sdata(:,41-10);
        W1min.rain_dur      = Sdata(:,42-10);
        W1min.rain_intens   = Sdata(:,43-10);
        W1min.hail_accum    = Sdata(:,44-10);
        W1min.hail_dur      = Sdata(:,45-10);
        W1min.hail_intens   = Sdata(:,46-10);

        [W1min.C1, W1min.C1_QA] = DeSpike(Sdata(:,49-10),W1minWin,Dspk.CO2STD,Dspk.CO2mn,Dspk.CO2mx,Dspk.CO2tpr);
        fluxW.C1_QA = DataQA (W1min.C1_QA, NSavg);
        [W1min.C2, W1min.C2_QA] = DeSpike(Sdata(:,50-10),W1minWin,Dspk.CO2STD,Dspk.CO2mn,Dspk.CO2mx,Dspk.CO2tpr);
        fluxW.C2_QA = DataQA (W1min.C2_QA, NSavg);

     end %~isempty(slow)
    
      
    if readfromDAT==1
        ds = dir([datadirectory dsfile num2str(year) '_' num2str(filemonth,'%02.0f') '_' num2str(fileday,'%02.0f') '*']);    % List all met files in directory
    else
        ds = dir([loaddir 'DsData/' site '_' num2str(year) '_' num2str(currentday,'%03.0f') '*']);
    end
    
    if ~isempty(ds)
        
        [numodf,~]= size(ds);
        if readfromDAT==1
            [DsData, Header.DS1, Header.DS2, Header.DS3, Header.DS4 ] = T0A2Matrix( [datadirectory ds(1).name] );
            [DsData] = MakeStandard( DsData , year , currentday , freqD ,24);
        else
            load([loaddir 'DsData/' site '_' num2str(year) '_' num2str(currentday,'%03.0f') '_DsData.mat']);
        end
     if freqD~=freqS
         S = DsData(:,4);
         D = DsData(:,5);
         X = -S.*sin(D*pi/180);
         Y = -S.*cos(D*pi/180);
         
         avInt = freqD/freqS;
         Xmod = nanmean(reshape(X,avInt,length(X)/avInt))';
         Ymod = nanmean(reshape(Y,avInt,length(X)/avInt))';
         
         YY   = nanmean(reshape(DsData(:,1),avInt,length(X)/avInt))';
         DD   = nanmean(reshape(DsData(:,2),avInt,length(X)/avInt))';
         HH   = nanmean(reshape(DsData(:,3),avInt,length(X)/avInt))';
         TT   = nanmean(reshape(DsData(:,6),avInt,length(X)/avInt))';
         
         DsData = nan(1440,size(DsData,2));
         DsData(:,1) = YY;
         DsData(:,2) = DD;
         DsData(:,3) = HH;
         DsData(:,4) = sqrt(Xmod.^2 + Ymod.^2);
         DsData(:,5) = mod(180/pi*atan2(Ymod,Xmod)+90,360);
         DsData(:,6) = TT;
     end
     
    [W1min.DS2_ubar, W1min.DS2_uQA] = DeSpike (DsData(:,4),W1minWin,Dspk.VWSTD,Dspk.VWmn,Dspk.VWmx,Dspk.VWtpr);               % DS2 sonic anemometer wind speed [m/s]
     W1min.DS2_dir = DsData(:,5);                                                                                           % DS2 direction off north, clockwise [deg]
            
     W1min.DS2_u = Xmod;
     W1min.DS2_v = Ymod;
%     [W1min.DS2_gust, W1min.DS2_gustQA] = DeSpike (DsData(:,35),W1minWin,Dspk.VWSTD,Dspk.VWmn,50,Dspk.VWtpr);                % DS2 sonic anemometer gust speed [m/s]
    end
    
    if CD~=1
        save([savedir 'W1min/' site '_' num2str(year) '_' num2str(currentday,'%03.0f') '_' 'W1min.mat'],'W1min');
    else
        save([savedir 'W1min/' site '_' num2str(year) '_' num2str(currentday,'%03.0f') '_2_' 'W1min.mat'],'W1min');
    end
    
    if readfromDAT==1
        if CD~=1
            save([savedir 'SlowData/' site '_' num2str(year) '_' num2str(currentday,'%03.0f') '_SlowData.mat'],'Sdata');
	    save([savedir 'dsData/'   site '_' num2str(year) '_' num2str(currentday,'%03.0f') '_dsData.mat'  ],'DsData');
        else
            save([savedir 'SlowData/' site '_' num2str(year) '_' num2str(currentday,'%03.0f') '_2_SlowData.mat'],'Sdata');
	    save([savedir 'dsData/'   site '_' num2str(year) '_' num2str(currentday,'%03.0f') '_2_dsData.mat'  ],'DsData');
        end
    end
    clear Sdata
    
    
    if readfromDAT==1
        fast = dir([ datadirectory fastfile num2str(year) '_' num2str(filemonth,'%02.0f') '_' num2str(fileday,'%02.0f') '*']);    % List all ts files in directory
    else
        fast = dir([loaddir 'FastData/' site '_' num2str(year) '_' num2str(currentday,'%03.0f') '*']);
    end
    
    if ~isempty(fast)
        [numoff,~]= size(fast);
        if readfromDAT==1
            [Fdata, Header.F1, Header.F2, Header.F3, Header.F4 ] = T0A2Matrix( [datadirectory fast(1).name] );
            
            [Fdata] = MakeStandard( Fdata , year , currentday , freqF ,24);
        else
            load([loaddir 'FastData/' site '_' num2str(year) '_' num2str(currentday,'%03.0f') '_FastData.mat']);
        end
     
        currentdayf = Fdata(1,2);
        
        %%%%%%%%%%%%%% Header of fast ts file %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %"TOA5","CR3000Olen","CR3000","3257","CR3000.Std.26","CPU:2013_08_20_ORWRP.CR3","36738","ts_data"
        %"TIMESTAMP","RECORD","UxRMY","UyRMY","UzRMY","sosRMY","diag_RMY","UxCSAT","UyCSAT","UzCSAT","TsCSAT","diag_CSAT","CO2","H2O","diag_irga","Ta_irga","press_irga","co2_signal","h2o_signal","LI7700Diag","LI7700_CH4_dens","LI7700_Press","RSSI"
        %   "TS"       "RN"    "m/s"   "m/s"   "m/s"    "m/s"  "unitless"  "m/s"     "m/s"    "m/s"   "DegC"  "unitless"  "mg/m3"g/m3" "unitless"   "DegC"      "kPa"     "fraction"       ""           ""         "mmol/m^3"         "kPa"        "%"
        %     1          2       4       5       6        7         8        9        10        11      12        13        14    15        16        17         18           19           20           21            22                23          24
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        % Adjust column numbers from data including column numbers of diagnostics
        DespikeW.Time = Fdata(:,1:3);
        
        [DespikeW.RMY_U, DespikeW.RMY_UQA] = DeSpike(Fdata(:,4),1200,Dspk.VWSTD,Dspk.VWmn,Dspk.VWmx,Dspk.VWtpr,'RMYOUNG',Fdata(:,8)); % RMYoung Ux [m/s]
        fluxW.RMY_UQA = DataQA (DespikeW.RMY_UQA, NFavg);
        [DespikeW.RMY_V, DespikeW.RMY_VQA]= DeSpike(Fdata(:,5),1200,Dspk.VWSTD,Dspk.VWmn,Dspk.VWmx,Dspk.VWtpr,'RMYOUNG',Fdata(:,8)); % RMYoung Uy [m/s]
        fluxW.RMY_VQA = DataQA (DespikeW.RMY_VQA, NFavg);
        [DespikeW.RMY_W, DespikeW.RMY_WQA]= DeSpike(Fdata(:,6),1200,Dspk.HWSTD,Dspk.HWmn,Dspk.HWmx,Dspk.HWtpr,'RMYOUNG',Fdata(:,8)); % RMYoung Uz [m/s]
        fluxW.RMY_WQA = DataQA (DespikeW.RMY_WQA, NFavg);
        tmpTEMPRMY = Fdata(:,7);                 % RMYoung Speed Of Sound - SOS
        tmpTEMPRMY = convertSOS_Tc( tmpTEMPRMY); % RMYoung convert SOS to Temp T [deg C]
        tmpTEMPRMY = Schotanus( tmpTEMPRMY,DespikeW.RMY_U,DespikeW.RMY_V,DespikeW.RMY_W); % RMYoung Corrected T (correcting error that is caused due to sonic structure) [deg C]
        [DespikeW.RMY_Tmp, DespikeW.RMY_TmpQA] = DeSpike(tmpTEMPRMY,3600,Dspk.TSTD,Dspk.Tmn,Dspk.Tmx, Dspk.Ttpr,'RMYOUNG',Fdata(:,8));% RMYoung temperature without spikes [deg C]
        fluxW.RMY_TmpQA = DataQA (DespikeW.RMY_TmpQA, NFavg);

        [DespikeW.CSAT_U2, DespikeW.CSAT_UQA2]= DeSpike(Fdata(:, 9),2400,Dspk.VWSTD,Dspk.VWmn,Dspk.VWmx,Dspk.VWtpr,'EC150',Fdata(:,13)); % CSAT Ux [m/s]
        fluxW.CSAT_UQA2 = DataQA (DespikeW.CSAT_UQA2, NFavg);
        [DespikeW.CSAT_V2, DespikeW.CSAT_VQA2] = DeSpike(Fdata(:,10),2400,Dspk.VWSTD,Dspk.VWmn,Dspk.VWmx,Dspk.VWtpr,'EC150',Fdata(:,13)); % CSAT Uy [m/s]
        fluxW.CSAT_VQA2 = DataQA (DespikeW.CSAT_VQA2, NFavg);
        [DespikeW.CSAT_W2, DespikeW.CSAT_WQA2] = DeSpike(Fdata(:,11),2400,Dspk.HWSTD,Dspk.HWmn,Dspk.HWmx,Dspk.HWtpr,'EC150',Fdata(:,13)); % CSAT Uz [m/s]
        fluxW.CSAT_WQA2 = DataQA (DespikeW.CSAT_WQA2, NFavg);
        [DespikeW.CSAT_Tmp2, DespikeW.CSAT_TmpQA2] = DeSpike(Fdata(:,12),3600,Dspk.TSTD,Dspk.Tmn,Dspk.Tmx, Dspk.Ttpr,'EC150',Fdata(:,13)); % CSAT Ts [Deg C]
        fluxW.CSAT_TmpQA2 = DataQA (DespikeW.CSAT_TmpQA2, NFavg);
        
        RSSIper=20;
        % LI7500 in mode 6 exporting CO2 & H2O [mmol/m^3], Pressure [kPa]
        % Diagnostic value (for interpetation use errorcodeLI7500.m).
        % (CRB command is: CS7500 (LI7500_CO2,1,7,6) )
	      
       
    end %~isempty(fast)
    
    if CD~=1
        save([savedir 'DespikeW/' site '_' num2str(year) '_' num2str(currentday,'%03.0f') '_DespikeW.mat'],'DespikeW');
        save([savedir 'FastData/' site '_' num2str(year) '_' num2str(currentday,'%03.0f') '_FastData.mat'],'Fdata');
        save([savedir 'DataHeader/' site '_' num2str(year) '_' num2str(currentday,'%03.0f') '_DataHeader.mat'],'Header');
    else
        save([savedir 'DespikeW/' site '_' num2str(year) '_' num2str(currentday,'%03.0f') '_2_DespikeW.mat'],'DespikeW');
        save([savedir 'FastData/' site '_' num2str(year) '_' num2str(currentday,'%03.0f') '_2_FastData.mat'],'Fdata');
        save([savedir 'DataHeader/' site '_' num2str(year) '_' num2str(currentday,'%03.0f') '_2_DataHeader.mat'],'Header');
    end
    
    clear Fdata Header tmpTEMPRMY
    
    for CW = 1:nw   % CW = Current Window ; nw = Number Of Windows per file
        
        fluxW=ORW_fluxW_process(fluxW,CW,NFavg,NSavg,ntsf,ntss,DespikeW,W1min,...
            MagneticCorrection, CSAT_angle,RMY_angle,...
            LI7500,Sdist_LI7500init,LI7700,Sdist_LI7700init,...
            Sonic_pathLength,Sonic_Height,AVGtime,dispH, freqF,nw, concatenating_files,...
            R_d);
        
    end %for CW (CW = Current Window)
    
    disp([savedir site '_' num2str(year) '_' num2str(currentday,'%03.0f') '_flux.mat; of raw file: ' fast.name])
    
    if CD~=1
        save([savedir 'flux/' site '_' num2str(year) '_' num2str(currentday,'%03.0f') '_flux.mat'],'fluxW');
    else
        save([savedir 'flux/' site '_' num2str(year) '_' num2str(currentday,'%03.0f') '_2_flux.mat'],'fluxW');
    end
    toc
end     % for the files in Data folder

if ~isempty(dir([savedir  'flux/' site '_' num2str(year) '_' num2str(currentday,'%03.0f') '_1_flux.mat'])) && ~isempty(dir([savedir 'flux/' site '_' num2str(year) '_' num2str(currentday,'%03.0f') '_2_flux.mat']))
    ORW126FluxJoiner;
end
