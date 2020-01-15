% Pre Allocating for all the variables in the process

% function [DespikeW,tmpTEMPRMY,fluxW,fluxWind,theta,alfa,u_CSAT,v_CSAT,w_CSAT,Fdata,Header,Tr_CSAT]...
%     =SOLPrealocatingVariables(nw,ntsf,NFavg,nfcol)
function [DespikeW,fluxW,fluxWind,theta,alfa,u_CSAT,v_CSAT,w_CSAT,Fdata,Header]...
    =SOLPrealocatingVariables(nw,ntsf,NFavg,nfcol)    
    
    NSavg = 30;
    
    % Measurments at 1/60 Hz
    %W1min.year   = nan(ntss,1);
    %W1min.doy    = nan(ntss,1);
    %W1min.dectime = nan(ntss,1);
    
    
    %W1min.Tm     = nan(ntss,1);       % HMP45C Temp [ Deg C ]
    %W1min.TmQA   = ones(ntss,1)*6;
    %W1min.Pvapor = nan(ntss,1);       % HMP45C Vapor pressure [ Pa ]
    %W1min.PvaporQA = ones(ntsf,1)*6;
    %W1min.Pvaporsat = nan(ntss,1);    % HMP45C Saturated Vapor Pressure [ ], Calculated. 611.2*exp(17.67*W1min.Tm./(243.51+W1min.Tm))
    %W1min.Hum    = nan(ntss,1);       % HMP45C Relative humidity [ % ]
    %W1min.HumQA  = ones(ntss,1)*6;
    %W1min.Q =      nan(ntss,1);       % HMP45C H2O Density [mol/m^3]  {Ideal gas law, n/V=Q=P/(RT)}

    %W1min.C1 = nan(ntss,1);  	      % ID200 CO2 reading at 2m [ppm]
    %W1min.C1QA = ones(ntss,1)*6;
    %W1min.C2 = nan(ntss,1);          % ID200 CO2 reading at 5m [ppm]
    %W1min.C2QA = ones(ntss,1)*6;
    %W1min.C3 = nan(ntss,1);         % ID200 CO2 reading at 10m [ppm]   
    %W1min.C3QA = nan(ntss,1)*6;
 
    %W1min.DS2_ubar = nan(ntss,1);        % DS2 wind speed [m/s]
    %W1min.DS2_uQA = nan(ntss,1)*6;
    %W1min.DS2_dir = nan(ntss,1);      % DS2 wind direction [deg off north, clockwise]
    %W1min.DS2_gust = nan(ntss,1);     % DS2 gust speed [m/s]
    %W1min.DS2_gustQA = nan(ntss,1)*6;
    
    %W1min.Albedo = nan(ntss,1);       % NR01 Albedo, calculated in calculate_radiation subrutine.
    %W1min.IRDown = nan(ntss,1);       % NR01 Net Long (Infra red) radiation down, calculated in calculate_radiation subrutine.
    %W1min.IRDownQA = ones(ntss,1)*6;
    %W1min.IRUp   = nan(ntss,1);       % NR01 Net Long (Infra red) radiation up, calculated in calculate_radiation subrutine.
    %W1min.IRUpQA   = ones(ntss,1)*6;
    %W1min.NetRs  = nan(ntss,1);       % NR01 Net radiation Short wave, calculated in calculate_radiation subrutine.
    %W1min.NetRl  = nan(ntss,1);       % NR01 Net radiation Long wave, calculated in calculate_radiation subrutine.
    %W1min.UpTot  = nan(ntss,1);       % NR01 Net radiation up, calculated in calculate_radiation subrutine.
    %W1min.DownTot= nan(ntss,1);       % NR01 Net radiation down, calculated in calculate_radiation subrutine.
    %W1min.NetRad = nan(ntss,1);       % NR01 Net radiation, calculated in calculate_radiation subrutine.
    %W1min.SWDown = nan(ntss,1);       % NR01 SR01Dn Short Radiation Down [ W/m^2 ]
    %W1min.SWDownQA = ones(ntss,1)*6;
    %W1min.SWUp   = nan(ntss,1);       % NR01 SR01UP Short Radiation Up [ W/m^2 ]
    %W1min.SWUpQA   = ones(ntss,1)*6;
    %W1min.Tcnr01 = nan(ntss,1);       % NR01 TCC Temp [ Deg C ]
    %W1min.Tcnr01QA = ones(ntss,1)*6;
    
    %W1min.Total_PAR  = nan(ntss,1);     % Direct/diffuse PAR sensor - Total PAR [mV]
    %W1min.Total_PARQA  = ones(ntss,1)*6;
    %W1min.Diffuse_PAR= nan(ntss,1);     % Direct/diffuse PAR sensor - Diffuse PAR  [mV]
    %W1min.Diffuse_PARQA= ones(ntsf,1)*6;
 

    %W1min.STat8cmOW_W1  = nan(ntss,1);     % 107-L Soil temp sensor [ Deg C ]
    %W1min.STat25cmOW_W1 = nan(ntss,1);     % 107-L Soil temp sensor [ Deg C ]
    %W1min.STat8cmOWr_W1 = nan(ntss,1);     % 107-L Soil temp sensor [ Deg C ]
    %W1min.STat8cmIM_W1  = nan(ntss,1);     % 107-L Soil temp sensor [ Deg C ]
    %W1min.STat25cmIM_W1 = nan(ntss,1);     % 107-L Soil temp sensor [ Deg C ]
    %W1min.STat8cmIMr_W1 = nan(ntss,1);     % 107-L Soil temp sensor [ Deg C ]
    %W1min.STat8cmUL_W1  = nan(ntss,1);     % 107-L Soil temp sensor [ Deg C ]
    %W1min.STat8cmOW_W2  = nan(ntss,1);     % 107-L Soil temp sensor [ Deg C ]
    %W1min.STat25cmOW_W2 = nan(ntss,1);     % 107-L Soil temp sensor [ Deg C ]
    %W1min.STat8cmOWr_W2 = nan(ntss,1);     % 107-L Soil temp sensor [ Deg C ]
    %W1min.STat8cmIM_W2  = nan(ntss,1);     % 107-L Soil temp sensor [ Deg C ]
    %W1min.STat25cmIM_W2 = nan(ntss,1);     % 107-L Soil temp sensor [ Deg C ]
    %W1min.STat8cmIMr_W2 = nan(ntss,1);     % 107-L Soil temp sensor [ Deg C ]
    %W1min.STat8cmUL_W2  = nan(ntss,1);     % 107-L Soil temp sensor [ Deg C ]
    
    %W1min.STat8cmOW_W1QA  = ones(ntss,1)*6;
    %W1min.STat25cmOW_W1QA = ones(ntss,1)*6; 
    %W1min.STat8cmOWr_W1QA = ones(ntss,1)*6; 
    %W1min.STat8cmIM_W1QA  = ones(ntss,1)*6;
    %W1min.STat25cmIM_W1QA = ones(ntss,1)*6;
    %W1min.STat8cmIMr_W1QA = ones(ntss,1)*6;
    %W1min.STat8cmUL_W1QA  = ones(ntss,1)*6;
    %W1min.STat8cmOW_W2QA  = ones(ntss,1)*6;
    %W1min.STat25cmOW_W2QA = ones(ntss,1)*6;
    %W1min.STat8cmOWr_W2QA = ones(ntss,1)*6;
    %W1min.STat8cmIM_W2QA  = ones(ntss,1)*6;
    %W1min.STat25cmIM_W2QA = ones(ntss,1)*6;
    %W1min.STat8cmIMr_W2QA = ones(ntss,1)*6;
    %W1min.STat8cmUL_W2QA  = ones(ntss,1)*6;
    
    %W1min.DS2_ubar = nan(ntss,1);
    %W1min.DS2_dir  = nan(ntss,1);
    %W1min.DS2_u    = nan(ntss,1);
    %W1min.DS2_v    = nan(ntss,1);
    
    %W1min.DS2_uQA = ones(ntss,1)*6;

    %W1min.WD_min = nan(ntss,1);
    %W1min.WD_bar = nan(ntss,1);
    %W1min.WD_max = nan(ntss,1);
    %W1min.WS_min  = nan(ntss,1);
    %W1min.WS_bar  = nan(ntss,1);
    %W1min.WS_max  = nan(ntss,1);
    %W1min.tair_2  = nan(ntss,1);
    %W1min.rH_2    = nan(ntss,1);
    %W1min.pres    = nan(ntss,1);
    %W1min.rain_accum  = nan(ntss,1);
    %W1min.rain_dur    = nan(ntss,1);
    %W1min.rain_intens = nan(ntss,1);
    %W1min.hail_accum  = nan(ntss,1);
    %W1min.hail_dur    = nan(ntss,1);
    %W1min.hail_intens = nan(ntss,1);
    
    %W1min.WS_min_QA = ones(ntss,1)*6;
    %W1min.WS_bar_QA = ones(ntss,1)*6;
    %W1min.WS_max_QA = ones(ntss,1)*6;
    %W1min.tair_2_QA = ones(ntss,1)*6;
    %W1min.rH_2_QA   = ones(ntss,1)*6;
    %W1min.pres_QA   = ones(ntss,1)*6;
    
    % Measurments at 10 Hz
    DespikeW.Time = nan(ntsf,3);
    DespikeW.ChkSumDataLog= nan(ntsf,1);% Check Sum for CR3000 DataLogger
    
    DespikeW.CSAT_U = nan(ntsf,1);       % CSAT Wind in U direction
    DespikeW.CSAT_UQA = ones(ntsf,1)*6;
    DespikeW.CSAT_V = nan(ntsf,1);       % CSAT Wind in V direction
    DespikeW.CSAT_VQA = ones(ntsf,1)*6;
    DespikeW.CSAT_W = nan(ntsf,1);       % CSAT Wind in W direction
    DespikeW.CSAT_WQA = ones(ntsf,1)*6;
    DespikeW.CSAT_Tmp = nan(ntsf,1);     % CSAT Sonic Temperature
    DespikeW.CSAT_TmpQA = ones(ntsf,1)*6;

%     DespikeW.CSAT_U2 = nan(ntsf,1);       % CSAT Wind in U direction
%     DespikeW.CSAT_UQA2 = ones(ntsf,1)*6;
%     DespikeW.CSAT_V2 = nan(ntsf,1);       % CSAT Wind in V direction
%     DespikeW.CSAT_VQA2 = ones(ntsf,1)*6;
%     DespikeW.CSAT_W2 = nan(ntsf,1);       % CSAT Wind in W direction
%     DespikeW.CSAT_WQA2 = ones(ntsf,1)*6;
%     DespikeW.CSAT_Tmp2 = nan(ntsf,1);     % CSAT Sonic Temperature
%     DespikeW.CSAT_TmpQA2 = ones(ntsf,1)*6;

    DespikeW.LI7500_C = nan(ntsf,1);    % LI7500 CO2 [mmol/mol]
    DespikeW.LI7500_CQA = ones(ntsf,1)*6;
    DespikeW.LI7500_Q = nan(ntsf,1);    % LI7500 H2O [mmol/mol]
    DespikeW.LI7500_QQA = ones(ntsf,1)*6;
    DespikeW.LI7500_P = nan(ntsf,1);    % LI7500 Pressure [ Pa ]
    DespikeW.LI7500_PQA = ones(ntsf,1)*6;

%     DespikeW.LI7500_C2 = nan(ntsf,1);    % LI7500 CO2 [mmol/mol]
%     DespikeW.LI7500_CQA2 = ones(ntsf,1)*6;
%     DespikeW.LI7500_Q2 = nan(ntsf,1);    % LI7500 H2O [mmol/mol]
%     DespikeW.LI7500_QQA2 = ones(ntsf,1)*6;
%     DespikeW.LI7500_P2 = nan(ntsf,1);    % LI7500 Pressure [ Pa ]
%     DespikeW.LI7500_PQA2 = ones(ntsf,1)*6;
    
    DespikeW.LI7700_M = nan(ntsf,1);    % LI7700 CH4
    DespikeW.LI7700_MQA = ones(ntsf,1)*6;
    DespikeW.LI7700_P = nan(ntsf,1);    % LI7700 Pressure [ Pa ]
    DespikeW.LI7700_PQA = ones(ntsf,1)*6;

    %fluxW
    
    SQAmat=zeros(nw,6);
    SQAmat(:,[1 6])=NSavg;
    FQAmat=zeros(nw,6);
    FQAmat(:,[1 6])=NFavg;
    
%     fluxW.tair  = nan(nw,1);            % HMP45C Mean of Temp in window period. [ Deg C ]
%     fluxW.tairQA = SQAmat;
%     fluxW.humbar= nan(nw,1);            % HMP45C Mean of humidity ( W1min.Hum ) [%]
%     fluxW.humbarQA=SQAmat;
%     fluxW.p_vapor_bar  = nan(nw,1);     % HMP45C Mean of Vapor Pressure ( W1min.Pvapor ) [Pa]
%     fluxW.p_vapor_barQA  = SQAmat;
%     fluxW.p_vaporSat_bar = nan(nw,1);   % HMP45C Saturated vapor pressure [Pa]

    fluxW.e = nan(nw,1);                % fluxW.pressure = mean(DspkW.LI7500_P), or default
                                        % HMP45C & LI7500 Saturated vapor pressure [Pa] from LI7500
                                        % e = mean(DespikeW.LI7500_Q) * 0.01802.* (tair+273.15)/ 0.21667 /10; 
                                        % Mean of LI7500 Vapor Pressure [Pa]
    fluxW.pressure = nan(nw,1);         % LI7500 mean Pressure. fluxW.pressure = mean( DespikeW.LI7500_P ) [Pa]
    fluxW.pressureQA = FQAmat;
    
    
%     fluxW.LWinbar  = nan(nw,1);         % NR01 Mean Long wave in radiation [ W/m^2 ]
%     fluxW.LWinbarQA  = SQAmat;
%     fluxW.LWoutbar = nan(nw,1);         % NR01 Mean Long wave out radiation [ W/m^2 ]
%     fluxW.LWoutbarQA = SQAmat;
%     fluxW.SWinbar  = nan(nw,1);         % NR01 Mean Short wave in radiation [ W/m^2 ]
%     fluxW.SWinbarQA  = SQAmat;
%     fluxW.SWoutbar = nan(nw,1);         % NR01 Mean Short wave out radiation [ W/m^2 ]
%     fluxW.SWoutbarQA = SQAmat;
%     fluxW.UpTot = nan(nw,1);            % NR01 Mean Radiation up [ W/m^2 ]
%     fluxW.DownTot =nan(nw,1);           % NR01 Mean Radiation down [ W/m^2 ]
%     fluxW.NetRad= nan(nw,1);            % NR01 Mean Net Radiation [ W/m^2 ]
%     fluxW.NetRs = nan(nw,1);            % NR01 Mean Net radiation short wave [ W/m^2 ]
%     fluxW.NetRl = nan(nw,1);            % NR01 Mean Net radiation long wave [ W/m^2 ]
%     fluxW.Albedo= nan(nw,1);            % NR01 Mean Albedo
% 
%     fluxW.Total_PAR     = nan(nw,1);    % BF5 mean Total PAR [mV]
%     fluxW.Total_PARQA   = SQAmat;
%     fluxW.Diffuse_PAR   = nan(nw,1);    % BF5 mean Total PAR [mV]
%     fluxW.Diffuse_PARQA = SQAmat;
%     
%     fluxW.STat8cmOW_W1 = nan(nw, 1);     % 107L soil temperature: 8cm into soil, in wetland 1, under open water   
%     fluxW.STat25cmOW_W1 = nan(nw, 1);    % 107L soil temperature: 25cm into soil, in wetland 1, under open water
%     fluxW.STat8cmOWr_W1 = nan(nw, 1);    % 107L soil temperature: 8cm into soil, in wetland 1, under open water replica
%     fluxW.STat8cmIM_W1 = nan(nw, 1);     % 107L soil temperature: 8cm into soil,in wetland 1, intermediate
%     fluxW.STat25cmIM_W1 = nan(nw, 1);    % 107L soil temperature: 25cm into soil, in wetland 1, intermediate
%     fluxW.STat8cmIMr_W1 = nan(nw, 1);    % 107L soil temperature: 8cm into soil, in wetland 1, intermediate replica
%     fluxW.STat8cmUL_W1 = nan(nw, 1);     % 107L soil temperature: 8cm into soil, in wetland 1, upland
%     fluxW.STat8cmOW_W2 = nan(nw, 1);     % 107L soil temperature: 8cm into soil, in wetland 2, under open water
%     fluxW.STat25cmOW_W2 = nan(nw, 1);    % 107L soil temperature: 25cm into soil, in wetland 2, under open water
%     fluxW.STat8cmOWr_W2 = nan(nw, 1);    % 107L soil temperature: 8cm into soil, in wetland 2, under open water replica
%     fluxW.STat8cmIM_W2 = nan(nw, 1);     % 107L soil temperature: 8cm into soil,in wetland 2, intermediate
%     fluxW.STat25cmIM_W2 = nan(nw, 1);    % 107L soil temperature: 25cm into soil, in wetland 2, intermediate
%     fluxW.STat8cmIMr_W2 = nan(nw, 1);    % 107L soil temperature: 8cm into soil, in wetland 2, intermediate replica
%     fluxW.STat8cmUL_W2 = nan(nw, 1);     % 107L soil temperature: 8cm into soil, in wetland 2, upland
%     fluxW.STat8cmOW_W1QA = SQAmat;
%     fluxW.STat25cmOW_W1QA = SQAmat;   
%     fluxW.STat8cmOWr_W1QA = SQAmat;   
%     fluxW.STat8cmIM_W1QA = SQAmat;    
%     fluxW.STat25cmIM_W1QA = SQAmat;   
%     fluxW.STat8cmIMr_W1QA = SQAmat;   
%     fluxW.STat8cmUL_W1QA = SQAmat;    
%     fluxW.STat8cmOW_W2QA = SQAmat;    
%     fluxW.STat25cmOW_W2QA = SQAmat;   
%     fluxW.STat8cmOWr_W2QA = SQAmat;   
%     fluxW.STat8cmIM_W2QA = SQAmat;    
%     fluxW.STat25cmIM_W2QA = SQAmat;   
%     fluxW.STat8cmIMr_W2QA = SQAmat;   
%     fluxW.STat8cmUL_W2QA = SQAmat;   
%     
%     fluxW.RMY_UQA=FQAmat;
%     fluxW.RMY_VQA=FQAmat;
%     fluxW.RMY_WQA=FQAmat;
%     fluxW.RMY_TmpQA=FQAmat;
    
    fluxW.CSAT_UQA=FQAmat;
    fluxW.CSAT_VQA=FQAmat;
    fluxW.CSAT_WQA=FQAmat;
    fluxW.CSAT_TmpQA=FQAmat;

%     fluxW.CSAT_UQA2=FQAmat;
%     fluxW.CSAT_VQA2=FQAmat;
%     fluxW.CSAT_WQA2=FQAmat;
%     fluxW.CSAT_TmpQA2=FQAmat;
    
%     fluxW.CO2_1 = nan(nw,1);         % ID200 CO2 sensor [ppm]
%     fluxW.CO2_2 = nan(nw,1);         % ID200 CO2 sensor [ppm]
%     fluxW.CO2_3 = nan(nw,1);         % ID200 CO2 sensor [ppm]
%     fluxW.C1QA = FQAmat;
%     fluxW.C2QA = FQAmat;
%     fluxW.C3QA = FQAmat;

    fluxW.LI7500_CQA = FQAmat;
    fluxW.LI7500_QQA = FQAmat;
    fluxW.LI7500_PQA = FQAmat;

%     fluxW.LI7500_CQA2 = FQAmat;
%     fluxW.LI7500_QQA2 = FQAmat;
%     fluxW.LI7500_PQA2 = FQAmat;
    
    fluxW.LI7700_MQA = FQAmat;
    fluxW.LI7700_PQA = FQAmat;
    
%     fluxW.ubar_2  = nan(nw,1);       % Lower sonic Mean U [m/s]
% %     fluxW.ubar_1  = nan(nw,1);       % Upper sonic Mean U [m/s]
%     fluxW.wbar_2  = nan(nw,1);       % Lower sonic Mean W [m/s]
% %     fluxW.wbar_1  = nan(nw,1);       % Upper sonic Mean W [m/s]
% %     fluxW.ustar = nan(nw,1);         % RMY or CSAT U star = (if RMY is nan than CSAT) (fluxW.uw^2+fluxW.vw^2)^0.25
%     fluxW.ustar_2 = nan(nw,1);       % Lower sonic U star = (fluxW.uw^2+fluxW.vw^2)^0.25
% %     fluxW.ustar_1 = nan(nw,1);       % Upper sonic U star = (fluxW.uw^2+fluxW.vw^2)^0.25
%     fluxW.uu_2 = nan(nw,1);          % Lower sonic Variance from mean u
% %     fluxW.uu_1    = nan(nw,1);       % Upper sonic Variance from mean u
%     fluxW.uw_2 = nan(nw,1);          % Lower sonic mean of (w.*un)
% %     fluxW.uw_1 = nan(nw,1);          % Upper sonic mean of (w.*un)
%     fluxW.u3_2 = nan(nw,1);          % Lower sonic third moment of u
% %     fluxW.u3_1 = nan(nw,1);          % Upper sonic moment of u
%     fluxW.vv_2 = nan(nw,1);          % Lower sonic Variance from mean v (After rotation mean u = 0)
% %     fluxW.vv_1 = nan(nw,1);          % Upper sonic Variance from mean v (After rotation mean u = 0)
%     fluxW.vw_2 = nan(nw,1);          % Lower sonic mean of (w.*v)
% %     fluxW.vw_1 = nan(nw,1);          % Upper sonic mean of (w.*v)
%     fluxW.v3_2 = nan(nw,1);          % Lower sonic third moment of v
% %     fluxW.v3_1 = nan(nw,1);          % Upper sonic moment of v
%     fluxW.ww_2 = nan(nw,1);          % Lower sonic Variance from mean w (After rotation mean w = 0)
% %     fluxW.ww_1 = nan(nw,1);          % Upper sonic Variance from mean w (After rotation mean w = 0)
%     fluxW.w3_2 = nan(nw,1);          % Lower sonic third moment of w
% %     fluxW.w3_1 = nan(nw,1);          % Upper sonic moment of w
%     fluxW.WD_2    = nan(nw,1);       % Lower sonic theta. calculated in RotateWind3D
%     fluxW.WDz_2   = nan(nw,1);       % Lower sonic alpha. Vertical rotation angle . Calculated in RotateWind3D
%     fluxW.WD_2_degN = nan(nw,1);     % RMY or CSAT (if CSAT is nan then RMY)
% %     fluxW.WD_1    = nan(nw,1);       % Upper sonic theta. calculated in RotateWind3D
% %     fluxW.WDz_1   = nan(nw,1);       % Upper sonic alpha. Vertical rotation angle . Calculated in RotateWind3D
%     fluxW.Sonic_Tmp_2 =nan(nw,1);    % Lower sonic Sonic Temperature
% %     fluxW.Sonic_Tmp_1 =nan(nw,1);    % Upper sonic Sonic Temperature
    
    
    
    fluxW.r    = nan(nw,1);         % HMP45C & LI7500 Mixing ratio r = 0.622 * p_vapor_bar/( pressure - p_vapor_bar),
    fluxW.tvair = nan(nw,1);            % HMP45C Mean of virtual temp calculated : fluxW.tair*(1+0.61*fluxW.r)
    fluxW.rho_dry_air = nan(nw,1);  % HMP45C & LI7500 Density of dry air
    fluxW.rho_moist_air = nan(nw,1);% HMP45C & LI7500 Density of moist air
    
    fluxW.Tr3    = nan(nw,1);           % Lower sonic & LI7500 moment(Tr,3)
    fluxW.Tr_bar  = nan(nw,1);          % Lower sonic & LI7500 Real Temperature. fluxW.Tr_bar= mean(Tr) | Tr = KaimalGaynor1990_10Hz( DespikeW.LI7500_P,Q,DespikeW.CSAT_Tmp,LI7500)
    fluxW.tt    = nan(nw,1);            % Lower sonic & LI7500 fluxW.tt = var(Tr). Tr = Real Temp (after KaimalGaynor)
    fluxW.rho_cp = nan(nw,1);       % HMP45C & LI7500 cp_moist*rho_v./1000 + cp_dry*rho_a./1000;
    fluxW.wtsonic= nan(nw,1);           % RMY or CSAT mean((Tson - Tson_mean).*w) Covariance between sonic t and w
    fluxW.wts   = nan(nw,1);            % Lower sonic & LI7500 fluxW.wts=mean(( Tr - fluxW.Tr_bar ).*w_Lower sonic)                                % calculate_rho_cp(tair ,pressure, p_vapor_bar, p_vaporSat_bar)
    fluxW.L     = nan(nw,1);            % Lower sonic & LI7500 & HMP45C Obukhov Length -(fluxW.ustar^3)/(0.4*9.81*(fluxW.wts/(fluxW.tvair+273.15)))
    fluxW.H     = nan(nw,1);            % Lower sonic & LI7500 fluxW.H=mean(w.*(Tr-nanmean(Tr))).*rho*1004.67 [rho, C, M] = WPL_correction_10Hz
    fluxW.Spectral_correction = nan(nw,1); % Lower sonic & LI7500 & HMP45C Spectral correction ratio, fluxes need to be divided by this 

    fluxW.lagq  = nan(nw,1);            % Lower sonic & LI7500 Lag between Sonic anemometer and CO2 from LI7500
    fluxW.lagc  = nan(nw,1);            % Lower sonic & LI7500 Lag between Sonic anemometer and H2O from LI7500
    fluxW.wqraw= nan(nw,1);             % Lower sonic & LI7500 Raw w'q' before wpl correction [mmol/m^2/s]
    fluxW.Q = nan(nw,1);
    fluxW.C = nan(nw,1);
    fluxW.qbar  = nan(nw,1);            % Lower sonic & LI7500 fluxW.qbar = mean( fluxW.qWPL )
    fluxW.qq    = nan(nw,1);            % Lower sonic & LI7500 fluxW.qq = var( fluxW.qWPL )
    fluxW.cbar  = nan(nw,1);            % Lower sonic & LI7500 fluxW.cbar = mean( fluxW.cWPL )
    fluxW.cc    = nan(nw,1);            % Lower sonic & LI7500 fluxW.cc = var( fluxW.cWPL )
    fluxW.cbarmpm3 = nan(nw,1);         % Lower sonic & LI7500 CO2 concentration mmol/m^3 after wpl 
    fluxW.wq    = nan(nw,1);            % Lower sonic & LI7500 mean(w.*qn)  ([u ,v, w ,alfa, theta] = RotateWind3D, qn = Q-fluxW.qbar )   fluxW.lagm  = nan(nw,1);            % Lower sonic & LI7700 Lag between Sonic anemometer and CH4 from LI7700
    fluxW.FH2O = nan(nw,1);
    fluxW.LE    = nan(nw,1);            % Lower sonic & LI7500 fluxW.LE=fluxW.wq*2440000 - unit change    
    fluxW.wc    = nan(nw,1);            % Lower sonic & LI7500 mean(w.*cn)  ( cn = C-fluxW.cbar )
    fluxW.Fc    = nan(nw,1);            % Lower sonic & LI7500 fluxW.Fc=fluxW.wc/44/1000*1000000 - unit change
    fluxW.c3    = nan(nw,1);            % Lower sonic & LI7500 fluxW.c3 = moment( fluxW.cWPL , 3 )
    fluxW.q3    = nan(nw,1);            % Lower sonic & LI7500 fluxW.q3 = moment( fluxW.qWPL , 3 )

    
%     fluxW.lagq3 = nan(nw,1);
%     fluxW.lagc3 = nan(nw,1);
%     fluxW.Q3 = nan(nw,1);
%     fluxW.C3 = nan(nw,1);
%     fluxW.wqraw3 = nan(nw,1);
%     fluxW.qbar3 = nan(nw,1);
%     fluxW.cbar3 = nan(nw,1);
%     fluxW.qq3 = nan(nw,1);
%     fluxW.cc3 = nan(nw,1);
%     fluxW.cbarmpm3_3 = nan(nw,1);
%     fluxW.wq_3 = nan(nw,1);
%     fluxW.FH2O_3 = nan(nw,1);
%     fluxW.LE_3 = nan(nw,1);
%     fluxW.wc_3 = nan(nw,1);
%     fluxW.Fc_3 = nan(nw,1);
% 
%    fluxW.ubar_3 = nan(nw,1);
%    fluxW.wbar_3 = nan(nw,1);
%    fluxW.uu_3 = nan(nw,1);
%    fluxW.vv_3 = nan(nw,1);
%    fluxW.ww_3 = nan(nw,1);
%    fluxW.u3_3 = nan(nw,1);
%    fluxW.v3_3 = nan(nw,1);
%    fluxW.w3_3 = nan(nw,1);
%    fluxW.uw_3 = nan(nw,1);
%    fluxW.vw_3 = nan(nw,1);
%    fluxW.ustar_3 = nan(nw,1);
%    fluxW.WD_3 = nan(nw,1);
%    fluxW.WDz_3 = nan(nw,1);
%    fluxW.Sonic_Tmp_3 = nan(nw,1);
%    fluxW.WD_3_degN = nan(nw,1);
   fluxW.CSAT_gust1 = nan(nw,1);
   fluxW.CSAT_gust2 = nan(nw,1);
   fluxW.CSAT_gust3 = nan(nw,1);
%    fluxW.DS2_u = nan(nw,1);
%    fluxW.DS2_dir = nan(nw,1);
%    fluxW.DS2_gust = nan(nw,1);
   fluxW.Spectral_correctionLI7500 = nan(nw,1);
%    fluxW.ubar_3 = nan(nw,1);
%    fluxW.wbar_3 = nan(nw,1);
%    fluxW.uu_3 = nan(nw,1);
%    fluxW.vv_3 = nan(nw,1);
%    fluxW.ww_3 = nan(nw,1);
%    fluxW.u3_3 = nan(nw,1);
%    fluxW.v3_3 = nan(nw,1);
%    fluxW.w3_3 = nan(nw,1);
%    fluxW.uw_3 = nan(nw,1);
%    fluxW.vw_3 = nan(nw,1);
%    fluxW.ustar_3 = nan(nw,1);
%    fluxW.WD_3 = nan(nw,1);
%    fluxW.WDz_3 = nan(nw,1);
%    fluxW.Sonic_Tmp_3 = nan(nw,1);
   fluxW.WD_3_degN = nan(nw,1);
   fluxW.CSAT_gust1 = nan(nw,1);
   fluxW.CSAT_gust2 = nan(nw,1);
   fluxW.CSAT_gust3 = nan(nw,1);
%    fluxW.DS2_u = nan(nw,1);
%    fluxW.DS2_dir = nan(nw,1);
%    fluxW.DS2_gust = nan(nw,1);
   fluxW.Spectral_correctionLI7500 = nan(nw,1);
   fluxW.Spectral_correctionLI7700 = nan(nw,1);
%    fluxW.Spectral_correctionEC150 = nan(nw,1);
   fluxW.Spectral_correctionLI7700 = nan(nw,1);
%    fluxW.Spectral_correctionEC150 = nan(nw,1);

    
    fluxW.wt2   = zeros(nw,1);          % Lower sonic & LI7500 fluxW.wt2 = fluxW.wt2 + (mean(w).*t)-mean(w).*mean(t))/J
    fluxW.wc2   = zeros(nw,1);          % Lower sonic & LI7500 fluxW.wc2 = fluxW.wc2 + (mean(w).*c)-mean(w).*mean(c))/J
    fluxW.wq2   = zeros(nw,1);          % Lower sonic & LI7500 fluxW.wq2 = fluxW.wq2 + (mean(w).*q)-mean(w).*mean(q))/J
    fluxW.wm2   = zeros(nw,1);          % Lower sonic & LI7700 fluxW.wm2 = fluxW.wm2 + (mean(w).*m)-mean(w).*mean(m))/J    
 
    % Storage terms - dx is averaged on a sub-interval (2D [sec]) at the beginning and the end of the interval window
    fluxW.dtdt  = nan(nw,1);            % Lower sonic Temp (convertSOS_Tc, Schotanus, De_spike3)
    fluxW.dqdt  = nan(nw,1);            % LI7500 CO2
    fluxW.dcdt  = nan(nw,1);            % LI7500 H2O
    fluxW.dmdt  = nan(nw,1);            % LI7700 CH4
    
    fluxW.lagm = nan(nw,1);
    fluxW.mA    = nan(nw,1);            % Lower sonic & LI7700 spectroscopic effects of temperature, pressure, and water vapor on methane density
    fluxW.mB    = nan(nw,1);            % Lower sonic & LI7700 spectroscopic corrections to the latent heat flux term for pressure and water vapor
    fluxW.mC    = nan(nw,1);            % Lower sonic & LI7700 spectroscopic corrections to the sensible heat flux term for temperature, pressure and water vapor
    fluxW.mbar  = nan(nw,1);            % Lower sonic & LI7700 it's the mean of CH4 flux after WPL and spectografic corrections
    fluxW.M = nan(nw,1);
    fluxW.wm    = nan(nw,1);            % Lower sonic & LI7700 mean(w.*mn)  ( mn=M-fluxW.mbar )
    fluxW.mm    = nan(nw,1);            % Lower sonic & LI7700 var(M)
    fluxW.FCH4 = nan(nw,1);
    fluxW.m3    = nan(nw,1);            % Lower sonic & LI7700 moment(M,3)

%     fluxW.WD_slow = nan(nw,1);          %WXT
%     fluxW.WS_slow = nan(nw,1);
%     fluxW.tair_2 = nan(nw,1);
%     fluxW.rH_2 = nan(nw,1);
%     fluxW.pres_2 = nan(nw,1);
%     fluxW.rain_accum = nan(nw,1);
%     fluxW.rain_dur = nan(nw,1);
%     fluxW.rain_intens = nan(nw,1);
    
    
    
    
    % Wind measurments
    fluxWind.ubar = nan(nw,1);
    fluxWind.wbar = nan(nw,1);
    fluxWind.uu = nan(nw,1);
    fluxWind.vv = nan(nw,1);
    fluxWind.ww = nan(nw,1);
%     fluxWind.u3 = nan(nw,1);
%     fluxWind.v3 = nan(nw,1);
%     fluxWind.w3 = nan(nw,1);
    fluxWind.uw = nan(nw,1);
    fluxWind.vw = nan(nw,1);
    fluxWind.ustar = nan(nw,1);
    theta = nan(nw,1);
    alfa = nan(nw,1);
    

    
    % Wind measurments for Lower Sonic
    u_CSAT = nan(NFavg,1);
    v_CSAT = nan(NFavg,1);
    w_CSAT = nan(NFavg,1);
    Tr_CSAT  = nan(nw,1);   % KaimalGaynor1990_10Hz( pressureUse, DespikeW.LI7500_Q(CFV), DespikeW.CSAT_Tmp(CFV), LI7500);  % CSAT & LI7500
    %Good_Data_CO2_H20 = nan(NFavg,1);
    qWPL  = nan(NFavg,1);               % Lower sonic & LI7500 Mean CO2 flux, after WPL correction.
    cWPL  = nan(NFavg,1);               % Lower sonic & LI7500 Mean H2O flux, after WPL correction.

% u = u_RMY;
%             v = v_RMY;
%             w = w_RMY;
%             Tr = Tr_RMY;
%             Tson = DespikeW.RMY_Tmp(CFV);
%             Tson_mean = fluxW.Sonic_Tmp_1(CW);
%             ubar = fluxW.ubar_1(CW);
%             Sdist_LI7500 = Sdist_LI7500init;
%             Sdist_LI7700 = Sdist_LI7700init;
% Q = nan(NFavg,1);
% C = nan(NFavg,1);
    
    mWPL  = nan(NFavg,1);               % Lower sonic & LI7700 Mean CH4 flux, after WPL correction.
    
    Header.S1 = nan;
    Header.S2 = nan;
    Header.S3 = nan;
    Header.S4 = nan;
    Header.F1 = nan;
    Header.F2 = nan;
    Header.F3 = nan;
    Header.F4 = nan;
    Header.DS1 = nan;
    Header.DS2 = nan;
    Header.DS3 = nan;
    Header.DS4 = nan;
    
 
    Fdata=nan(ntsf,nfcol);              %Fast data raw file
    
%     fluxW.WS_min_QA = SQAmat;
%     fluxW.WS_bar_QA = SQAmat;
%     fluxW.WS_max_QA = SQAmat;
%     fluxW.tair_2_QA = SQAmat;
%     fluxW.rH_2_QA   = SQAmat;
%     fluxW.pres_QA   = SQAmat;
    
%     fluxW.CSAT_UQA2 = FQAmat;
%     fluxW.CSAT_VQA2 = FQAmat;
%     fluxW.CSAT_WQA2 = FQAmat;
%     fluxW.CSAT_TmpQA2 = FQAmat;
%     fluxW.LI7500_CQA2 = FQAmat;
%     fluxW.LI7500_QQA2 = FQAmat;
%     fluxW.LI7500_PQA2 = FQAmat;
    
    
    
    
end


 
