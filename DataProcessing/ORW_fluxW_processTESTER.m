function fluxW=ORW_fluxW_process(fluxW,CW,NFavg,NSavg,ntsf,ntss,DespikeW,W1min,...
    MagneticCorrection, CSAT_angle,RMY_angle,...
    LI7500,Sdist_LI7500init,LI7700,Sdist_LI7700init,...
    Sonic_pathLength,Sonic_Height,AVGtime,dispH, freqF,nw, concatenating_files,...
    R_d,...
    varargin)
CFV = ( (CW-1)*NFavg+1 : min( CW*NFavg , ntsf ))'; % CFV = Current Fast Values, Make a vector of the fast values to work on for this Current Window
CSV = ( (CW-1)*NSavg+1 : min( NSavg*CW , ntss ))'; % CSV = Current Slow Values, Make a vector of the slow values to work on for this Current Window

LI7500path=0.125;                %m
LI7500dia=0.015;                 %m
LI7700path=0.475;                %m
LI7700dia=0.075;                   %m
EC150path=0.144;                 %m
EC150dia=0.015;                  %m


tmpTEMPRMY = nan(ntsf,1);
SQAmat=zeros(nw,6);
SQAmat(:,[1 6])=NSavg;
% Wind measurments
fluxWind.ubar = nan(nw,1);
fluxWind.wbar = nan(nw,1);
fluxWind.uu = nan(nw,1);
fluxWind.vv = nan(nw,1);
fluxWind.ww = nan(nw,1);
fluxWind.u3 = nan(nw,1);
fluxWind.v3 = nan(nw,1);
fluxWind.w3 = nan(nw,1);
fluxWind.uw = nan(nw,1);
fluxWind.vw = nan(nw,1);
fluxWind.ustar = nan(nw,1);
theta = nan(nw,1);
alfa = nan(nw,1);
w = nan(NFavg,1);
Tr_1  = nan(nw,1);   % KaimalGaynor1990_10Hz( pressureUse, DespikeW.LI7500_Q(CFV), DespikeW.CSAT_Tmp(CFV), LI7500);  % CSAT & LI7500
Tr_2  = nan(nw,1);
Tson_1 = nan(NFavg,1);
Tson_2 = nan(NFavg,1);
%Good_Data_CO2_H20 = nan(NFavg,1);
qWPL  = nan(NFavg,1);               %#ok<*NASGU> % Lower sonic & LI7500 Mean CO2 flux, after WPL correction.
cWPL  = nan(NFavg,1);               % Lower sonic & LI7500 Mean H2O flux, after WPL correction.
mWPL  = nan(NFavg,1);               % Lower sonic & LI7700 Mean CH4 flux, after WPL correction.
    

%%% Reading data from sonics
if isempty(varargin)
    upper_son='RMY';
    lower_son='CSAT';
elseif max(size(varargin))==1
    disp('only one sonic input')
    return
else
    upper_son=varargin(1);
    lower_son=varargin(2);
end


% from Upper Sonic
if ~strcmp(upper_son,'None') 
    if strcmp(upper_son,'RMY') % upper_son==RMY
        Good_Data_T_W_RMY = find( ~isnan( DespikeW.RMY_Tmp( CFV )) & ~isnan( DespikeW.RMY_W( CFV )));
        
        if length(Good_Data_T_W_RMY) > 0.5 * NFavg   % NFavd = Number of fast values that are meaned into one value per window
            
            [fluxWind, ~, ~, w_2, alfa, theta] = ...
                WindData( DespikeW.RMY_U(CFV), DespikeW.RMY_V(CFV), DespikeW.RMY_W(CFV),'RMY', MagneticCorrection, RMY_angle);
            
            % Wind measurments for Upper Sonic
            fluxW.ubar_2(CW) = fluxWind.ubar;       % RMYoung Mean U
            fluxW.wbar_2(CW) = fluxWind.wbar;       % RMYoung Mean W
            fluxW.uu_2(CW) = fluxWind.uu;           % RMYoung Variance from mean u
            fluxW.vv_2(CW) = fluxWind.vv;           % RMYoung Variance from mean v (After rotation mean u = 0)
            fluxW.ww_2(CW) = fluxWind.ww;           % RMYoung Variance from mean w (After rotation mean w = 0)
            fluxW.u3_2(CW) = fluxWind.u3;           % RMYoung moment of u
            fluxW.v3_2(CW) = fluxWind.v3;           % RMYoung moment of v
            fluxW.w3_2(CW) = fluxWind.w3;           % RMYoung moment of w
            fluxW.uw_2(CW) = fluxWind.uw;           % RMYoung mean of (w.*un)
            fluxW.vw_2(CW) = fluxWind.vw;           % RMYoung mean of (w.*v)
            fluxW.ustar_2(CW) = fluxWind.ustar;     % RMYoung U star = (fluxW.uw^2+fluxW.vw^2)^0.25
            fluxW.WD_2(CW)    = theta;      % RMYoung theta. calculated in RotateWind3D
            fluxW.WDz_2(CW)   = alfa;       % RMYoung alpha. Vertical rotation angle . Calculated in RotateWind3D
            fluxW.Sonic_Tmp_2(CW) = nanmean( DespikeW.RMY_Tmp(CFV)); % RMYoung Mean Sonic temp [Deg C]
            fluxW.WD_2_degN = mod(fluxW.WD_2*180/pi,360);
        else
            w_2=nan(NFavg,1);
        end
    end
end

% from Lower Sonic
if strcmp(lower_son,'CSAT') % lower_son==CSAT 
    Good_Data_T_W = find( ~isnan( DespikeW.CSAT_Tmp( CFV )) & ~isnan( DespikeW.CSAT_W( CFV )));
    
    if length(Good_Data_T_W) > 0.5 * NFavg   % NFavd = Number of fast values that are meaned into one value per window
        
        [fluxWind, ~, ~, w_1, alfa, theta] = ...
            WindData( DespikeW.CSAT_U(CFV), DespikeW.CSAT_V(CFV), DespikeW.CSAT_W(CFV),'CSAT', MagneticCorrection, CSAT_angle);
        
        % Wind measurments for Lower Sonic
        fluxW.ubar_1(CW) = fluxWind.ubar;      % Lower Sonic Mean U
        fluxW.wbar_1(CW) = fluxWind.wbar;      % Lower Sonic Mean W
        fluxW.uu_1(CW) = fluxWind.uu;          % Lower Sonic Variance from mean u
        fluxW.vv_1(CW) = fluxWind.vv;          % Lower Sonic Variance from mean v (After rotation mean u = 0)
        fluxW.ww_1(CW) = fluxWind.ww;          % Lower Sonic Variance from mean w (After rotation mean w = 0)
        fluxW.u3_1(CW) = fluxWind.u3;          % Lower Sonic moment of u
        fluxW.v3_1(CW) = fluxWind.v3;          % Lower Sonic moment of v
        fluxW.w3_1(CW) = fluxWind.w3;          % Lower Sonic moment of w
        fluxW.uw_1(CW) = fluxWind.uw;          % Lower Sonic mean of (w.*un)
        fluxW.vw_1(CW) = fluxWind.vw;          % Lower Sonic mean of (w.*v)
        fluxW.ustar_1(CW) = fluxWind.ustar;    % Lower Sonic U star = (fluxW.uw^2+fluxW.vw^2)^0.25
        fluxW.WD_1(CW)    = theta;      % Lower Sonic theta. calculated in RotateWind3D
        fluxW.WDz_1(CW)   = alfa;       % Lower Sonic alpha. Vertical rotation angle . Calculated in RotateWind3D
        fluxW.Sonic_Tmp_1(CW) = nanmean( DespikeW.CSAT_Tmp(CFV)); % Lower Sonic Mean Sonic temp [Deg C]
        fluxW.WD_1_degN = mod(fluxW.WD_1*180/pi,360);
    end
elseif strcmp(lower_son,'RMY') % lower_son==RMY
    Good_Data_T_W_RMY = find( ~isnan( DespikeW.RMY_Tmp( CFV )) & ~isnan( DespikeW.RMY_W( CFV )));
    
    if length(Good_Data_T_W_RMY) > 0.5 * NFavg   % NFavd = Number of fast values that are meaned into one value per window
        
        [fluxWind, ~, ~, w_1, alfa, theta] = ...
            WindData( DespikeW.RMY_U(CFV), DespikeW.RMY_V(CFV), DespikeW.RMY_W(CFV),'RMY', MagneticCorrection, RMY_angle);
        fluxW.ubar_1(CW) = fluxWind.ubar;      % Lower Sonic Mean U
        fluxW.wbar_1(CW) = fluxWind.wbar;      % Lower Sonic Mean W
        fluxW.uu_1(CW) = fluxWind.uu;          % Lower Sonic Variance from mean u
        fluxW.vv_1(CW) = fluxWind.vv;          % Lower Sonic Variance from mean v (After rotation mean u = 0)
        fluxW.ww_1(CW) = fluxWind.ww;          % Lower Sonic Variance from mean w (After rotation mean w = 0)
        fluxW.u3_1(CW) = fluxWind.u3;          % Lower Sonic moment of u
        fluxW.v3_1(CW) = fluxWind.v3;          % Lower Sonic moment of v
        fluxW.w3_1(CW) = fluxWind.w3;          % Lower Sonic moment of w
        fluxW.uw_1(CW) = fluxWind.uw;          % Lower Sonic mean of (w.*un)
        fluxW.vw_1(CW) = fluxWind.vw;          % Lower Sonic mean of (w.*v)
        fluxW.ustar_1(CW) = fluxWind.ustar;    % Lower Sonic U star = (fluxW.uw^2+fluxW.vw^2)^0.25
        fluxW.WD_1(CW)    = theta;      % Lower Sonic theta. calculated in RotateWind3D
        fluxW.WDz_1(CW)   = alfa;       % Lower Sonic alpha. Vertical rotation angle . Calculated in RotateWind3D
        fluxW.Sonic_Tmp_1(CW) = nanmean( DespikeW.RMY_Tmp(CFV)); % Lower Sonic Mean Sonic temp [Deg C]
        fluxW.WD_1_degN = mod(fluxW.WD_1*180/pi,360);
    end
end

if isfield(W1min,'DS2_u')
    CSATubar = sqrt(DespikeW.CSAT_U(CFV).^2 + DespikeW.CSAT_V(CFV).^2);
    fluxW.CSAT_gust1(CW) = max(CSATubar);
    fluxW.CSAT_gust2(CW) = prctile(CSATubar,90);
    fluxW.CSAT_gust3(CW) = max(moving(CSATubar,10));
    fluxW.DS2_u(CW) = nanmean(W1min.DS2_ubar(CSV));
    DS2_u=nanmean(W1min.DS2_u(CSV));
    DS2_v=nanmean(W1min.DS2_v(CSV));
    fluxW.DS2_dir(CW) = mod(180/pi*atan2(DS2_v,DS2_u)+360,360);
    fluxW.DS2_gust(CW) = nanmean(W1min.DS2_gust(CSV));
end

% Data from CO2 sensors
fluxW.CO2_1(CW) = nanmean(W1min.C1(CSV));			% Mean CO2 concentration at lower CO2 sensor
fluxW.CO2_2(CW) = nanmean(W1min.C2(CSV));                       % Mean CO2 concentration at mid CO2 sensor
fluxW.CO2_3(CW) = nanmean(W1min.C3(CSV));                       % Mean CO2 concentration at upper CO2 sensor

% Data from BF5
fluxW.Total_PAR(CW)  = nanmean(W1min.Total_PAR(CSV));             % Mean Total PAR [mV]
fluxW.Diffuse_PAR(CW)  =  nanmean(W1min.Diffuse_PAR(CSV));        % Direct/diffuse PAR sensor - Total PAR [mV]

% Data from Soil temp probes
fluxW.STat8cmOW_W1 (CW) = nanmean(W1min.STat8cmOW_W1(CSV));
fluxW.STat25cmOW_W1 (CW) = nanmean(W1min.STat25cmOW_W1(CSV));
fluxW.STat8cmOWr_W1 (CW) = nanmean(W1min.STat8cmOWr_W1(CSV));
fluxW.STat8cmIM_W1 (CW) = nanmean(W1min.STat8cmIM_W1(CSV));
fluxW.STat25cmIM_W1 (CW) = nanmean(W1min.STat25cmIM_W1(CSV));
fluxW.STat8cmIMr_W1 (CW) = nanmean(W1min.STat8cmIMr_W1(CSV));
fluxW.STat8cmUL_W1 (CW) = nanmean(W1min.STat8cmUL_W1(CSV));
fluxW.STat8cmOW_W2 (CW) = nanmean(W1min.STat8cmOW_W1(CSV));
fluxW.STat25cmOW_W2 (CW) = nanmean(W1min.STat25cmOW_W2(CSV));
fluxW.STat8cmOWr_W2 (CW) = nanmean(W1min.STat8cmOWr_W2(CSV));
fluxW.STat8cmIM_W2 (CW) = nanmean(W1min.STat8cmIM_W2(CSV));
fluxW.STat25cmIM_W2 (CW) = nanmean(W1min.STat25cmIM_W2(CSV));
fluxW.STat8cmIMr_W2 (CW) = nanmean(W1min.STat8cmIMr_W2(CSV));
fluxW.STat8cmUL_W2 (CW) = nanmean(W1min.STat8cmUL_W2(CSV));


%%% T Air
fluxW.tair(CW) = nanmean( W1min.Tm( CSV ));         % HMP45C mean air temp [Deg C]
% Defalt for temperature measurement
if isnan(fluxW.tair(CW))
    fluxW.tair(CW)=fluxW.Sonic_Tmp_1(CW);           % Default to CSAT temp if HMP45C temp is missing
    if isnan(fluxW.tair(CW))
        fluxW.tair(CW)=fluxW.Sonic_Tmp_2(CW);           % Default to RMY temp if CSAT temp is missing
    end
end

%%% Vapor Pressure - Mean of LI7500 Vapor Pressure [Pa]
if sum( ~isnan( DespikeW.LI7500_Q( CFV )))>0.5 * NFavg
    fluxW.e( CW )= nanmean(DespikeW.LI7500_Q(CFV)) * 0.01802.* (fluxW.tair(CW)+273.15)/ 0.21667 /10;
end

%%% Saturated Vapor Pressure - Mean of Saturation Vapor Pressure [Pa]
fluxW.p_vaporSat_bar( CW ) = nanmean(W1min.Pvaporsat(CSV));
if isnan(fluxW.p_vaporSat_bar( CW )) % Retries with a possibly updated temperature value [Pa]
    fluxW.p_vaporSat_bar(CW) = 611.2*exp(17.67*fluxW.tair(CW)/(243.51+fluxW.tair(CW)));
end

%%% Vapor Pressure - Mean of Vapor Pressure [Pa]
fluxW.p_vapor_bar( CW ) = nanmean(W1min.Pvapor(CSV));
if isnan(fluxW.p_vapor_bar( CW ))
    % fluxW.p_vapor_bar( CW )=nanmean(DespikeW.LI7500_Q(CFV)/1000)*R_g*(fluxW.tair(CW)+273.16);
    % From ideal gas law. Output in [Pa]
    fluxW.p_vapor_bar( CW )=fluxW.e(CW);
    if fluxW.p_vapor_bar( CW )>fluxW.p_vaporSat_bar( CW )
        fluxW.p_vapor_bar( CW )=fluxW.p_vaporSat_bar( CW );
    end
end

%%% Relative Humidity
fluxW.humbar(CW) = nanmean( W1min.Hum( CSV ));                  % HMP45C mean of hummidity [%]
if isnan(fluxW.humbar(CW))
    fluxW.humbar(CW)=fluxW.p_vapor_bar( CW )/fluxW.p_vaporSat_bar( CW )*100;
    %Comes in over 100% sometimes
    if fluxW.humbar(CW)>100 && fluxW.humbar(CW)<110
        fluxW.humbar(CW)=100;
    end
end

%%% Atmospheric Pressure
fluxW.pressure( CW ) = nanmean( DespikeW.LI7500_P( CFV ));      % LI7500 Pressure [Pa]
pressureUse = DespikeW.LI7500_P( CFV );
if isnan(fluxW.pressure( CW ))
    fluxW.pressure( CW ) = nanmean( DespikeW.LI7700_P( CFV ));
    pressureUse = DespikeW.LI7700_P( CFV );
end
if isnan(fluxW.pressure( CW ))
    pressureUse = nan(length(CFV),1);
    if CW > 1
        fluxW.pressure( CW ) = fluxW.pressure(CW-1);
        pressureUse(:) = fluxW.pressure(CW-1);
    else
        fluxW.pressure( CW ) = 9.9116e+04;
        pressureUse(:) = 9.9116e+04;
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fluxW.r(CW) = 0.622 * fluxW.p_vapor_bar( CW )/( fluxW.pressure( CW ) - fluxW.p_vapor_bar( CW )); % HMP45C & LI7500 Mixing ratio [unitless]

% Virtual temperature of a moist air parcel is the temperature at which a theoretical dry air parcel would have a total pressure and density equal to the moist parcel of air.
fluxW.tvair(CW) = fluxW.tair( CW ) * ( 1 + 0.61*fluxW.r(CW));  % HMP45C & LI7500 mean of virtual air temp [Deg C]

fluxW.rho_dry_air(CW)=fluxW.pressure(CW)/R_d/(fluxW.tair(CW)+273.15); % Dry Air density [kg/m^3]
% Rd = Specific gas constant for dry air, 287.058 J/(kg·K)
fluxW.rho_moist_air(CW)=(fluxW.pressure(CW)-0.3780*fluxW.p_vapor_bar(CW))/R_d/(fluxW.tair(CW)+273.15); % Moist air density [kg/m^3]
% Rv = Specific gas constant for water vapor, 461.495 J/(kg·K)

fluxW.SWinbar(CW) = nanmean(W1min.SWUp(CSV));    % NR01 mean Short wave down radiation
fluxW.SWoutbar(CW) = nanmean(W1min.SWDown(CSV));     % NR01 mean Short wave up radiation
fluxW.LWinbar(CW) = nanmean(W1min.IRUp(CSV));    % NR01 mean Long wave down radiation
fluxW.LWoutbar(CW) = nanmean(W1min.IRDown(CSV));     % NR01 mean Long wave up radiation
fluxW.NetRs(CW) = nanmean(W1min.NetRs(CSV));       % NR01 mean Short wave net radiation
fluxW.NetRl(CW)= nanmean(W1min.NetRl(CSV));        % NR01 mean Long wave net radiation
fluxW.Albedo(CW)= nanmean(W1min.Albedo(CSV));      % NR01 mean albedo
fluxW.UpTot(CW)= nanmean(W1min.UpTot(CSV));        % NR01 mean Total up radiation
fluxW.DownTot(CW)= nanmean(W1min.DownTot(CSV));    % NR01 mean Total down radiation
fluxW.NetRad(CW)= nanmean(W1min.NetRad(CSV));      % NR01 mean Total radiation

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Correction of Temperature from sonic. Error is caused due to H2O and Pressure.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if strcmp(lower_son,'RMY') % lower_son==RMY
    [Tr_1] = KaimalGaynor1990_10Hz( pressureUse, DespikeW.LI7500_Q(CFV), DespikeW.RMY_Tmp(CFV), LI7500);  % RMY & LI7500
    Tson_1 = DespikeW.RMY_Tmp(CFV);
else
    [Tr_1] = KaimalGaynor1990_10Hz( pressureUse, DespikeW.LI7500_Q(CFV), DespikeW.CSAT_Tmp(CFV), LI7500);  % CSAT & LI7500
    Tson_1 = DespikeW.CSAT_Tmp(CFV);
end
[Tr_2] = KaimalGaynor1990_10Hz( pressureUse, DespikeW.LI7500_Q(CFV), DespikeW.RMY_Tmp(CFV), LI7500);  % RMY & LI7500
Tson_2 = DespikeW.RMY_Tmp(CFV);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% If no wind data from lower sonic and there is an upper sonic, use upper
% sonic for wind data
if ~isnan( fluxW.ustar_1(CW) )
    fluxW.ustar(CW) = fluxW.ustar_1(CW);
    w = w_1;
    Tr = Tr_1;
    Tson = Tson_1;
    Tson_mean = fluxW.Sonic_Tmp_1(CW);
    ubar = fluxW.ubar_1(CW);
    Sdist_LI7500 = Sdist_LI7500init;
    Sdist_LI7700 = Sdist_LI7700init;
elseif strcmp(upper_son,'RMY') 
    fluxW.ustar(CW) = fluxW.ustar_2(CW);
    w = w_2;
    Tr = Tr_2;
    Tson = Tson_2;
    Tson_mean = fluxW.Sonic_Tmp_2(CW);
    ubar = fluxW.ubar_2(CW);
    Sdist_LI7500 = Sdist_LI7500init+2;
    Sdist_LI7700 = Sdist_LI7700init+1.6;
else
    Tr = nan(NFavg,1);
    w = nan(NFavg,1);
    Sdist_LI7500 = Sdist_LI7500init;
    Sdist_LI7700 = Sdist_LI7700init;
    ubar = nan;
    Tson = nan(NFavg,1);
    Tson_mean = nan;
end

fluxW.Tr3(CW)  = nanmoment( Tr, 3 );    % Sonic & LI7500 Real Temperature moment
fluxW.Tr_bar(CW) = nanmean( Tr ) ;      % Sonic & LI7500 Real Temperature
fluxW.tt( CW ) = nanvar( Tr);           % Sonic & LI7500 Variance of Tr = Real Temp (after KaimalGaynor)

fluxW.rho_cp(CW) = calculate_rho_cp(fluxW.tair(CW) ,fluxW.pressure( CW ),fluxW.p_vapor_bar( CW ),fluxW.p_vaporSat_bar( CW ));
fluxW.wtsonic(CW) = nanmean((Tson - Tson_mean).*w);
fluxW.wts(CW) = nanmean(( Tr - fluxW.Tr_bar(CW) ).*w);    % Sonic & LI7500 Temperature flux
fluxW.L(CW)= -(fluxW.ustar(CW)^3) / (0.4*9.81*( fluxW.wts(CW) / (fluxW.tair(CW)+273.15)));    % Sonic & LI7500 & HMP45C Obukov lengh

fluxW.H(CW) = fluxW.wts(CW).* fluxW.rho_cp(CW) ; % Sonic & LI7500 & HMP45C sensible Heat flux [W/m2]
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% --SCALARS--------------------------------
wsc=(w-nanmean(w))/nanstd(w);
wsc(isnan(wsc))=0;

fluxW.Spectral_correctionLI7500(CW) = opirgafreq(ubar,Sonic_pathLength,Sonic_Height,AVGtime,(Sonic_Height-dispH)/fluxW.L(CW),LI7500path,LI7500dia);
fluxW.Spectral_correctionLI7700(CW) = opirgafreq(ubar,Sonic_pathLength,Sonic_Height,AVGtime,(Sonic_Height-dispH)/fluxW.L(CW),LI7700path,LI7700dia);
fluxW.Spectral_correctionEC150(CW) = opirgafreq(ubar,Sonic_pathLength,Sonic_Height,AVGtime,(Sonic_Height-dispH)/fluxW.L(CW),EC150path,EC150dia);

% LI7500 : CO2 & H2O
Good_Data_CO2_H20 = find( ~isnan( DespikeW.LI7500_Q( CFV )) & ~isnan( DespikeW.LI7500_C( CFV )) & ~isnan( w ));

if length( Good_Data_CO2_H20 ) > 0.5 * NFavg
    
    % find lag time between LI7500 and Sonic
    
    Q=DespikeW.LI7500_Q(CFV);
    Qsc=(Q-nanmean(Q))/nanstd(Q);
    Qsc(isnan(Q))=0;
    C=DespikeW.LI7500_C(CFV);
    Csc=(C-nanmean(C))/nanstd(C);
    Csc(isnan(C))=0;
    
    fluxW.lagq( CW ) = lag_test( wsc, Qsc, LI7500, Sdist_LI7500, ubar, freqF);   % Calculate lag between CO2 and Wind
    fluxW.lagc( CW ) = lag_test( wsc, Csc, LI7500, Sdist_LI7500, ubar, freqF);   % Calculate lag between H2O and Wind
    
    %                 if abs(fluxW.lagq(CW)-fluxW.lagc(CW))>max_lag_diff_LI7500
    C = nan( NFavg , 1 ) ;
    Q = nan( NFavg , 1 ) ;
    % Calculate indices for vector shift CO2 & H2O
    % nw = number of windows per file [number/file]
    % ntsf =  # of timesteps per file - fast data
    [V_lag_1, V_lag_end , V_1 , V_end] = Lag_Shift_index (fluxW.lagc(CW) , concatenating_files , CW , ntsf/nw , nw );
    C( V_lag_1 : V_lag_end ) = DespikeW.LI7500_C( V_1 : V_end );
    [V_lag_1, V_lag_end , V_1 , V_end] = Lag_Shift_index (fluxW.lagq(CW) , concatenating_files , CW , ntsf/nw , nw );
    Q( V_lag_1 : V_lag_end ) = DespikeW.LI7500_Q( V_1 : V_end );
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    fluxW.Q(CW) = nanmean(Q);
    fluxW.C(CW) = nanmean(C);
    fluxW.wqraw(CW) = nanmean((Q-fluxW.Q(CW)).*w);
    
    % Aapply WPL correction on raw data. Error is caused by changes
    % in air volume deu to heat.
    %add condition regardind LI7700 isnan and if exist do wpl 4
    %both
    [~, qWPL, cWPL, ~] = WPL_correction_10Hz...
        ( Q, C,fluxW.tair(CW),Tson, DespikeW.LI7500_P(CFV), nan, LI7500, DespikeW.LI7700_P(CFV) );
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
    % evaluete statistics
    
    fluxW.qbar( CW ) = nanmean( qWPL ) ;   % Sonic & LI7500 mean of H2O concentration [kg/m^3]
    fluxW.cbar( CW ) = nanmean( cWPL ) ;   % Sonic & LI7500 mean of CO2 concentration [kg/m^3]
    fluxW.qq( CW ) = nanvar( qWPL );   % Sonic & LI7500 Variance of H2O concentration
    fluxW.cc( CW ) = nanvar( cWPL );   % Sonic & LI7500 Variance of CO2 concentration
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    fluxW.cbarmpm3( CW ) =  fluxW.cbar( CW )*1000*1000/44.0095 ;   % Mean of CO2 concentration [mmol/m^3]
    fluxW.wq( CW ) = nanmean( w.*( qWPL - fluxW.qbar( CW )) )/fluxW.Spectral_correctionLI7500(CW);   % H2O flux in :[kg/m2/s] (Water Vapor flux in : [kg/m2/s] units)
    fluxW.FH2O( CW )=fluxW.wq( CW )/(18.0153/1000)*1000; % [mmol/m2/s]
    % mol H2O = 18.0153 gr H2O     O=15.9994 [gr/mol] H=1.00794 [gr/mol]
    
    fluxW.LE( CW ) = fluxW.wq( CW )*2440000;	% --- (Latent Heat flux in : [W/m2] units)
    fluxW.wc( CW ) = nanmean( w.*( cWPL - fluxW.cbar( CW )))/fluxW.Spectral_correctionLI7500(CW);   % CO2 flux in : [kg/m2/s] units
    fluxW.Fc( CW ) = fluxW.wc( CW )/(44.0095/1000)*1000000;  % CO2 flux in : [umol/m2/s] units
    % mol CO2 = 44.0095 gr CO2     C=12.0107 [gr/mol] O=15.9994 [gr/mol]
    
    fluxW.q3( CW ) = nanmoment(qWPL,3); % moment( qWPL, 3 );    % Sonic & LI7500 Moment of H2O concentration
    fluxW.c3( CW ) = nanmoment(cWPL,3); % moment( cWPL, 3 );    % Sonic & LI7500 Moment of CO2 concentration
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % Compute fluxes in specific interval (Sp_Int) average for stationary test
    % NFavg = Number of fast values that are meaned into one value per window
    Sp_Int = 2 ;
    for j = 1:Sp_Int
        Block_Int = NFavg/Sp_Int*( j-1 ) + 1 : NFavg/Sp_Int*j ;
        fluxW.wt2( CW ) = fluxW.wt2( CW ) + ...
            ( nanmean( w( Block_Int ) .*Tr ( ( Block_Int ))) - ...
            nanmean( w( Block_Int )).*nanmean( Tr( ( Block_Int ))))/Sp_Int;
        fluxW.wq2( CW ) = fluxW.wq2( CW ) + ...
            ( nanmean( w( Block_Int ) .* qWPL( Block_Int )) -...
            nanmean( w( Block_Int )) .* nanmean( qWPL( Block_Int )))/Sp_Int;
        fluxW.wc2( CW ) = fluxW.wc2( CW ) + ...
            ( nanmean( w( Block_Int ) .* cWPL( Block_Int )) -...
            nanmean( w( Block_Int )) .* nanmean( cWPL( Block_Int )))/Sp_Int;
    end
    
    % Storage terms - dx is averaged on a sub-interval (2D [sec]) at the beginning and the
    % end of the interval window
    D=120/1/freqF;
    if CW == 1          % In the case of the first window
        fluxW.dtdt( CW ) = ...
            ( nanmean( DespikeW.CSAT_Tmp( CFV(end)-D : CFV(end)+D )) - nanmean( DespikeW.CSAT_Tmp( CFV(1) : CFV(1)+D )))/AVGtime;
        fluxW.dqdt( CW ) = ...
            ( nanmean( DespikeW.LI7500_Q( CFV(end)-D : CFV(end)+D )) - nanmean( DespikeW.LI7500_Q( CFV(1) : CFV(1)+D )))/AVGtime;
        fluxW.dcdt( CW ) = ...
            ( nanmean( DespikeW.LI7500_C( CFV(end)-D : CFV(end)+D )) - nanmean( DespikeW.LI7500_C( CFV(1) : CFV(1)+D )))/AVGtime;
    elseif  CW == nw    % In the case of the last window
        fluxW.dtdt( CW ) = ...
            ( nanmean( DespikeW.CSAT_Tmp( CFV(end)-D : CFV(end) )) - nanmean( DespikeW.CSAT_Tmp( CFV(1)-D : CFV(1)+D )))/AVGtime;
        fluxW.dqdt( CW ) = ...
            ( nanmean( DespikeW.LI7500_Q( CFV(end)-D : CFV(end) )) - nanmean( DespikeW.LI7500_Q( CFV(1)-D : CFV(1)+D )))/AVGtime;
        fluxW.dcdt( CW ) = ...
            ( nanmean( DespikeW.LI7500_C( CFV(end)-D : CFV(end) )) - nanmean( DespikeW.LI7500_C( CFV(1)-D : CFV(1)+D )))/AVGtime;
    else
        fluxW.dtdt( CW ) = ...
            ( nanmean( DespikeW.CSAT_Tmp( CFV(end)-D : CFV(end)+D )) - nanmean( DespikeW.CSAT_Tmp( CFV(1)-D : CFV(1)+D )))/AVGtime;
        fluxW.dqdt( CW ) = ...
            ( nanmean( DespikeW.LI7500_Q( CFV(end)-D : CFV(end)+D )) - nanmean( DespikeW.LI7500_Q( CFV(1)-D : CFV(1)+D )))/AVGtime;
        fluxW.dcdt( CW ) = ...
            ( nanmean( DespikeW.LI7500_C( CFV(end)-D : CFV(end)+D )) - nanmean( DespikeW.LI7500_C( CFV(1)-D : CFV(1)+D )))/AVGtime;
    end
else
    %sprintf ('Not enough good CO2 or H2O or wind data in window number %d ', CW )
    C = nan( NFavg , 1 ) ;
    Q = nan( NFavg , 1 ) ;
end % if Good_Data_CO2_H2O

Good_Data_CH4 = find( ~isnan( DespikeW.LI7700_M( CFV )) & ~isnan( w ));

if length( Good_Data_CH4 ) > 0.5 * NFavg
    % find lag time between LI7700 and Sonic
    M=DespikeW.LI7700_M(CFV);
    Msc=(M-nanmean(M))/nanstd(M);
    Msc(isnan(M))=0;
    
    fluxW.lagm( CW ) = lag_test( wsc, Msc, LI7700, Sdist_LI7700, ubar, freqF);   % Calculate lag between CH4 and Wind
    
    M = nan( NFavg , 1 ) ;
    % Calculate indices for vector shift CH4
    [V_lag_1, V_lag_end , V_1 , V_end] = Lag_Shift_index (fluxW.lagm(CW) , concatenating_files , CW , NFavg , nw );
    % nw = number of windows per file [number/file]
    % ntsf =  # of timesteps per file - fast data
    M( V_lag_1 : V_lag_end ) = DespikeW.LI7700_M( V_1 : V_end );   % Sonic & LI7700
    
    if length( Good_Data_CO2_H20 ) > 0.5 * NFavg
        % Aapply WPL correction on raw data. Error is caused by changes in air volume deu to heat/water vapor/pressure.
        [~, ~, ~, mWPL, mA, mB, mC] = WPL_correction_10Hz...
            ( Q, C, fluxW.tair(CW), Tson, DespikeW.LI7500_P(CFV), nan, LI7500, DespikeW.LI7700_P(CFV), ...
            M, LI7700, DespikeW.LI7700_P(CFV));
        
        fluxW.mA( CW )= nanmean(mA);
        fluxW.mB( CW )= nanmean(mB);
        fluxW.mC( CW )= nanmean(mC);
        fluxW.mbar(CW) = nanmean(mWPL);    % Sonic & LI7500 & LI7700 CH4 concentration mean
        fluxW.M(CW) = nanmean(M);
        fluxW.mm(CW) = nanvar( mWPL ); % Sonic & LI7500 & LI7700 CH4 concentration variance
        
        fluxW.wm(CW) = nanmean( w .* ( mWPL - fluxW.mbar(CW) ))/fluxW.Spectral_correctionLI7700(CW); % Sonic & LI7500 & LI7700 CH4 Flux in [ Kg/m2/s] Units.
        fluxW.FCH4(CW) = fluxW.wm(CW)/(16.0425/1000)*1000000; % [umol/m2/s]
        % mol CH4 = 16.0425 gr CH4     H=1.00794 [gr/mol] C=12.0107 [gr/mol]
        
        fluxW.m3(CW) = nanmoment(mWPL,3); % = moment( mWPL, 3 ); % Sonic & LI7500 & LI7700 CH4 concentration moment
        
        % Compute fluxes in specific interval (Sp_Int) average for stationary test
        Sp_Int = 2 ;
        for j = 1:Sp_Int
            Block_Int = NFavg/Sp_Int*(j-1)+1:NFavg/Sp_Int*j;
            fluxW.wm2(CW) = fluxW.wm2(CW) + ...
                ( nanmean( w( Block_Int ).*mWPL( Block_Int ))-nanmean( w( Block_Int )).*nanmean( mWPL( Block_Int )))/Sp_Int;
        end
    end
    % Storage terms - dx is averaged on a sub-interval (2D [sec]) at the beginning and the
    % end of the interval window
    D=120/1/freqF;
    if ( CW > 1 ) && ( CW < nw )
        fluxW.dmdt( CW )= ( nanmean( DespikeW.LI7700_M( CFV( end ) - D:CFV( end ) + D )) - ...
            nanmean( DespikeW.LI7700_M( CFV( 1 ) - D:CFV( 1 ) + D )))/AVGtime;
        
    elseif CW == 1
        fluxW.dmdt( CW ) = ( nanmean( DespikeW.LI7700_M( CFV( end ) - D:CFV( end ) + D )) - ...
            nanmean( DespikeW.LI7700_M( CFV( 1 ):CFV( 1 ) + D )))/AVGtime;
        
    elseif CW == nw
        fluxW.dmdt( CW ) = ( nanmean( DespikeW.LI7700_M( CFV( end ) - D:CFV( end ))) - ...
            nanmean( DespikeW.LI7700_M( CFV( 1 ) - D:CFV( 1 ) + D )))/AVGtime;
    end
else
    % sprintf ('Not enough good CH4 data in window number %d ', CW )
end % if Good_Data_CH4

end
