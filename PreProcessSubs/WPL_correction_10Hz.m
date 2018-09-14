function  [rho,qWPL,cWPL,mWPL,A,B,C] = WPL_correction_10Hz(Q,CO2,T,Tr,IRGA_P,IRGA_T,IRGA_C,Pbar_in,M,IRGA_M,Mp)

% WPL correction - After Webb-Pearman-Leuning (Webb et al.,1980) is used so
% the gas flux itself that is calculated from the changes in the gas 
% density is not atributed to the changes in density caused by changes in 
% air temprature, air humidity(water vapor dilution), and air pressure. 
% these air density changes can lead to false readings of gas (H2O/CO2/CH4)
% flux.
%
% Q - Water vapor fluxuations. For Open path [mmol/m^3] for Closed path
% [mmol/mol] for KH20 (KH20 is an Hygrometer that only measures water
% vapor) [g/m^3] RAW DATA.
%
% C - CO2 fluxuations. For Open path [mmol/m^3] for Closed path [umol/mol]
% RAW DATA.
%
% T - Temperature. [Deg C] from SONIC anamometer.
%
% Tr - Rreal temperature. [Deg C] from Kaimal & Gaynor process
%
% IRGA_P - IRGA sensor pressure, [Pa] RAW DATA
%
% IRGA_T - IRGA sensor temperature,[Deg C] Only with Closed path, make NAN
% for Open path. RAW DATA.
%
% IRGA_C - sensor type code.  1 for OPEN PATH (LI7500 or LI7700)   2 or 3 for CLOSED
% PATH (LI6262 or LI7000)   4 for KH20 Hygrometer.
%
% pbar - air pressure from another sensor, [Pa] This will be used if
% there isn't a measurment from a faster sensor, make NaN if unavailable.
% 
% M - Methane fluxuations. [mmol/m^3] from LI 7700 , do NOT input if no
% methane sensor exists, also can be used for any arbitrary scalar flux,
% CHECK UNITS!!!.
% 
% IRGA_M [methane sensor handling code]     1 for OPEN PATH (LI7700)
% 2 for CLOSED PATH WITH simultaneous measurement of water vapor inside
% chamber      3 for CLOSED PATH WITHOUT simultaneous measurement of water
% vapor in the chamber.
% 
% Mp -methane sensor pressure, [Pa] is pressure from methane sensor
%
% Mt -methane sensor temperature, [Deg C] is temperature from closed path
% methane sensor, make NaN if methane is obtained with OPEN path sensor.
%
% rho - Air density (for open path IRGA) / Chamber's air density (for closed path IRGA) [mol/m^3].
% qWPL - WPL corrected water vapor pertubations (minus the bar value) [kg/m^3]
% cWPL - WPL corrected CO2 flux pertubations (minus the bar value) [kg/m^3]
% mWPL - WPL corrected Methane flux [kg/m^3]

 
Pbar = nanmean(Pbar_in);

    if IRGA_C ~=1 && IRGA_C ~=2 && IRGA_C ~=3 && IRGA_C ~=4 
        disp('BAD IRGA INPUT')
        return
    end
    
%%% Constant
    ma=29.002/1000;                             % Molar mass of dry air [Kg/mol]
    mv=18/1000;                                 % Molar mass of water vapor [Kg/mol]
    mc=44/1000;                                 % Molar mass of co2 [Kg/mol]
    mm=16.042/1000;                             % Molar mass of methane [Kg/mol]
    R=8.314510;                                 % Universal gas constant [J/K/mol or m^3*Pa/K/mol]

    TK=T+273.16;                                % Air temperature [Deg K] from temp probe
    TKbar=nanmean(TK);                          % Mean of Air temp [Deg K]
    TrK=Tr+273.15;                              % Real air temperature [Deg K] from Kaimal & Gaynor process
    TrKbar=nanmean(TrK);                        % Mean of real air temperature [Deg K]

    if ~isnan(IRGA_T)
        TIbar=nanmean(IRGA_T)+273.16;           % IRGA_C temperature [Deg K]
    else
        TIbar=nanmean(TK);
    end

    if sum(~isnan(IRGA_P))>=0.6*sum(~isnan(Q))
        PIbar=nanmean(IRGA_P);
    elseif sum(~isnan(Mp))>=0.6*sum(~isnan(Q))
        PIbar=nanmean(Mp);
    elseif ~isnan(Pbar)
        PIbar=Pbar;
    else
        PIbar=101000;
    end

    if IRGA_C == 2 || IRGA_C==3                     % Convert from mixing ratio to molar density for closed path (LI6262, LI7000).
        rho=PIbar/(R*TIbar);                    % Chamber air density  [mol/m^3]
        qKG=Q*rho*mv/1000;                      % [Kg/m^3]
        cKG=CO2*rho*mc/1000000;                   % [Kg/m^3]
        rhov=nanmean(Q)*mv/1000*rho;            % Vapor density [Kg/m^3]
        rhoc=nanmean(CO2)*mc/1000000*rho;         % CO2 density [Kg/m^3]
        rhoa=(rho-rhov/mv-rhoc/mc)*ma;          % Dry air density  [Kg/m^3]        
    elseif IRGA_C ==1                             % Open Path with wc & wq units mmol/m^2/s
        rho=PIbar/(R*TKbar);                    % Air density of LI7500 [mol/m^3]
        qKG=Q*mv/1000;                          % [Kg/m^3]
        cKG=CO2*mc/1000;                          % [Kg/m^3]  
        rhov=nanmean(Q)*mv/1000;                % Vapor density [Kg/m^3]
        rhoc=nanmean(CO2)*mc/1000;                % CO2 density [Kg/m^3] 
        rhoa=(rho-rhov/mv-rhoc/mc)*ma;          % Dry air density  [Kg/m^3]
    elseif IRGA_C == 4                            % KH20 Hygrometer
        rho=Pbar/(R*nanmean(TrK));              % Air density  [mol/m^3]
        qKG=Q*1000;                             % [Kg/m3]
        rhov=nanmean(Q)*1000;                   % Vapor density [Kg/m^3]
        rhoc=nan;
        rhoa=(rho-rhov/mv)*ma;                  % Dry air density  [Kg/m3^]
    end 


    %%% Detto and Katul BLM, 2006
    mu=ma/mv;
    si=rhov/rhoa;

    if IRGA_C == 1 % Open Path Correction

        qWPL=qKG + (qKG-nanmean(qKG))*mu*si + nanmean(qKG)*(1+mu*si)*(TrK-TrKbar)/TKbar;  % [kg/m^3]
        cWPL= cKG + mu*rhoc/rhoa*(qKG-nanmean(qKG)) + rhoc*(1+mu*si)*(TrK-TrKbar)/TKbar;  % [kg/m^3]
        
    elseif IRGA_C == 2 || IRGA_C == 3 % Closed Path Correction 

        qWPL = Pbar.*TIbar./(TK.*PIbar).*(1+mu*si).*(qKG-nanmean(qKG)) + nanmean(qKG);                              % [kg/m^3]
        cWPL = Pbar.*TIbar./(TK.*PIbar).*((cKG-nanmean(cKG))+ mu*rhoc/rhoa.*(qKG-nanmean(qKG))) + nanmean(cKG);     % [kg/m^3]

    elseif IRGA_C == 4 % KH20 Hygrometer

        qWPL=qKG + (qKG-nanmean(qKG))*mu*si + nanmean(qKG)*(1+mu*si)*(TrK-TrKbar)/TKbar;  % [kg/m^3]
        cWPL = nan;

    end

    
   A=nan;
   B=nan;
   C=nan;
   
    
    if nargin >= 9
        if IRGA_M == 1 % Open Path Correction
            
            if sum(~isnan(Mp))>=0.6*sum(~isnan(M))
                P=Mp;
            elseif (sum(~isnan(IRGA_P))>=0.6*sum(~isnan(M)))
                P=IRGA_P;
            elseif ~isnan(Pbar) 
                P=Pbar;
            else
                P=101000;
            end

            %rhoMeth = M; %PMbar/(R*TrKbar);        % M = Air density of LI7700 [mmol/m^3]
            mKg = M*mm/1000;                        % [mmol CH4] = 0.001 [mol CH4] = 0.016042 [Kg]
            rhom=nanmean(mKg);                      % Mean of methane density [Kg/m^3] 
            
            [A, B, C] = Calc_SpectComp( Tr, P/1000, (nanmean(qWPL)/mv*1000)*R*TKbar/(1000*nanmean(P)) );

            mWPL= A.*( mKg + B .*( mu*rhom/rhoa*(qWPL-nanmean(qWPL))) + C.*( rhom*(1+mu*si)*(TrK-TrKbar)/TKbar));  % [kg/m^3]
            % A B and C are corrections to counter affect the of the spectrogoraphy 

        end
    else mWPL = nan;    
    end
end
