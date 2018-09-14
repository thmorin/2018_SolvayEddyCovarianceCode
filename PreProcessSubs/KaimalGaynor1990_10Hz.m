function [Tr,e2] = KaimalGaynor1990_10Hz(P,Q,T,IRGA)

% P from pressure sensor [Pa]. If no pressure sensor, use pressure from open path IRGA
% Q for Open path [mmol/m^3]. For Closed path [mmol/mol] for KH20 [g/m3] RAW DATA
% T from SONIC [Deg C]
% IRGA      1 for OPEN PATH (LI7500 & LI7700)       2 or 3 for CLOSED PATH (LI6262 or LI7000)       4 for KH20 Hygrometer   

TK = T+273.16;                              % temp in [Deg K]
R=8.314;                                    % Universal gas constant [J/K/mol]
    
P(isnan(P))=101000;

if IRGA==1                                  % Open path [mmol/m3]
    Q=Q/1000;                               % from [mmol/m^3] to [mol/m^3]
    rho=P./(R.*TK);
    
elseif IRGA==2 || IRGA==3                   % Closed path [mmol/mol]
    rho=P./(R.*TK);                         % air density in [mol/m^3]
    Q=Q/1000.*rho;                          % from [mmol/mol] to [mol/m^3]
        
elseif IRGA==4                              % KH20
    Q=Q/18;                                 % from g/m3 to [mol/m^3]
    rho=P./(R.*TK);
end

e2 = R*Q.*TK;                               % Q must be in [mol/m^3]

%SH = Q./(Q+1);                              % Specific humidity in decimal form
%Hum = SH./0.622.*((rho-Q)./Q)./100;         % RH in decimal form

%Pvaporsat = 611.2*exp(17.67*T./(243.51+T)); % Saturated Vapor Pressure
%Pvapor = Hum.*Pvaporsat;                    % Vapor Pressure

%e2corr = Pvapor/nanmean(e2).*e2; 

%Tr = TK./(1 + 0.32.*e2corr./P)-273.16;      % Real Temperature [Deg C] Kaimal and Gaynor 1990 Eq.3
Tr = TK./(1 + 0.32.*e2./P)-273.16; 
end