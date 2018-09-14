% Function to calculate the sensible heat flux coefficient



function [rho_cp]   =   calculate_rho_cp(t,p,e,E)

% t : air temperatue in [deg Celcius]
% p : pressure in [Pa]
% e : partial vapur pressure [Pa]
% E : saturation vapor pressure [Pa]

R_d     =   8.31451; % Unit [Nm mol^{-1} K^{-1}] according to Foken (2005), Springer

% density of dry air
rho_a   =   29.002*p./(R_d*(t+273.15));  % [g m-3]

% density of water vapour
rho_v   =   18.01*e./(R_d*(t+273.15));              % [g m-3]

% Relative humidity
rH      =   e./E *100;                                   % [1]

% specific heat for dry air
cp_dry  =   1005 + ((t+23.12).^2./3364);                 % [J kg-1 K-1]

% specific heat of moist air
cp_moist=   1859 + 0.13*rH + t*(0.193+0.00569*rH)+t.^2*(0.001+0.00005*rH); % [J kg-1 K-1]

% rho times cp
rho_cp  =   cp_moist*rho_v./1000 + cp_dry*rho_a./1000;

end