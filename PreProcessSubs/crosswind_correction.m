% Function correcting the sensible heat flux for the effect of water
% vapour, i.e. converting buoyancy flux into sensible heat flux, and
% performing the cross wind correction
% Method used after Liu et al. 2001, BLM

% Developed and written by Christoph Thomas, 
% Dept. of Forest Science, Oregon State University, 2006
% last update 20-Nov-2006


function [wT_corr]   =  crosswind_correction(wT,T,p,a_mean,sonic,cov_m,stat_uvw)

%   wT          uncorrected covariance of sonic temperature and vertical wind [K m s-1]
%   T           mean dry bulb temperature [deg C]
%   p           barometric pressure [Pa]
%   cov_m       covariance vector for shear stress components, contains u'w' and v'w' [m2 s-1]
%   stat_uvw    matrix containing statistics of wind components u,v, and w (all 4 moments)
%   stat_s      matrix containing statistics of scalars components sonic Temp, and other scalars if present (all 4 moments)
%   sonic       type of sonic anemometer, used to determine sensor specific coefficients; '1' Gill R2, '2' Metek USAT, '3' Gill R3-50, '4' Young 81000, '5' Campbell CSAT3
%   wa          covariance of absolute humidity (a) and vertical wind, i.e. vertical latent heat flux [mmol m-2 s-1]
%   a_mean      mean absolute humidity [mmol m-3]



    % Extract needed variables from input matrices
    u_mean          =   stat_uvw(1);
    v_mean          =   stat_uvw(2);
    uw              =   cov_m(1);
    vw              =   cov_m(2);

    % Calculating mean specific humidity
    a_mean          =   a_mean / 1000 * 0.01802;                    % conversion into [kg m-3]
    e_mean          =   a_mean * (T+273.15) / 0.21667 /10;          % partial pressure of water vapour [kPa]
    q_mean          =   0.62198 * e_mean / (p - 0.378*e_mean);      % specific humidity [kg kg-1]

    % Calculate speed of sound (squared) from sonic Ts_mean
    % c2              =   (Ts_mean + 273.15)*403;
    c2              =   403.36 * (T+273.15)*(1 + 0.51 * q_mean);    % [m2 s-1]

    % Determine coeffs A and B according to sonic anemometer for cross-wind correction 
    if     sonic == 1
        A           =   0.5;
        B           =   1;
    elseif sonic == 2
        A           =   3/4;
        B           =   3/4;
    elseif sonic == 3
        A           =   1 - 0.5*cos(deg2rad(45)).^2;
        B           =   1 - 0.5*cos(deg2rad(45)).^2;
    elseif sonic == 4 %Better use scotaneous here instead of cross wind
        A           =   3/4;            % wild guess
        B           =   3/4;            % wild guess
    elseif sonic == 5  %CSAT does not need crosswind corrcection
        A           =   7/8;
        B           =   7/8;
    elseif sonic == 6
        A           =   1 - 0.5*cos(deg2rad(45)).^2;
        B           =   1 - 0.5*cos(deg2rad(45)).^2;
    end;

    % Converting buoynacy into sensible heat flux and applying cross-wind corrrection
    wT_corr         =   wT + (2*(T+273.15)/c2*(u_mean*uw*A + v_mean*vw*B));  % [K m s-1]


end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function arad = deg2rad(adeg)

arad = adeg/180*pi;
end