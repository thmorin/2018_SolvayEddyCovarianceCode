% Function correcting the sensible heat flux for the effect of water
% vapour, i.e. converting buoyancy flux into sensible heat flux, and
% performing the cross wind correction
% Method used after Liu et al. 2001, BLM

% Developed and written by Christoph Thomas, 
% Dept. of Forest Science, Oregon State University, 2006
% last update 20-Nov-2006


function [wT_corr]   =  calculate_wT(wT,T,p,wa)

%   wT          uncorrected covariance of sonic temperature and vertical wind [K m s-1]
%   T           mean dry bulb temperature [deg C]
%   p           barometric pressure [Pa]
%   cov_m       covariance vector for shear stress components, contains u'w' and v'w' [m2 s-1]
%   stat_uvw    matrix containing statistics of wind components u,v, and w (all 4 moments)
%   stat_s      matrix containing statistics of scalars components sonic Temp, and other scalars if present (all 4 moments)
%   sonic       type of sonic anemometer, used to determine sensor specific coefficients; '1' Gill R2, '2' Metek USAT, '3' Gill R3-50, '4' Young 81000, '5' Campbell CSAT3
%   wa          covariance of absolute humidity (a) and vertical wind, i.e. vertical latent heat flux [mmol m-2 s-1]
%   a_mean      mean absolute humidity [mmol m-3]

% Converting the vertical latent heat flux into vertical specific latent heat flux
wa              =   wa / 1000 * 0.01802;                        % [kg m-2 s-1]
we              =   wa * (T+273.15) / 0.21667 / 10;             % [kPa m-2 s-1]
wq              =   0.62198 * we / (p/1000);                           % [kg m kg-1 s-1]
wT_corr         =   wT - (0.51*wq*(T+273.15)); % no crosswind correction for CSAT3



end

