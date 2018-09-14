function [ cOut ] = PpmToumolPerM3( cIn, T, P )
%function [ cOut ] = PpmToumolPerM3( cIn )
% Written by Tim Morin 07/2014
%
% Converts gas concentrations given in [ppm] to [umol m^-3] mixing ratio 
% using ideal gas law
%
% cIn  = nX1 vector of gas concentrations in units                         [ppm]
% T    = nX1 vector of temperatures                                        [C]
% P    = nX1 vector of pressures                                           [Pa]
% cOut = nX1 vector of gas mixing ratios in units                          [umol m^-2 s^-1]

%%
%%-------------------------Constants---------------------------------------
R=8.314;                                       %Universal gas constant     [m^3 Pa K^-1 mol^-1]
T_k = T+273;                                   %Unit conversion            [K]

%%
%%-------------------------Convert concentration---------------------------
cOut = bsxfun(@times, bsxfun(@rdivide, cIn/R ,T_k) , P); %Conc. conversion [umol/m^3]