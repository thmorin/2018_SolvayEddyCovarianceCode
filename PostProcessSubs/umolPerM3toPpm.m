function [ cOut ] = umolPerM3toPpm( cIn, T, P )
%function [ cOut ] = ppm2umolPerM3( cIn )
% Written by Tim Morin 07/2014
%
% Converts gas concentrations given in [ppm] to [umol/m^3] mixing ratio 
% using ideal gas law
%
% cIn  = nX1 vector of gas concentrations in units                         [umol/m^3]
% T    = nX1 vector of temperatures                                        [C]
% P    = nX1 vector of pressures                                           [Pa]
% cOut = nX1 vector of gas mixing ratios in units                          [ppm]

%%
%%-------------------------Constants---------------------------------------
R=8.314;                                       %Universal gas constant     [m^3 Pa K^-1 mol^-1]
T_k = T+273;                                   %Unit conversion            [K]

%%
%%-------------------------Convert concentration---------------------------
cOut = cIn*R.*T_k./P;                            %Concentration conversion   [umol/m^3]
