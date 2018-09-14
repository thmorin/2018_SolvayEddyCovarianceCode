function [DOY,DECDAY,DECYear,C,Fc,H,L,LE,LWinbar,LWoutbar,net_lw,net_sw,pressure,...
    p_vapor_bar,p_vaporSat_bar,Q,Resp,Resp_Error,Resp_Reconstructed,rho_cp,rH,RNET,SoilTemp,Spectral_correction,SWinbar,SWoutbar,tair,ustar_1,...
    wq,wts,ww,vv_1,wind_dir_1,wind_speed_1,ustar_2,vv_2,wind_dir_2,wind_speed_2,Total_PAR,Diffuse_PAR,C1,C2]=...
    InitializeArrays(totallength,loops,nargin)

DOY=nan(totallength*48,1); %[days]
DECDAY=nan(totallength*48,1); %[days]
DECYear=nan(totallength*48,1); %[years]

Total_PAR=nan(totallength*48,1);
Diffuse_PAR=nan(totallength*48,1);


C=nan(totallength*48,1); %[mmol/m^3]
C1=nan(totallength*48,1); %[mmol/m^3]
C2=nan(totallength*48,1); %[mmol/m^3]

Fc=nan(totallength*48,1); %[umol/m^2/s]

H=nan(totallength*48,1); %[W/m^2]

L=nan(totallength*48,1); %[m]

LE=nan(totallength*48,1); %[W/m^2]

LWinbar=nan(totallength*48,1); %[W/m^2]
LWoutbar=nan(totallength*48,1); %[W/m^2]

net_lw=nan(totallength*48,1); %[W/m^2]

net_sw=nan(totallength*48,1); %[W/m^2]

pressure=nan(totallength*48,1); %[Pa]

p_vapor_bar=nan(totallength*48,1); %[Pa]

p_vaporSat_bar=nan(totallength*48,1); %[Pa]

Q=nan(totallength*48,1); %[mmol/m^3]

Resp=nan(totallength*48,loops); %[umol/m^2/s]
Resp_Error=nan(totallength*48,loops); %[umol/m^2/s]
Resp_Reconstructed=nan(totallength*48,loops); %[umol/m^2/s]

rho_cp=nan(totallength*48,1); 

rH=nan(totallength*48,1); %[100]

RNET=nan(totallength*48,1); %[W/m^2]

SoilTemp=nan(totallength*48,1); %[C]

Spectral_correction=nan(totallength*48,1);

SWinbar=nan(totallength*48,1); %[W/m^2]
SWoutbar=nan(totallength*48,1); %[W/m^2]

tair=nan(totallength*48,1); %[C]

ustar_1=nan(totallength*48,1); 

ustar_2=nan(totallength*48,1); 

vv_1=nan(totallength*48,1);

vv_2=nan(totallength*48,1);

wind_dir_1=nan(totallength*48,1); %[rad]
wind_dir_2=nan(totallength*48,1); %Converted to [deg] in wetland_build_arrays.m

wind_speed_1=nan(totallength*48,1);
wind_speed_2=nan(totallength*48,1);

wq=nan(totallength*48,1); %[kg/m2/s]

wts=nan(totallength*48,1);

ww=nan(totallength*48,1);
end