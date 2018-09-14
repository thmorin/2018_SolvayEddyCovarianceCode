function [rat]=opirgafreq(u,lso,z_geom,TB,zbyl,dirga,lirga)
% new function which will correct opened path IRGA
% attenuation using massmans analytical methodology
% inputs required are 
% u - average horizontal wind speed [m/s]...
% lso - line averaging for sonic -length [m]
% z_geom = geometric sampling height [m] above ground
% TB-block averaging period in seconds
% zbyl  - z/L Monin -Obukhov length scale
% usage [rat]=opirgafreq(u,lso,z_geom,TB,zbyl);
% no lateral or vertical seperation is specified as the sensors are alinged
% carefully
% no linear detrending or high pass filtering 


if (zbyl>0 && zbyl<=2); % stable atmospheric conditions closed path-slow response (0.1-0.3)
    %nx = 2-(1.95/zbyl); % Hollinger (2004), 10-1689-1706 Global change biology
    nx=2-1.915/(1+0.5*zbyl); % from Horst BLM(1997) 82: 219-233
elseif (zbyl<= 0);
    nx=0.085;
else
    nx=NaN;
end
fx=(u*nx)/z_geom;

%lirga is the line averaging in the sensor-for LI7500 it is 0.15m
 % time constants  
 t1=lirga/(4*u); % line averaging for open path  (Horst 1997 & 2000)
 t2=lso/(8.4*u); % scalar flux line averaging
 t3=(0.2+0.4*dirga/lirga)*(lirga/u);%volume averaging-see massman 2000
 t4=0.016;% open path band width (Horst 1997, Boundary layer meteorology , 82:219-233)
 %t5=llat/(1.1*u); %for lateral seperation llat=lateral separation in (m)
 

 te=sqrt((t1^2)+(t2^2)+(t3^2)+(t4^2));
 tb=TB/2.8; %time constant for block averaging  
 %presenting the  b and p terms

 b=2*pi*fx*tb;
 p=2*pi*fx*te; 
 alf=0.925; % shape of the cospectra factor
 if zbyl<=2 && b>=0 && p>=0
    rat=((b^alf)/((b^alf)+1))*((b^alf)/((b^alf)+(p^alf)))*(1/((p^alf)+1)); % For all conditions
 else
    rat=1;
 end
 if isnan(rat)
     rat=1;
 end;
