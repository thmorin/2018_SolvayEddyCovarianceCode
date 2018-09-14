 function [rat]=cpirgafreq(u,lso,zbyl,TB,tair,lpm,a,sflg,irflg,Lt,z_geom)
% new function which will correct closed path IRGA
% attenuation using massmans analytical methodology for sonic line
% averaging, IRGA line averaging, IRGA volume averaging, tube attenuation and
% block averaging
% inputs required are 
% u - average horizontal wind speed (m/sec)...
% fs - data acqusition frequency
% TB-block averaging period in seconds
% lso - line averaging for sonic -length(m)
% zbyl  - z/L Monin -Obukhov length scale
% a - inner diamter of the tube in meters (sampling tube)
% lpm - flow rate through the tube in liters per minute
% tair = air temperature in °C
% sflg - is the flag for scalar 1 for CO2 and 2 for H2O
% Lt length of tubemetric sample height
% z_geom is the geometric sampling height
% irflg=irga flag 1-7000 2-6262
% usage [rat]=cpirgafreq(u,lso,zbyl,tavg,tair,lpm,a,sflg,irflg,Lt,z_geom);

tair = tair + 273.15; % convert air temperature to kelvin

if irflg==1
    Vchamb=10.86*10^-6; %(units = m^3)
    lirg=0.1524; %units of meter
elseif irflg ==2
    Vchamb=1.19*10^-5;  %(units = m^3)
    lirg=0.152;% units in m
else
    Vchamb=1.609*10^-5;  %(units = m^3)
    lirg=0.125;% units in m
end

%using a theoretical defenition of fx-frequency of maximum cospectra
if (zbyl>0 && zbyl<=2); % stable atmospheric conditions closed path-slow response (0.1-0.3)
    %nx = 2-(1.95/zbyl); % Hollinger (2004), 10-1689-1706 Global change biology
    nx=2-1.915/(1+0.5*zbyl); % from Horst BLM(1997) 82: 219-233
elseif (zbyl<= 0);
    nx=0.085;
else
    nx=NaN;
end
fx=(u*nx)/z_geom;

 if sflg == 1
 D=0.1381*10^-4; %m^2 s^-1 molecular diffusivity of CO2
 else
 D=0.2178*10^-4;% m^2 s^-1 molecular diffusivity of H2O
 end

 v=(-1.1555e-14*tair^3) + (9.5728E-11*tair^2) + (3.7604e-08*tair) - 3.4484e-06;
 %(polynomial relationship with air temperature)
 %you can also use-from Raupachs relationship
 %v=((0.0916*(tair-273.15))+13.195)*10^(-6);


 %tube attenuation
 Ut=lpm/(1000*60*pi*(a/2)^2);
 %v=(3.719815-(0.0358737*Pr))+((0.0007342*Pr^2)/100000.0);
 re=(2*(a/2)*Ut)/v; % reynolds number
 %This part adapted from (Leuning and Judd, 1992)
 % a good paper to look will be Takanori shimizu-: Boundary-layer meteorology
 % 2007 122:417:438 -Practical applicability of high frequency correction
 % theories to CO2 flux measured by a closed-path system.
 if re<2300
     del=(0.0104*v*re)/D;
 elseif (re>= 2300 && sflg ==1)
     del=0.693/((8*pi^2)*(log(0.75*(re^0.04)))^2); % for CO2
 elseif (re>=2300  && sflg ==2)
     del=0.693/((8*pi^2)*(log(0.76*(re^0.039)))^2); % for H2O
 else
     del=NaN;
 end
 
 %volume averaging
 tvol=Vchamb/(pi*(a/2)^2*Ut);
 
 %lirga for closed path is 0.1524 m  
 % time constants  
 t1=lirg/(4.0*u); % for scalar means closed path
 t2=lso/(8.4*u); % scalar flux line averaging for sonic
 t3=(sqrt((del*(a/2))/Lt)/0.83)*(Lt/Ut); % tube attenuation
 t4=0.3*tvol;       %volume averaging
 te=sqrt((t1^2)+(t2^2)+(t3^2)+(t4^2));
 tb=TB/2.8;   
 %presenting the  b and p terms
 %a=2*pi*fs*th; is used only if linear detrend or a recursive filter is
 %used for determining the covariance
 b=2*pi*fx*tb;
 p=2*pi*fx*te;
 alf=0.925; % shape of the cospectra factor

 
     
 if (zbyl>0 && zbyl<=2); % stable atmospheric conditions closed path-slow response (0.1-0.3)
     rat=((b^alf)/((b^alf)+1))*((b^alf)/((b^alf)+(p^alf)))*(1/((p^alf)+1))*((1+(0.9*(p^alf)))/(1+(p^alf)));
 elseif (zbyl<= 0); %unstable atmospheric conditions closed path-slow response (0.1-0.3)
     rat=((b^alf)/((b^alf)+1))*((b^alf)/((b^alf)+(p^alf)))*(1/((p^alf)+1));
 else
     rat=1; % no correction applied if conditions are not met
 end
 
 
 if isnan(rat)
     rat=1;
 end
