function [spikeflag, ustk, ustfilt]=despikeandfiltercor3(c_flux,u_star,daynight,t_ave)
%This program first runs a despiking algorithm after Papale et al. (2006)
%It then runs an algorithm to determine ustar thresholds modified from Reichstein et al (2005)

%Updated May 2010 by K. Maurer and G. Bohrer

%Inputs
%Fc - Carbon flux...the algorithm should not be sensitive to units for NEE
%ustar - friction velocity m/s
%daynight - a variable that flags daytime fluxes as 1, nighttime fluxes as 0
%Ta - air temperature  degrees C
%year - year
%DOY - day of year

%Outputs
%spikeflag - the output from the despiking routine.  1 = good data, 0 = bad data
%ustk - acceptable ustar values for the roughly 2-month periods bounded by cdoy > mindoy & cdoy < maxdoy 
%mindoy, maxdoy - the days bounding the roughly 2-month perios for which 

Fc=c_flux;
Ta=t_ave;
ustar=u_star;
spikeflag=zeros(size(Fc));

%Despiking algorithm ... determines spikes separately for daytime and nighttime data

wndw=1:672:length(Fc); 
cc=length(wndw); 
wndw(cc)=length(Fc); %672 is equal to two weeks of half hourly averages

for i=1:(length(wndw)-1)
    mask=wndw(i):1:wndw(i+1);
    Fci=Fc(mask);
    
    Fcdayi=Fci(daynight(mask)==1);
    Fcnighti=Fci(daynight(mask)==0);
%    sdayi=ones(size(Fcdayi));
%    snighti=ones(size(Fcnighti));
    
    medday=nanmedian(Fcdayi); 
    mednight=nanmedian(Fcnighti);
    threshday=nanmedian(abs(Fcdayi-medday));
    threshnight=nanmedian(abs(Fcnighti-mednight));
    % updated to z value of 5.5 as of Papale et al 2006, KS 10/20/09
    minday=medday-5.5*threshday/0.6745;
    maxday=medday+5.5*threshday/0.6745;
    minnight=mednight-5.5*threshnight/0.6745;
    maxnight=mednight+5.5*threshnight/0.6745;
    
    spikeflag(mask)= (((Fci>minday | Fci<maxday) & daynight(mask)==1) | ((Fci>minnight | Fci<maxnight) & daynight(mask)==0)); 
    
end

%Run Reichstein Ustar routine
[ustk, ustfilt]=GPF3(Fc,spikeflag,Ta,ustar,daynight); 

        
    
            
        

    