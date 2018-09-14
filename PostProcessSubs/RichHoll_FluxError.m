function [MeanError,StdError] = RichHoll_FluxError(flux,par,temp,ubar)
% Program calculates the random error of the CO2 flux for a given year
% Inputs:
% flux: flux that you want seasonal error calculated for, filtered for u*
% par: PAR, or equivalent radiation variable
% temp: temperature
% ubar: average wind speed
% season: vector to specify seasons
%
% Outputs:
% MeanError: mean error based on similar time periods difference in flux
% StdError: standard deviation of error based on similar time periods difference in flux
% Written by: KD Maurer, July 2011
% Modified by: T Morin, Oct 2012
 
F=flux;
P=par;
T=temp;
U=ubar;

MeanError=nan(size(flux));
StdError=nan(size(flux));

Delta=[];

for i=1:length(F)-48
    if abs(P(i)-P(i+48))<=75 && abs(T(i)-T(i+48))<=3 && abs(U(i)-U(i+48)<=1) && ~isnan(F(i)) && ~isnan(F(i+48))
        Delta=[Delta;abs((F(i)-F(i+48))/sqrt(2))];
    end
end

MeanError=nanmean(Delta);
StdError=nanstd(Delta);

end