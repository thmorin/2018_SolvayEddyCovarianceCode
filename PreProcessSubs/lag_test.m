function lag = lag_test( w, v, IRGA, dist, windspeed, freqF  )

%lag_test.m
%
%   lag_test, find the time lag between two vectors using the maximum
%   correlation method
%   w and v are the two vectors,
%   w = wind from sonic,
%   v = IRGA concentration measurement (CO2/H2O/CH4).
%   IRGA = type determains the span. The lag can only be in the [-span span] range
%   a positive lag means that the signal c is delayed respect w
%   dist = Distance between IRGA and Sonic [m]
%   windspeed = mean wind speed at Current Window
%  Aug 1: Added condition regarding mean wind speed and separation distance
% (distance between IRGA and Sonic) for IRGA=1 (Open path).
if IRGA==1      % Open Path
    if nargin==6 && windspeed > 0.2
        span = floor(freqF*2*dist/windspeed)+2;
        lag0 = 0;
        low = -span;
%         disp(['span is:' num2str(span) 'and WS is: ' num2str(windspeed)...
%             'dist is:' num2str(dist) ]);
    elseif nargin==6  
        span = 5;
        lag0 = 0;
        low = -span;
    else
        span=0;
        lag0=0;
        low=0;
        disp('nargin for lag shift is not 6');
    end
elseif IRGA==2  % Closed Path (AF 46m & FASET 32m)
    span = 120;
    lag0 = 60;
    low = 30;
elseif IRGA==3  % Closed Path (AF 34m)
    span = 100;
    lag0 = 45;
    low = 25;
elseif IRGA==4  % KH20 Hygrometer
    span = 50;
    lag0 = 5;
    low = -span;
end

[ wv_xcov, lags ] = xcov( w, v, span);

% acceptable = -lags > low ;

[ ~ , lag_max_index] = max( abs( wv_xcov( -lags > low )));

%  n=1;
%  while lags<low
%     lag=-lags(lag_max(n));
%     n=n+1;
%  end

lag = -lags( lag_max_index );

if abs( lag ) > span
    lag = lag0 ;
end
end