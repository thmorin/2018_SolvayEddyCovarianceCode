%':::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
%These functions calculate sunrise and sunset times for any
% given latitude and longitude. 
%:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
function [localTrise, localTset]= Suntime(latitude, longitude, localOffset, dst, year, DOYstart) 

if DOYstart<=31
    month=1;
    day=DOYstart;
elseif DOYstart<=59
    month=2;
    day=DOYstart-31;
elseif DOYstart<=90
    month=3;
    day=DOYstart-59;
elseif DOYstart<=120
    month=4;
    day=DOYstart-90;
elseif DOYstart<=151
    month=5;
    day=DOYstart-120;
elseif DOYstart<=181
    month=6;
    day=DOYstart-151;
elseif DOYstart<=212
    month=7;
    day=DOYstart-181;
elseif DOYstart<=243
    month=8;
    day=DOYstart-212;
elseif DOYstart<=273
    month=9;
    day=DOYstart-243;
elseif DOYstart<=304
    month=10;
    day=DOYstart-273;
elseif DOYstart<=334
    month=11;
    day=DOYstart-304;
elseif DOYstart<=366
    month=12;
    day=DOYstart-334;
else 
    disp('Please check if DOYstart falls within 1 through 366');
end

date = [year month day];

% An zenith altitude of 30 minutes i.e. 0.833 radians is generally accepted
% as the angle of the sun at which sunrise/sunset occurs. It is not exactly 
% zero because of refraction effects of the earth's atmosphere. 

zenith = 90+50/60;

%% determine the sunrise and sunset time

% ii. Since 2007, Daylight Saving Time (DST) runs second sunday of March 
% through first sunday of November

for i = 1:7
    dsmatch = strmatch(datestr([date(1) 3 7+i 0 0 0],'ddd'),'Sun');
    dematch = strmatch(datestr([date(1) 11 i 0 0 0],'ddd'),'Sun');
    if dsmatch
        dststart = [date(1) 3 7+i 0 0 0];
        N1 = floor(275* dststart(2)/9);
        N2 = floor((dststart(2)+9)/12);
        N3 = (1 + floor((dststart(1) - 4*floor(dststart(1)/4)+2)/3));
        Nds = N1 - N2*N3 + dststart(3) - 30;
    end
    if dematch
        dstend = [date(1) 11 i 0 0 0];
        N1 = floor(275* dstend(2)/9);
        N2 = floor((dstend(2)+9)/12);
        N3 = (1 + floor((dstend(1) - 4*floor(dstend(1)/4)+2)/3));
        Nde = N1 - N2*N3 + dstend(3) - 30;
    end
end

% 1. calculate the day of the year

N1 = floor(275* date(2)/9);
N2 = floor((date(2)+9)/12);
N3 = (1 + floor((date(1) - 4*floor(date(1)/4)+2)/3));
N = N1 - N2*N3 + date(3) - 30;

% 2. convert the longitude to hour value and calculate an approximate time

lngHour = longitude/15;
trise = N + (6-lngHour)/24;
tset = N + (18-lngHour)/24;

% 3. calculate the Sun's mean anomaly
	
Mrise = 0.9856*trise - 3.289;
Mset = 0.9856*tset - 3.289;

% 4. calculate the Sun's true longitude
	
Lrise = mod(Mrise + 1.916*sind(Mrise) + 0.020*sind(2*Mrise) + 282.634,360);
Lset = mod(Mset + 1.916*sind(Mset) + 0.020*sind(2*Mset) + 282.634,360);

% 5a. calculate the Sun's right ascension
	
RArise = atand(0.91764*tand(Lrise));
RAset = atand(0.91764*tand(Lset));

% 5b. right ascension value needs to be in the same quadrant as L

Lrisequadrant  = floor(Lrise/90)*90;
RArisequadrant = floor(RArise/90)*90;
RArise = RArise + (Lrisequadrant - RArisequadrant);

Lsetquadrant  = floor(Lset/90)*90;
RAsetquadrant = floor(RAset/90)*90;
RAset = RAset + (Lsetquadrant - RAsetquadrant);


% 5c. right ascension value needs to be converted into hours

RArise = RArise/15;
RAset = RAset/15;

% 6. calculate the Sun's declination

sinDecrise = 0.39782 * sind(Lrise);
cosDecrise = cosd(asind(sinDecrise));

sinDecset = 0.39782 * sind(Lset);
cosDecset = cosd(asind(sinDecset));

% 7a. calculate the Sun's local hour angle
	
cosHrise = (cosd(zenith)-(sinDecrise*sind(latitude)))/(cosDecrise*cosd(latitude));
cosHset = (cosd(zenith)-(sinDecset*sind(latitude)))/(cosDecset*cosd(latitude));

if cosHrise>1
    % the sun never rises on this location (on the specified date)
    Trise = 6;
else
    % 7b. finish calculating H and convert into hours if rising time is
    % desired:
    Hrise = 360 - acosd(cosHrise);
    % if setting time is desired: H = acos(cosH)
    Hrise = Hrise/15;

    % 8. calculate local mean time of rising/setting
    Trise = Hrise + RArise - 0.06571*trise - 6.622;
end

if cosHset<-1
    % the sun never sets on this location (on the specified date)
    Tset = 12+6;
else
    % 7b. finish calculating H and convert into hours if rising time is
    % desired: H = 360 - acos(cosH) if setting time is desired:
    Hset = acosd(cosHset);
    Hset = Hset/15;

    % 8. calculate local mean time of rising/setting
    Tset = Hset + RAset - 0.06571*tset - 6.622;
end


% 9. adjust back to UTC

UTrise = Trise - lngHour;
UTset = Tset - lngHour;

% 10. convert UT value to local time zone of latitude/longitude

localTrise = mod(UTrise + localOffset+dst*(N>=Nds && N<Nde), 24);
localTset = mod(UTset + localOffset+dst*(N>=Nds && N<Nde), 24);

end