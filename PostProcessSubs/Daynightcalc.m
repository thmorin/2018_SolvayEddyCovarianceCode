function [daynight]=Daynightcalc(DECYear,par,lat,lon,GMToffset)

marker = true(size(par));
daynight = par>10;

i=1;
while i < length(DECYear)
    year = floor(DECYear(i));
    currentday = floor((DECYear(i)-year)*(365+leapyear(year)));
    
    daynight(i:i+10)=0;
    daynight(i+42:i+47)=0;

    daynight(i+16:i+34)=1;

    [localTrise, localTset] = Suntime(lat,lon,GMToffset,0,year,currentday);

    marker(i:i+round(localTrise*2))=0;
    marker(i+round(localTset*2):i+47)=0;
    
    i=i+48;
end
daynight(isnan(daynight))=marker(isnan(daynight));

end