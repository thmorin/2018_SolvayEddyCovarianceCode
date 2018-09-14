function [DayFormat , IsLeap] = IsLeapYear(year)

% Function get year and return Day format as a vector with the number of
% days of each month ralated to start of years in days,
% and returns 1=leap year 0=Non leap year

if mod(year,400)==0 
    DayFormat=[1 32 61 92 122 153 183 214 245 275 306 336];
    IsLeap = true;
elseif mod(year,100)==0
    DayFormat=[1 32 60 91 121 152 182 213 244 274 305 335];
    IsLeap = false;
elseif mod(year,4)==0
    DayFormat=[1 32 61 92 122 153 183 214 245 275 306 336];
    IsLeap = true;
else
    DayFormat=[1 32 60 91 121 152 182 213 244 274 305 335];
    IsLeap = false;
end

end