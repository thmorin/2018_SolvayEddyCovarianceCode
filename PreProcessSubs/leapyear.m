function [days]=leapyear(year)
if mod(year,4)==0 && (mod(year,100)~=0 || (mod(year,400)==0))
    days=366;
else
    days=365;
end
end