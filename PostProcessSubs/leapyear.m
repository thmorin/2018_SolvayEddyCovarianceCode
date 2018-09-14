function ly=leapyear(year)
if (mod(year,4)==0 && ~(mod(year,100)==0)) || mod(year,400)==0
    ly=true;
else
    ly=false;
end
end