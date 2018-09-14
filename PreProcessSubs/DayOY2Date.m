function [ month , day ] = DayOY2Date( DayOY , year )

%function [ month , day ] = DayOY2Date( DayOY , year ) take as input the Day Of year
% (depending if it's leap year or not) and return month and Day of month.

IfLeap = [ (1:366)' , ...
[ ones(31,1) ; ones(29,1)*2 ; ones(31,1)*3 ; ones(30,1)*4 ; ones(31,1)*5 ; ones(30,1)*6 ; ones(31,1)*7 ; ones(31,1)*8 ; ones(30,1)*9 ; ones(31,1)*10 ; ones(30,1)*11 ; ones(31,1)*12 ] , ...
[ (1:31)' ; (1:29)' ; (1:31)' ; (1:30)' ; (1:31)' ; (1:30)' ; (1:31)' ; (1:31)' ; (1:30)' ; (1:31)' ; (1:30)' ; (1:31)' ] ];

IfNotLeap = [ (1:365)' , ...
[ ones(31,1) ; ones(28,1)*2 ; ones(31,1)*3 ; ones(30,1)*4 ; ones(31,1)*5 ; ones(30,1)*6 ; ones(31,1)*7 ; ones(31,1)*8 ; ones(30,1)*9 ; ones(31,1)*10 ; ones(30,1)*11 ; ones(31,1)*12] , ...
[ (1:31)' ; (1:28)' ; (1:31)' ; (1:30)' ; (1:31)' ; (1:30)' ; (1:31)' ; (1:31)' ; (1:30)' ; (1:31)' ; (1:30)' ; (1:31)' ] ];

[ ~ , IsLeap] = IsLeapYear(year) ;

if IsLeap == 1
    month = IfLeap(DayOY,2) ;
    day = IfLeap(DayOY,3) ;
else
    month = IfNotLeap(DayOY,2) ;
    day = IfNotLeap(DayOY,3) ;
end

end


