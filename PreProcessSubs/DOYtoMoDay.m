function [month day]=DOYtoMoDay(CD,year)
    leap=leapyear(year)-365;
    if CD<1+31
        month=1;
        day=CD;
    elseif CD<1+31+28+leap
        month=2;
        day=CD-31;
    elseif CD<91+leap
        month=3;
        day=CD-59-leap;
    elseif CD<121+leap
        month=4;
        day=CD-90-leap;
    elseif CD<152+leap
        month=5;
        day=CD-120-leap;
    elseif CD<182+leap;
        month=6;
        day=CD-151-leap;
    elseif CD<213+leap
        month=7;
        day=CD-181-leap;
    elseif CD<244+leap
        month=8;
        day=CD-212-leap;
    elseif CD<274+leap
        month=9;
        day=CD-243-leap;
    elseif CD<305+leap
        month=10;
        day=CD-273-leap;
    elseif CD<335+leap
        month=11;
        day=CD-304-leap;
    else
        month=12;
        day=CD-334-leap;
    end
end