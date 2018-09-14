function Dspk=getDspkPrmtrs(year,Doy,site) 
if nargin<3
    site=char.empty(0,0);
end
% Parameterd can be vectors of 1 2 3 4 6 or 12 parameters
% importing the paramaters from DespikeLimits_(year).m

eval(['Dspk = ' site 'DeSpikeLimitsTESTER();']);
dspkvarname=fieldnames(Dspk);
for i=1:length(dspkvarname)
    eval(['l = length (Dspk.' dspkvarname{i} ');']);
    if l == 12
        [ month , ~ ] = DayOY2Date( Doy , year ) ;
        eval([' Dspk.' dspkvarname{i} '= Dspk.' dspkvarname{i} '(' num2str(month) ');']);
    elseif l==4
        [ month , ~ ] = DayOY2Date( Doy , year ) ;
        if (month==12 || month==1 || month==2)
            eval([' Dspk.' dspkvarname{i} '= Dspk.' dspkvarname{i} '(1);']);
        elseif (month==3 || month==4 || month==5)
            eval([' Dspk.' dspkvarname{i} '= Dspk.' dspkvarname{i} '(2);']);
        elseif (month==6 || month==7 || month==8)
            eval([' Dspk.' dspkvarname{i} '= Dspk.' dspkvarname{i} '(3);']);    
        else 
            eval([' Dspk.' dspkvarname{i} '= Dspk.' dspkvarname{i} '(4);']); 
        end
    elseif l==1
    else
        disp('number of variables in despike limits is wrong');
    end
end
end
