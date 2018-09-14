function [UstLimits,USTflag]=ustarfiltdataprocess(c_flux,u_star,t_ave,daynight,MinUstar)
USTflag=true(size(u_star));

[~,~,UstLimits]=despikeandfiltercor3(c_flux,u_star,daynight,t_ave);

if UstLimits<MinUstar
    disp(['UstLimits defaulted to MinUstar value of ' num2str(MinUstar) ' from value of ' num2str(UstLimits)]);
    UstLimits=MinUstar;
end

USTflag(u_star>UstLimits)=0;

end