function [errorout] = Error_calc(fluxes,par,temp,ubar,knownerror)

if nargin==5
    knownerror=nan(size(fluxes));
end

errorout=nan(size(fluxes));

for i=1:size(fluxes,2)
    [errorout(:,i),~] = RichHoll_FluxError(fluxes(:,i),par,temp,ubar);
end

errorout(~isnan(knownerror))=knownerror(~isnan(knownerror));

end