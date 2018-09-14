function [fluxOut]=FluxJoiner(flux1,flux2)
fluxOut=flux1;
Vars=fieldnames(flux1);
for i=1:length(Vars)
    eval(['fluxOut.' Vars{i} '(isnan(fluxOut.' Vars{i} '))=flux2.' Vars{i} '(isnan(flux1.' Vars{i} ')); '])
end
end