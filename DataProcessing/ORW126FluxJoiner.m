load(([savedir 'flux/OWR_2013_126_1_flux.mat']));flux1=fluxW;
load(([savedir 'flux/OWR_2013_126_2_flux.mat']));flux2=fluxW;
clear fluxW
Vars=fieldnames(flux1);
BreakPt=24;

fluxW=flux1;
for i=1:length(Vars)
    eval(['fluxW.' Vars{i} '(BreakPt+1:end)=flux2.' Vars{i} '(BreakPt+1:end);'])
end

movefile([savedir 'flux/OWR_2013_126_1_flux.mat'],[savedir 'flux/orig/OWR_2013_126_1_flux.mat']);
movefile([savedir 'flux/OWR_2013_126_2_flux.mat'],[savedir 'flux/orig/OWR_2013_126_2_flux.mat']);
save([savedir 'flux/OWR_2013_126_flux.mat'],'fluxW');