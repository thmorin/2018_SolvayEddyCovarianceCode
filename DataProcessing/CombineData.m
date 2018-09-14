function CombineData(year)
OrigDir=['/home/morin.37/poolA/wetland_data/Matlab_Tim/DataProcessing/Processed' num2str(year) '/'];
SubDirs=dir([OrigDir num2str(year) '*']);
FinalDir=['/home/morin.37/poolA/wetland_data/Matlab_Tim/DataProcessing/Processed' num2str(year) '/flux/'];
delete([FinalDir 'OWR*']);

for i=1:length(SubDirs)
    files=dir([OrigDir num2str(SubDirs(i).name) '/flux/OWR*']);
    for j=1:length(files)
        load([OrigDir num2str(SubDirs(i).name) '/flux/' files(j).name]);
        if isempty(dir([FinalDir files(j).name]))
            save([FinalDir files(j).name],'fluxW');
        else
            flux1=fluxW;
            load([FinalDir files(j).name]);
            flux2=fluxW;
            clear fluxW;
            fluxW=MergeFluxFiles(flux1,flux2);
            save([FinalDir files(j).name],'fluxW');
        end
    end
end
end




function fluxOut = MergeFluxFiles(flux1,flux2)

fields=fieldnames(flux1);
for i=1:length(fields)
    fluxOut.(fields{i})=flux1.(fields{i});
    fluxOut.(fields{i})(isnan(fluxOut.(fields{i})))=flux2.(fields{i})(isnan(fluxOut.(fields{i})));
end
end