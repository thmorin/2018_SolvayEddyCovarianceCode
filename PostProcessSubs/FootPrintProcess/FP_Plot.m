year=2011;

% Path to map of ORWRP. Name of file is ORWRP_Map.mat name of matrix is WetlandMap
addpath /home/naor-azrieli.1/OlenData/OlentangyDataProcessMatlab/FootPrintProcess

% path to FP_process_ORWRP.m
directory = '/home/naor-azrieli.1/OlenData/OlentangyDataProcessMatlab/FootPrintProcess';

% Path to subroutins De_spike_gf and Hsieh2DRotate
addpath ( [ directory, '/Footprint/GapfillFootprint/'] )

% Path to IsLeapYear
addpath /home/naor-azrieli.1/OlenData/OlentangyDataProcessMatlab/Complete_Process_From_Gil/ProcessingSubRoutines/;

eval( ['load FootSum' num2str(year)] )

Days4Plotting = [78];

plotYN=1;
HH4Plotting = [30];
%footsum_day = nan(NumValFile,length(PatchIndex(WetlandMap)));
for day=Days4Plotting
    for Halfhour= HH4Plotting
        %[foot,xx,yy,footRotate,Fxprime,Fyprime,footsum]=Hsieh2DRotate
        [~,~,~,~,~,~,footsum]=Hsieh2DRotate_(ustar(day*NumValFile-NumValFile+Halfhour),...
            Lo(day*NumValFile-NumValFile+Halfhour),...
            sv(day*NumValFile-NumValFile+Halfhour),...
            zo,zH, d, windir(day*NumValFile-NumValFile+Halfhour),...
            plotYN, Fx1, Fy1, WetlandMap, day, Halfhour, year);
    end
end