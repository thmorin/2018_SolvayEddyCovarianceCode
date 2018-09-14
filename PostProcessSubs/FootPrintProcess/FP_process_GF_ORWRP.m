clear all
%

% This process is computing the foot print at the ORWRP in half hour intervals.
% The foot print is "where the wind is coming from"
% Inputs for this process are : flux data computed through DataProces201*.m
%                               roughness length and displacment hight
%                               (computed in the RoughLengthCalc.m)
% Before trying to process the flux data, copy all *flux.mat into a folder
% using the following comand in terminal:
% cp /home/naor-azrieli.1/OlenData/OlentangyDataProcessMatlab/Complete_Process_From_Gil/processed2011_03/*[f]* /home/naor-azrieli.1/OlenData/OlentangyDataProcessMatlab/Complete_Process_From_Gil/processed2011
% cp /home/naor-azrieli.1/OlenData/OlentangyDataProcessMatlab/Complete_Process_From_Gil/processed2011_08/*[f]* /home/naor-azrieli.1/Desktop/processed2011/
%
% The output of this process is a matrix the atributes percent of "where the wind is coming from"
% to the different patch type with the following code:
% Patch type code to patch canopy properties data in ForestCanopy_data.m 
% Type code in ForestCanopy_data -> Type code in CompisiteMeters.txt
% Code for ORWRP map: 
% 1 = Water       -> 333 Water
% 2 = Pavement    -> 309 Pavement
% 3 = Grass       -> 338 Grass
% 4 = Small trees -> 205 Trees 5 to 10 meters
% 5 = Tall trees  -> 210 Trees 10 to 15 meters (tall trees  -> 215 Trees 15 to 20 meters ;-> 220 Trees 20 to 25 meters)
% 6 = Buildings   -> 314 Buildings
% 7 =
% 8 = Wetland 1 Open Water
% 9 = Wetland 2 Open water
% 10= Location of Tower

% Path to map of ORWRP. Name of file is ORWRP_Map.mat name of matrix is WetlandMap
addpath /home/naor-azrieli.1/OlenData/OlentangyDataProcessMatlab/FootPrintProcess

% process_Period = '_06' ;  % If I want to process by periods with same
% data structure (same columns)

% path to FP_process_ORWRP.m
directory = '/home/naor-azrieli.1/OlenData/OlentangyDataProcessMatlab/FootPrintProcess';

% Path to subroutins De_spike_gf and Hsieh2DRotate
addpath ( [ directory, '/Footprint/GapfillFootprint/'] )

% Path to IsLeapYear
addpath /home/naor-azrieli.1/OlenData/OlentangyDataProcessMatlab/Complete_Process_From_Gil/ProcessingSubRoutines/;

year=2011;

% Days and half-hours for plotting
% which means 12:00 on Day of year 75

PlotDay = 75; 
PlotTime = 24;

% Path to Outcome
ProcessedDataDir = [directory,'/ProcessedFP' num2str(year) '/'];

% Path to *flux.mat files folder by year
% DataDirectory = ['/home/naor-azrieli.1/OlenData/OlentangyDataProcessMatlab/...
%  Complete_Process_From_Gil/processed' num2str(year) process_Period '/'];
DataDirectory = [directory '/flux' num2str(year) '/'];


site = 'OWR';                           % Olentangy Wetland Research
FluxFile='flux.mat';
% load([datadirectory site '_' num2str(year) '_' num2str(currentday,'%03.0f') '_flux.mat'],'fluxW');

NumODays = dir([ DataDirectory '*' FluxFile ]);    % List all ts files in directory
[numoff,~]= size(NumODays);
DaysInFolder = nan*(1:numoff);

for DF=1:numoff
    a1=regexprep(NumODays(DF).name,'.mat','');
    a=isletter(a1);
    b=NumODays(DF).name(not(a));
    c=regexprep(b,'_','');
    d=regexprep(c,num2str(year),'');
    DaysInFolder(DF)=str2double(d);
end
    

[DayFormat , IsLeap] = IsLeapYear(year);
DaysInYear=365+IsLeap; % Isleap=1 for leap year and 0 otherwise.

NumValFile = 48;

%  fluxW
flux = {'ustar','ww_9_6m','L','WD_9_6_degN'};

%flux = {'tair','tvair','humbar','NetRad','ubar_9_6m','wbar_9_6m','ustar_9_6m','ustar','pressure','p_vapor_bar',...
%    'r','qbar','cbar','tbar','LE','H','Fc','wts','rho','wq','wc'};

for i=1:length(flux)
    
     eval([flux{i} 'Vector = ' num2str(nan) '( length( ' num2str(1) ':' num2str(DaysInFolder(end)-DaysInFolder(1)+1) ' ) *' num2str(NumValFile) ',1 );' ]);
    %eval([flux{i} 'Vector = ' num2str(nan) '(' num2str(NumValFile) ', length( ' num2str(1) ...
      %  ':' num2str(DaysInFolder(end)-DaysInFolder(1)+1) ' ) );' ]);
end

for CD = DaysInFolder
    % It loads the flux.mat  - Name format is : OWR_2011_091_flux.mat
    load([DataDirectory site '_' num2str(year) '_' num2str(CD,'%03.0f') '_flux.mat'],'fluxW');
    CV = ( (CD-DaysInFolder(1))*NumValFile+1 : (CD-DaysInFolder(1)+1)*NumValFile )'; % CFV = Current Fast Values, Make a vector of the fast values to work on for this Current Window
    for i=1:length(flux)
        eval([flux{i} 'Vector(' num2str((CD-DaysInFolder(1))*NumValFile+1) ':' num2str((CD-DaysInFolder(1)+1)*NumValFile) ',1)' '=fluxW.' flux{i} ';']);
        %eval([flux{i} 'Vector( : ,' num2str(CD) '-' ...
            %num2str(DaysInFolder(1)) '+1 )' '=fluxW.' flux{i} ';']);
    end
        
end

% ustar = ustarVector ;
% Lo = LVector ;
% sv = sqrt(ww_9_6mVector) ;
% windir = WD_9_6_degNVector ;



[~ , ustarVector] = De_spike_gf(ustarVector,(DaysInFolder(end)-DaysInFolder(1)+1)*48,6,48);
[~ , LVector] = De_spike_gf(LVector,(DaysInFolder(end)-DaysInFolder(1)+1)*48,6,48);
[~ , svVector] = De_spike_gf(sqrt(ww_9_6mVector),(DaysInFolder(end)-DaysInFolder(1)+1)*48,6,48);
[~ , WD_9_6_degNVector] = De_spike_gf(WD_9_6_degNVector,(DaysInFolder(end)-DaysInFolder(1)+1)*48,6,48);

ustar = nan(DaysInYear*NumValFile,1);
Lo = nan(DaysInYear*NumValFile,1);
sv = nan(DaysInYear*NumValFile,1);
windir = nan(DaysInYear*NumValFile,1);

ustar((DaysInFolder(1)-1)*NumValFile+1:DaysInFolder(end)*NumValFile,1)=ustarVector;
Lo((DaysInFolder(1)-1)*NumValFile+1:DaysInFolder(end)*NumValFile,1)=LVector;
sv((DaysInFolder(1)-1)*NumValFile+1:DaysInFolder(end)*NumValFile,1)=svVector;
windir((DaysInFolder(1)-1)*NumValFile+1:DaysInFolder(end)*NumValFile,1)=WD_9_6_degNVector;

CLEARVARS -EXCEPT windir sv Lo ustar DaysInFolder ...
    NumValFile DaysInYear PlotDay PlotTime year
H=8; % Canopy height
zo=0.15*H; zH=9.6; d=0.67*H; plotYN=0;
%       Lo      = Obukhov stability length [m]
%       zo      = momentum roughness length [m]
%       zH      = height of the measurement [m]
%       d - displacement height

load('ORWRP_Map.mat')


FOOTSUM=nan(NumValFile,length(PatchIndex(WetlandMap)),DaysInYear);
footsum_day = nan(NumValFile,length(PatchIndex(WetlandMap)));
for day=DaysInFolder
    for Halfhour=1:NumValFile
        %[foot,xx,yy,footRotate,Fxprime,Fyprime,footsum]=Hsieh2DRotate(ustar(runs),Lo(runs),sv(runs),zo,zH, d, windir(runs), plotYN, Fx1, Fy1, WetlandMap);
        if PlotDay==day && PlotTime==Halfhour
            plotYN=1;
        else 
            plotYN=0;
        end
        [~,~,~,~,~,~,footsum]=Hsieh2DRotate_(ustar(day*NumValFile-NumValFile+Halfhour),...
            Lo(day*NumValFile-NumValFile+Halfhour),...
            sv(day*NumValFile-NumValFile+Halfhour),...
            zo,zH, d, windir(day*NumValFile-NumValFile+Halfhour),...
            plotYN, Fx1, Fy1, WetlandMap, Day, Halfhour, year);
        footsum_day(Halfhour,:) = footsum;
        clear foot xx yy footRotate Fxprime Fyprime footsum
    end
    FOOTSUM(:,:,day)=footsum_day;
    
    if any(1:10:365==day)
        disp(day)
    end
end
eval( ['save FootSum' num2str(year)] )
% figure(101)
% pcolor(XX,YY,FF);shading('interp')
% hold on
% contour(Fxprime,Fyprime,F,[0 0],'LineColor',[1 1 1])
% axis([-200 500 -500 500])

Days4Plotting = [78];
plotYN=1;
HH4Plotting = [30];
FOOTSUM=nan(NumValFile,length(PatchIndex(WetlandMap)),length(DaysInFolder)+1);
footsum_day = nan(NumValFile,length(PatchIndex(WetlandMap)));
for day=Days4Plotting
    for Halfhour= HH4Plotting
        %[foot,xx,yy,footRotate,Fxprime,Fyprime,footsum]=Hsieh2DRotate
        [~,~,~,~,~,~,footsum]=Hsieh2DRotate_(ustar(day*NumValFile-NumValFile+Halfhour),...
            Lo(day*NumValFile-NumValFile+Halfhour),...
            sv(day*NumValFile-NumValFile+Halfhour),...
            zo,zH, d, windir(day*NumValFile-NumValFile+Halfhour),...
            plotYN, Fx1, Fy1, WetlandMap, Day, Halfhour, year);
    end
end