function FOOTSUM=FP_process_ORWRPbackup(year,ustar,ww_Vector,Lo,windir,zo,zH,d,FX,FY,Map,PlotDay,PlotTime)

% This process is computing the foot print at the ORWRP in half hour intervals.
% The foot print is "where the wind is coming from"
% Inputs for this process are : year
%                               ustar data from fluxW file, year long format
%                               ww data
%                               Obukhov length
%                               wind direction data
%                               roughness length and displacment hight
%                               (computed in the RoughLengthCalc.m)
%                               Map data, FX, FY, and the map itself
%                               Optional: PlotDay and PlotTime. Will plot the wind profile at date and time specified
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

if nargin==11
    PlotDay = nan; 
    PlotTime = nan;
end

NumValFile=48;
sv=sqrt(ww_Vector);

ly=leapyear(year);
DaysInYear=365+ly;

FOOTSUM=nan(NumValFile,length(PatchIndex(Map)),DaysInYear);
NanVec=nan(NumValFile,length(PatchIndex(Map)),DaysInYear);
for day=1:DaysInYear
    footsum_day = nan(NumValFile,length(PatchIndex(Map)));
    for Halfhour=1:NumValFile
        if PlotDay==day && PlotTime==Halfhour
%         if any(1:10:365==day)%BadFoot(day*NumValFile-NumValFile+Halfhour)==1
            plotYN=1;
        else 
            plotYN=0;
        end
        if isnan(ustar(day*NumValFile-NumValFile+Halfhour)*Lo(day*NumValFile-NumValFile+Halfhour)*...
                sv(day*NumValFile-NumValFile+Halfhour)*windir(day*NumValFile-NumValFile+Halfhour)*...
                zo(day*NumValFile-NumValFile+Halfhour)*zH(day*NumValFile-NumValFile+Halfhour)*d(day*NumValFile-NumValFile+Halfhour))
            
            NanVec(day*NumValFile-NumValFile+Halfhour)=1;
        else
            [~,~,~,~,~,~,footsum]=Hsieh2DRotate(ustar(day*NumValFile-NumValFile+Halfhour),...
                Lo(day*NumValFile-NumValFile+Halfhour),...
                sv(day*NumValFile-NumValFile+Halfhour),...
                zo(day*NumValFile-NumValFile+Halfhour),...
                zH(day*NumValFile-NumValFile+Halfhour),...
                d(day*NumValFile-NumValFile+Halfhour),...
                windir(day*NumValFile-NumValFile+Halfhour),...
                plotYN, FX, FY, Map);
            footsum_day(Halfhour,:) = footsum;
        end
    end
    FOOTSUM(:,:,day)=footsum_day;
%     if any(1:10:365==day)
        disp(['Footprint complete for day ' num2str(day)])
%     end
end
end