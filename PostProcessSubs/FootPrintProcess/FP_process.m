function FOOTSUM=FP_process_ORWRP(ustar,vv_Vector,Lo,windir,zo,zH,d,FX,FY,Map,PlotVec)
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
%                               Optional: PlotVec. 0 for don't plot this
%                               half hour, 1 for plot it
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

if nargin==10
    PlotVec=zeros(size(ustar));
end

% NumValFile=48;
sv=sqrt(vv_Vector);

% ly=leapyear(year);
% DaysInYear=365+ly;

FOOTSUM=nan(length(ustar),length(PatchIndex(Map)));
parfor i=1:length(ustar)
    if ~isnan(ustar(i)*Lo(i)*sv(i)*zo*zH*d*windir(i))
        [~,~,~,~,~,~,FOOTSUM(i,:)]=Hsieh2DRotate(ustar(i),Lo(i),sv(i),zo,zH,d,windir(i),PlotVec(i), FX, FY, Map);
    end
    if mod(i,48*10)==0    
        disp(['Footprint complete for day ' num2str(floor(i/48))])
    end
end
end