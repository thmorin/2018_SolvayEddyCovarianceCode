function [extra_vars,DOY,DECDAY,DECYear,tair,tair_gf,pressure,pressure_gf,Fc,Fc_gf,LE,LE_gf,SWinbar,SWinbar_gf,SWoutbar,SWoutbar_gf,net_sw,net_sw_gf,LWinbar,LWinbar_gf,LWoutbar,LWoutbar_gf,net_lw,net_lw_gf,Q,Q_gf,...
    Total_PAR,RNET,RNET_gf,rH,rH_gf,vpd,vpd_gf,L,L_gf,ww,ww_gf,vv_1,vv_1_gf,H,H_gf,wind_speed_1,wind_speed_1_gf,wind_dir_1,wind_dir_1_gf,ustar_1,ustar_1_gf,...
    wind_speed_2,wind_speed_2_gf,wind_dir_2,wind_dir_2_gf,ustar_2,ustar_2_gf,vv_2,vv_2_gf,Spectral_correction,C]=....
    Build_arrays(datadirectory,NumWindMeasures,loops,Dspk,start_year,start_DOY,end_year,end_DOY,varargin)
%Sylvania_build_arrays
%Written by Tim Morin
%Last edited July, 2013
%
%datadirectory: where fluxW files are stored. Must in format of 'OWR_year_DOY_flux.mat'
%Code must be changed if files do not appear in this format
%
%NumWindMeasures: 1 or 2. Determines if data from an upper wind sensor
%should be found or not. Also allows for gapfilling of lower sensor from
%upper sensor data and vice versa
%
%start_year: year at the beginning of the season of interest
%start_DOY: DOY at the beginning of the season of interest (eg. 270)
%end_year: year at the end of the season of interest
%end_DOY: DOY at the end of hte season of interest
%
%nargin: non standard variables that should also be created. Format: 'NameOfVariable1','NameOfVariable2',etc...

extra_vars=struct;

if start_year==end_year
    totallength=end_DOY-start_DOY+1;
    DateList=nan(totallength,2);
    DateList(:,1)=start_year;
    DateList(:,2)=start_DOY:end_DOY;
else    
    ly1=leapyear(start_year);
    year1total=(365+leapyear(start_year)-start_DOY)+1;
    year2total=end_DOY;
    totallength=year1total+year2total;
    DateList=nan(totallength,2);
    DateList(1:year1total,1)=start_year;
    DateList(year1total+1:end,1)=end_year;
    DateList(1:year1total,2)=start_DOY:365+ly1;
    DateList(year1total+1:end,2)=1:end_DOY;
end

[DOY,DECDAY,DECYear,C,Fc,H,L,LE,LWinbar,LWoutbar,net_lw,net_sw,pressure,...
    p_vapor_bar,p_vaporSat_bar,Q,rho_cp,rH,RNET,~,Spectral_correction,SWinbar,SWoutbar,tair,ustar_1,...
    ~,wts,ww,vv_1,wind_dir_1,wind_speed_1,ustar_2,vv_2,wind_dir_2,wind_speed_2,Total_PAR,Diffuse_PAR]=...
    InitializeArrays(totallength,loops,NumWindMeasures);

if ~isempty(varargin)
    for i=1:length(varargin)
       eval(['extra_vars.' varargin{i} '= nan(totallength,1);']);
    end
end

for i=1:length(DateList)
    year=DateList(i,1);
    currentday=DateList(i,2);
    files=dir([datadirectory 'Processed' num2str(year) '/flux/*' num2str(currentday,'%03.0f') '_flux.mat']);
    if ~isempty(files)
        load([datadirectory 'Processed' num2str(year) '/flux/' files(1).name])
        
        Spectral_correction((i-1)*48+1:i*48) = fluxW.Spectral_correction;

        Fc((i-1)*48+1:i*48) = fluxW.Fc;

        if isfield(fluxW,'ubar_1')
            wind_speed_1((i-1)*48+1:i*48) = fluxW.ubar_1;
        elseif isfield(fluxW,'ubar')
            wind_speed_1((i-1)*48+1:i*48) = fluxW.ubar;
        end

        LE((i-1)*48+1:i*48) = fluxW.LE;
        if isfield(fluxW,'ww_1')
            ww((i-1)*48+1:i*48) = fluxW.ww_1;
        elseif isfield(fluxW,'ww')
            ww((i-1)*48+1:i*48) = fluxW.ww;
        end

        if isfield(fluxW,'WD_1_degN')
            if ~all(isnan(fluxW.WD_1_degN))
                wind_dir_1((i-1)*48+1:i*48) = fluxW.WD_1_degN; %[deg]
            else
                wind_dir_1((i-1)*48+1:i*48) = fluxW.WD_1*180/pi;
            end
        elseif isfield(fluxW,'WDdegN')    
            wind_dir_1((i-1)*48+1:i*48) = fluxW.WDdegN; %[deg]
        else
            wind_dir_1((i-1)*48+1:i*48) = fluxW.WD_1*180/pi;
        end
        
        if isfield(fluxW,'vv_1')
            vv_1((i-1)*48+1:i*48) = fluxW.vv_1;
        elseif isfield(fluxW,'vv')
            vv_1((i-1)*48+1:i*48) = fluxW.vv;
        end
        
        if isfield(fluxW,'Total_PAR')
            Total_PAR((i-1)*48+1:i*48) = fluxW.Total_PAR;
        end
        
        if isfield(fluxW,'Diffuse_PAR')
            Diffuse_PAR((i-1)*48+1:i*48) = fluxW.Diffuse_PAR;
        end
        
        Spectral_correction((i-1)*48+1:i*48) = fluxW.Spectral_correction;
        
        tair((i-1)*48+1:i*48) = fluxW.tair;
        
        p_vaporSat_bar((i-1)*48+1:i*48) = fluxW.p_vaporSat_bar; %Used to calculate vpd
        rH((i-1)*48+1:i*48) = fluxW.humbar;
        p_vapor_bar((i-1)*48+1:i*48) = fluxW.p_vapor_bar;
        
        pressure((i-1)*48+1:i*48) = fluxW.pressure;
        
        if isfield(fluxW,'ustar_1')
            ustar_1((i-1)*48+1:i*48) = fluxW.ustar_1;
        elseif isfield(fluxW,'ustar')
            ustar_1((i-1)*48+1:i*48) = fluxW.ustar;
        end
                
        H((i-1)*48+1:i*48) = fluxW.H;
        
        LWinbar((i-1)*48+1:i*48) = fluxW.LWinbar;
        LWoutbar((i-1)*48+1:i*48) = fluxW.LWoutbar;
        
        net_sw((i-1)*48+1:i*48) = fluxW.NetRs;
        net_lw((i-1)*48+1:i*48) = fluxW.NetRl;
                
        RNET((i-1)*48+1:i*48) = fluxW.NetRad;
        
        rho_cp((i-1)*48+1:i*48) = fluxW.rho_cp;
        wts((i-1)*48+1:i*48) = fluxW.wts;        
        
        L((i-1)*48+1:i*48) = fluxW.L;
        
        SWinbar((i-1)*48+1:i*48) = fluxW.SWinbar;
        SWoutbar((i-1)*48+1:i*48) = fluxW.SWoutbar;
        
        Q((i-1)*48+1:i*48) = fluxW.Q;
        
        C((i-1)*48+1:i*48) = fluxW.C;
        
        wts((i-1)*48+1:i*48) = fluxW.wts;
        
        if ~isempty(varargin)
            for j=1:length(varargin)
                if isfield(fluxW,varargin{j})
                    %disp(['Making variable ' varargin{j} ' for ' num2str(currentday)]);
                    eval(['extra_vars.' varargin{j} '((i-1)*48+1:i*48) = fluxW.' varargin{j} ';'])
                else
                    warning(['Field ' varargin{j} ' does not exist in fluxW.']);
                end
            end
        end
        if NumWindMeasures==2
            if isfield(fluxW,'WD_2_degN')
                wind_dir_2((i-1)*48+1:i*48) = fluxW.WD_2_degN; %[deg]
            elseif isfield(fluxW,'WD_2') 
                wind_dir_2((i-1)*48+1:i*48) = fluxW.WD_2*180/pi; %[deg]
            end

            if isfield(fluxW,'vv_2')
                vv_2((i-1)*48+1:i*48) = fluxW.vv_2;
            end

            if isfield(fluxW,'ustar_2')
                ustar_2((i-1)*48+1:i*48) = fluxW.ustar_2;
            end
            if isfield(fluxW,'wind_speed_2')
                wind_speed_2((i-1)*48+1:i*48) = fluxW.ubar_2;
            end
        end
    end
    DOY((i-1)*48+1:i*48) = currentday;
    DECDAY((i-1)*48+1:i*48) = currentday:1/48:currentday+47/48;
    DECYear((i-1)*48+1:i*48)=year+DECDAY((i-1)*48+1:i*48)/(365+leapyear(year));
end

%%%%%%%%%%%%%%%VERY TEMPORARY FIX%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
wind_dir_2=mod(-wind_dir_2-90,360);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

SWinbar=DeSpike(SWinbar,Dspk.Interval,Dspk.STD,Dspk.SWinbar.min,Dspk.SWinbar.max,Dspk.trimp);
if sum(~isnan(SWinbar))~=0
	[~,SWinbar_gf]=gapfill(SWinbar,totallength,48);
else
	SWinbar_gf=nan(size(SWinbar));
end

SWoutbar=DeSpike(SWoutbar,Dspk.Interval,Dspk.STD,Dspk.SWoutbar.min,Dspk.SWoutbar.max,Dspk.trimp);
if sum(~isnan(SWoutbar))~=0
	[~,SWoutbar_gf]=gapfill(SWoutbar,totallength,48);
else
	SWoutbar_gf=nan(size(SWoutbar));
end

LWinbar=DeSpike(LWinbar,Dspk.Interval,Dspk.STD,Dspk.LWinbar.min,Dspk.LWinbar.max,Dspk.trimp);
[~,LWinbar_gf]=gapfill(LWinbar,totallength,48);

LWoutbar=DeSpike(LWoutbar,Dspk.Interval,Dspk.STD,Dspk.LWoutbar.min,Dspk.LWoutbar.max,Dspk.trimp);
[~,LWoutbar_gf]=gapfill(LWoutbar,totallength,48);

if sum(~isnan(Q))~=0
	[~,Q_gf]=gapfill(Q,totallength,48);
else
	Q_gf=nan(size(Q));
end

if sum(~isnan(Fc))~=0
	[~,Fc_gf]=gapfill(Fc,totallength,48);
else
	Fc_gf=nan(size(Fc));
end

pressure=DeSpike(pressure,Dspk.Interval,Dspk.STD,Dspk.pressure.min,Dspk.pressure.max,Dspk.trimp);
if sum(~isnan(pressure))~=0
	[~,pressure_gf]=gapfill(pressure,totallength,48);
else
	pressure_gf=nan(size(pressure));
end

tair=DeSpike(tair,Dspk.Interval,Dspk.STD,Dspk.tair.min,Dspk.tair.max,Dspk.trimp);
[~, tair_gf]=gapfill(tair,totallength,48);

net_sw=DeSpike(net_sw,Dspk.Interval,Dspk.STD,Dspk.net_sw.min,Dspk.net_sw.max,Dspk.trimp);
[~, net_sw_gf]=gapfill(net_sw,totallength,48);

net_lw=DeSpike(net_lw,Dspk.Interval,Dspk.STD,Dspk.net_lw.min,Dspk.net_lw.max,Dspk.trimp);
[~, net_lw_gf]=gapfill(net_lw,totallength,48);

wind_speed_1=DeSpike(wind_speed_1,Dspk.Interval*2,Dspk.STD,Dspk.wind_speed_1.min,Dspk.wind_speed_1.max,Dspk.trimp);
if sum(~isnan(wind_speed_1))~=0
    [~, wind_speed_1_gf]=gapfill(wind_speed_1,totallength,48);
else
    wind_speed_1_gf=nan(size(wind_speed_1));
end

if sum(~isnan(wind_speed_2))~=0
    wind_speed_2=DeSpike(wind_speed_2,Dspk.Interval*2,Dspk.STD,Dspk.wind_speed_2.min,Dspk.wind_speed_2.max,Dspk.trimp);
    [~, wind_speed_2_gf]=gapfill(wind_speed_2,totallength,48);
else
    wind_speed_2_gf=nan(size(wind_speed_2));
end

p_vaporSat_bar=DeSpike(p_vaporSat_bar,Dspk.Interval*2,Dspk.STD,Dspk.p_vaporSat_bar.min,Dspk.p_vaporSat_bar.max,Dspk.trimp);
[~, p_vaporSat_bar_gf]=gapfill(p_vaporSat_bar,totallength,48);

wind_dir_1_gf=wind_dir_1;wind_dir_1_gf(isnan(wind_dir_1))=wind_dir_2(isnan(wind_dir_1));
wind_dir_2_gf=wind_dir_2;wind_dir_2_gf(isnan(wind_dir_2))=wind_dir_1(isnan(wind_dir_2));
if sum(~isnan(wind_dir_1_gf))~=0
    [~, wind_dir_1_gf]=gapfill(wind_dir_1_gf,totallength,48);
else
    wind_dir_1_gf=nan(size(wind_dir_1));
end

if sum(~isnan(wind_dir_2_gf))~=0
	[~, wind_dir_2_gf]=gapfill(wind_dir_2_gf,totallength,48);
end

ustar_1=DeSpike(ustar_1,Dspk.Interval,Dspk.STD,Dspk.ustar_1.min,Dspk.ustar_1.max,Dspk.trimp);
ustar_1_gf=ustar_1;

if sum(~isnan(ustar_2))~=0
    ustar_2=DeSpike(ustar_2,Dspk.Interval,Dspk.STD,Dspk.ustar_2.min,Dspk.ustar_2.max,Dspk.trimp);
    ustar_2_gf=ustar_2;
    ustar_1_gf(isnan(ustar_1_gf))=ustar_2(isnan(ustar_1_gf));
    ustar_2_gf(isnan(ustar_2_gf))=ustar_1(isnan(ustar_2_gf));
    [~,ustar_2_gf]=gapfill(ustar_2_gf,totallength,48);
else
    ustar_2_gf=ustar_2;
end

if sum(~isnan(ustar_1_gf))~=0
	[~,ustar_1_gf]=gapfill(ustar_1_gf,totallength,48);
end
wts=DeSpike(wts,Dspk.Interval,Dspk.STD,Dspk.wts.min,Dspk.wts.max,Dspk.trimp);
if sum(~isnan(wts))~=0
	[~,wts_gf]=gapfill(wts,totallength,48);
else
	wts_gf=nan(size(wts));
end

rho_cp=DeSpike(rho_cp,Dspk.Interval,Dspk.STD,Dspk.rho_cp.min,Dspk.rho_cp.max,Dspk.trimp);
[~,rho_cp_gf]=gapfill(rho_cp,totallength,48);

H =DeSpike(H,Dspk.Interval,Dspk.STD,Dspk.H.min,Dspk.H.max,Dspk.trimp);
H_gf=H;
H_gf(isnan(H_gf))=wts_gf(isnan(H_gf)).*rho_cp_gf(isnan(H_gf));

LE=DeSpike(LE,Dspk.Interval,Dspk.STD,Dspk.LE.min,Dspk.LE.max,Dspk.trimp);
if sum(~isnan(LE))~=0
	[~,LE_gf]=gapfill(LE,totallength,48);
else
	LE_gf=nan(size(LE));
end

RNET=DeSpike(RNET,Dspk.Interval,Dspk.STD,Dspk.RNET.min,Dspk.RNET.max,Dspk.trimp);
RNET_gf=RNET;
RNET_gf(isnan(RNET_gf))=net_sw_gf(isnan(RNET_gf))+net_lw_gf(isnan(RNET_gf));

p_vapor_bar=DeSpike(p_vapor_bar,Dspk.Interval/2,Dspk.STD,Dspk.p_vapor_bar.min,Dspk.p_vapor_bar.max,Dspk.trimp);
p_vapor_bar_gf=p_vapor_bar;
p_vapor_bar_gf(isnan(p_vapor_bar_gf))=611.2*exp(17.67*tair_gf(isnan(p_vapor_bar_gf))...
                                                      ./(243.51+tair_gf(isnan(p_vapor_bar_gf))));

rH=DeSpike(rH,Dspk.Interval,Dspk.STD,Dspk.rH.min,Dspk.rH.max,Dspk.trimp);
rH_gf=rH;
rH_gf(isnan(rH_gf))=p_vapor_bar(isnan(rH_gf))./p_vaporSat_bar(isnan(rH_gf));
[~,rH_gf]=gapfill(rH_gf,totallength,48);
rH_gf(rH_gf>100)=100;
rH_gf(rH_gf<0)=0;
                                                  
vpd=p_vaporSat_bar - p_vapor_bar;
vpd_gf=p_vaporSat_bar_gf - p_vapor_bar_gf;
vpd_gf(vpd_gf<0)=0;

L=DeSpike(L,Dspk.Interval,Dspk.STD,Dspk.L.min,Dspk.L.max,Dspk.trimp);
L_gf=L;
L_gf(isnan(L))= -(ustar_1_gf(isnan(L)).^3) ./ (0.4*9.81*( wts_gf(isnan(L)) ./ (tair_gf(isnan(L))+273.15)));
if sum(~isnan(LE_gf))~=0
	[~,L_gf]=gapfill(L_gf,totallength,48);
end

ww=DeSpike(ww,Dspk.Interval,Dspk.STD,Dspk.ww.min,Dspk.ww.max,Dspk.trimp);
if sum(~isnan(ww))~=0
    [~,ww_gf]=gapfill(ww,totallength,48);
else
    ww_gf=nan(size(ww));
end

vv_1=DeSpike(vv_1,Dspk.Interval,Dspk.STD,Dspk.vv.min,Dspk.vv.max,Dspk.trimp);

vv_1_gf=vv_1;vv_2_gf=vv_2;
if sum(~isnan(vv_1_gf))~=0
[~,vv_1_gf]=gapfill(vv_1_gf,totallength,48);
end

if sum(~isnan(vv_2))~=0
	vv_1_gf(isnan(vv_1_gf))=vv_2(isnan(vv_1_gf));
	vv_2=DeSpike(vv_2,Dspk.Interval,Dspk.STD,Dspk.vv.min,Dspk.vv.max,Dspk.trimp);
	vv_2_gf(isnan(vv_2_gf))=vv_1(isnan(vv_2_gf));
	[~,vv_2_gf]=gapfill(vv_2_gf,totallength,48);
end

if ~isempty(varargin)
    
    for j=1:length(varargin)
        if isfield(fluxW,varargin{j}) && isfield(Dspk,varargin{j})
            eval(['extra_vars.' varargin{j} '=DeSpike(extra_vars.' varargin{j} ',Dspk.Interval,Dspk.STD,Dspk.' varargin{j} '.min,Dspk.' varargin{j} '.max,Dspk.trimp);'])
            eval(['Checker=sum(~isnan(extra_vars.' varargin{j} '));']);
            if Checker>5
                eval(['[~,extra_vars.' varargin{j} '_gf]=gapfill(extra_vars.' varargin{j} ',totallength,48);']);
            end
        end
    end
end

end
