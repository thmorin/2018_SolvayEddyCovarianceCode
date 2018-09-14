% process data
% 
%clear all

% when switching between 46m to 34m, changes need to be made for
% datadirectory, procdir, site, filelists, data2 read file (twice), fast read
% file

% CHECK "YEAR" IN ALL READ FILES

addpath('/old-drive/Processdata/ProcessingSubRoutines/')

load('listoffiles_faset')

IRGA = 2;                               % 1='Open', 2='Closed'(AF46m), 3='Closed'(FASET32m) VERY IMPORTANT$&$&$&$&$&$&$$&$&$&$&$$&$
irgflg = 2;                             % 1=LI-7000, 2=LI-6262
TubeD = 0.004725;                       % IRGA tube inner diameter [m]
LT = 41.5;                              % IRGA tube length diameter [m]
z_geom = 32;                            % Sensor height [m]
lso = 0.3;                              % Line averaging length sonic [m]

year = 2011;
rowID = year-2006;                      %tells program which row to look for files in
preclevel = 2;                          %preclevel 1= milisec; 2=sec;3=min
nhead = 0;                              %number of header lines
datadirectory = '/home/maurer.189/FASET/';
procdir = ['/home/maurer.189/FASET/' num2str(year) 'Proc/'];
site = 'Faset';                         % 46mAF for 46m and 34mAF for 34m
nhours = 24;                            %number of hours per day
nhoursinfile = 1;                       %number of hours per file
nhalfhours = 48;                        %number of half hours per day
frec = 10;                              %data frequency in Hz
AVGtime = 1800 ;                        %averging window in seconds


nts = nhours*3600*frec;                 % # of timesteps per day
ntsf = nhoursinfile*3600*frec;          % # of timesteps per file

fs = frec;                              % frequency of data measurements [Hz]
dt=1/frec;                              % data interval [sec]
nn = frec*AVGtime;                      % window size 
DT=AVGtime;                             % averaging interval [sec]
nw=floor(nhours*3600/AVGtime);          % number of overlapping windows (and spectra)      

fluxFA.ubar = ones(nw,1)*nan;
fluxFA.wbar= ones(nw,1)*nan;
fluxFA.tsvbar= ones(nw,1)*nan;
fluxFA.uu= ones(nw,1)*nan;
fluxFA.vv= ones(nw,1)*nan;
fluxFA.ww= ones(nw,1)*nan;
fluxFA.uw= ones(nw,1)*nan;
fluxFA.vw= ones(nw,1)*nan;
fluxFA.u3= ones(nw,1)*nan;
fluxFA.v3= ones(nw,1)*nan;
fluxFA.w3= ones(nw,1)*nan;
fluxFA.ustar= ones(nw,1)*nan;
fluxFA.L= ones(nw,1)*nan;
fluxFA.WD = ones(nw,1)*nan;
fluxFA.tair = ones(nw,1)*nan;
fluxFA.tvair  = ones(nw,1)*nan;
fluxFA.pbar  = ones(nw,1)*nan;
fluxFA.humbar = ones(nw,1)*nan;
fluxFA.Pvaporbar = ones(nw,1)*nan;
fluxFA.Pvaporsatbar = ones(nw,1)*nan;
fluxFA.r2 = ones(nw,1)*nan;
fluxFA.lagq  = ones(nw,1)*nan;
fluxFA.lagc  = ones(nw,1)*nan;
fluxFA.tbar  = ones(nw,1)*nan;
fluxFA.cbar  = ones(nw,1)*nan;
fluxFA.qbar  = ones(nw,1)*nan;
fluxFA.tt  = ones(nw,1)*nan;
fluxFA.qq  = ones(nw,1)*nan;
fluxFA.cc  = ones(nw,1)*nan; 
fluxFA.wq  = ones(nw,1)*nan;
fluxFA.wc  = ones(nw,1)*nan;
fluxFA.wt2  = ones(nw,1)*nan;
fluxFA.wc2  = ones(nw,1)*nan;
fluxFA.wq2 = ones(nw,1)*nan;
fluxFA.t3 = ones(nw,1)*nan;
fluxFA.q3 = ones(nw,1)*nan;clear all
year=2011;
DATA=[];
hrmin=(0:0.5:23.5)';
YEAR=year*ones(48,1);

%2007 2008 2009 2010 2011
day1=[126 137 138 124 142];
day2=[309 311 295 294 303];

for d=day1(year-2006)+1:day2(year-2006)
    load(['/home/maurer.189/FASET/' num2str(year) 'Proc/Faset_' num2str(year) '_' num2str(d,'%03.0f') '_flux.mat'])
    load(['/home/maurer.189/FASET/' num2str(year) 'Proc/Faset_' num2str(year) '_' num2str(d,'%03.0f') '_data2.mat'])
    
    dectime=d+hrmin/24;
    doy=d*ones(48,1);
    ustar = fluxFA.ustar(:,1);
    tair = fluxFA.tair(:,1);
    WD = fluxFA.WD(:,1);
    ubar = fluxFA.ubar(:,1);
    Htr = fluxFA.Htr(:,1);
    le = fluxFA.LE(:,1);
    humbar = fluxFA.humbar(:,1);
    pbar = fluxFA.pbar(:,1);
    cbar = fluxFA.cbar(:,1);
    qbar = fluxFA.qbar(:,1);
    fc = fluxFA.wcWPL(:,1);
    YEAR=year*ones(48,1);
    
    FA10min.SWin(FA10min.SWin==7999)=nan;
    FA10min.SWout(FA10min.SWout==7999)=nan;
    FA10min.LWin(FA10min.LWin==7999)=nan;
    FA10min.LWout(FA10min.LWout==7999)=nan;
    j=1;
    for i=1:3:142       
        SWin(j,1)=nanmean(FA10min.SWin(i:i+2));
        SWout(j,1)=nanmean(FA10min.SWout(i:i+2));
        LWin(j,1)=nanmean(FA10min.LWin(i:i+2));
        LWout(j,1)=nanmean(FA10min.LWout(i:i+2));
        NetRadbar(j,1)=SWin(j,1)+SWout(j,1)+LWin(j,1)+LWout(j,1);
        j=j+1;
    end
    
        
    DATA=[DATA; YEAR doy hrmin dectime ustar tair WD ubar Htr le humbar pbar cbar qbar NetRadbar SWin SWout LWin LWout fc];
end
fluxFA.c3 = ones(nw,1)*nan;
fluxFA.dtdt = ones(nw,1)*nan;
fluxFA.dqdt = ones(nw,1)*nan;
fluxFA.dcdt = ones(nw,1)*nan;
fluxFA.clag = ones(nw,1)*nan;
fluxFA.NetRadbar = ones(nw,1)*nan; 
fluxFA.PPFDbar = ones(nw,1)*nan; 
fluxFA.SWRadbar = ones(nw,1)*nan;  
fluxFA.rhobar = ones(nw,1)*nan;
fluxFA.rhovbar = ones(nw,1)*nan;
fluxFA.wqWPL = ones(nw,1)*nan;
fluxFA.wcWPL = ones(nw,1)*nan;
fluxFA.spechum = ones(nw,1)*nan;
fluxFA.LE = ones(nw,1)*nan;
fluxFA.wtr= ones(nw,1)*nan; 
fluxFA.Htr = ones(nw,1)*nan; 
fluxFA.WSbar = ones(nw,1)*nan;
fluxFA.WC = ones(nw,1)*nan;
fluxFA.WQ= ones(nw,1)*nan; 
fluxFA.LPM = ones(nw,1)*nan; 
fluxFA.RATc = ones(nw,1)*nan;
fluxFA.RATq = ones(nw,1)*nan;
    
f=3385;%1;
g=142;%1;
h=0;%0; 
currentday=142;%1;

if (year==2000 || year==2004 || year==2008 || year==2012)
   nonleap=0;
else
   nonleap=24;
end
   
while f<=7538;%(length(fastfile(1,:))-nonleap)

    fast = dir([datadirectory char(fastfile(rowID,f))]);    %Checks to see if file slot is present or empty
    
    if isempty(fast)
        [file] = zeros(size(data))*nan;
    else    
        [file] = readdayfile_faset([datadirectory char(fastfile(rowID,f))],year,frec,nhead);   %%%%%%%%%%%%% CHANGE READDAYFILE %%%%%%%%%%%%%
    end

    if f==3385;%1;
        [dataA] = file;
        [dataB] = readdayfile_faset([datadirectory char(fastfile(rowID,f+1))],year,frec,nhead);   %%%%%%%%%%%%% CHANGE READDAYFILE %%%%%%%%%%%%%
        [dataC] = readdayfile_faset([datadirectory char(fastfile(rowID,f+2))],year,frec,nhead);   %%%%%%%%%%%%% CHANGE READDAYFILE %%%%%%%%%%%%%           
        [data2A] = readfaset_slow([datadirectory char(slowfileA(rowID,g))],year);   %%%%%%%%%%%%% CHANGE READDAYFILE %%%%%%%%%%%%% 
        [data2B] = readfaset_slow([datadirectory char(slowfileA(rowID,g+1))],year);   %%%%%%%%%%%%% CHANGE READDAYFILE %%%%%%%%%%%%%
        [data2C] = readfaset_slow([datadirectory char(slowfileA(rowID,g+2))],year);   %%%%%%%%%%%%% CHANGE READDAYFILE %%%%%%%%%%%%%
        
        [data2] = findholesSLOW([data2A; data2B; data2C],year,currentday);
            
        FA10min.year = data2(:,1);
        FA10min.doy = data2(:,2);
        FA10min.dectime = data2(:,3);            
        FA10min.Tm = De_spike_short (144,data2(:,4),4,6,-25,40); %[Celcius]
        FA10min.Hum = De_spike_short (144,data2(:,5),4,6,0,100); %[ % ]
        FA10min.PPFD = data2(:,6); %[micro-mol/(m^2*s)]
        FA10min.lpm = De_spike_short (144,data2(:,7),4,6,5,12); %[lpm]
        FA10min.P = De_spike_short (144,data2(:,8)*100,4,6,90000,120000); %[Pa]
        FA10min.WS = data2(:,9); %[m/s]
        FA10min.WD = data2(:,10);
        FA10min.SWin = data2(:,11); FA10min.SWin(FA10min.SWin==7999)=nan;%[W/m^2] 
        FA10min.SWout = data2(:,12); FA10min.SWout(FA10min.SWout==7999)=nan;%[W/m^2]
        FA10min.NR01_T = data2(:,16); FA10min.NR01_T(FA10min.NR01_T==7999)=nan;
        FA10min.LWin = data2(:,13); FA10min.LWin(FA10min.LWin==7999)=nan;%[W/m^2]
        FA10min.LWin = FA10min.LWin +(567e-10)*(NR01_T+273.15)^4;
        FA10min.LWout = data2(:,14); FA10min.LWout(FA10min.LWout==7999)=nan;%[W/m^2]
        FA10min.LWout = FA10min.LWout +(567e-10)*(NR01_T+273.15)^4;
        FA10min.NetRad = (FA10min.SWin + FA10min.LWin) - (FA10min.SWout + FA10min.LWout); %[W/m^2]             
            
            Pvaporsat = 611.2*exp(17.67*FA10min.Tm./(243.51+FA10min.Tm)); % Saturated Vapor Pressure
            Pvapor = FA10min.Hum.*Pvaporsat./100; % Vapor Pressure
            rho = (FA10min.P./286.9./(FA10min.Tm+273.15)).*(1+0.62197.*(Pvapor./(FA10min.P-Pvapor)))./(1+0.62197.*(Pvapor./(FA10min.P-Pvapor))*461.5/286.9);
            rhov = 0.62197*(Pvapor./(FA10min.P-Pvapor));
        
        f=f+2; %f=f+2 because on the next cycle, we want to read the 4th file, and there is always an f=f+1 at end of cycle
        g=g+2;
        
    elseif any(f==(4:3:length(fastfile(1,:))-2))
        [dataA] = file;
            
    elseif any(f==(5:3:length(fastfile(1,:))-1))
        [dataB] = file;
        
    elseif any(f==(6:3:length(fastfile(1,:))))
        [dataC] = file;
    end
    
    [data dataprelag datapostlag] = findholes3([dataA; dataB; dataC],year,currentday,h,frec,preclevel);
    

    if any(f==(3411:24:length(fastfile(1,:))))%(27:24:length(fastfile(1,:))))
            
            save([procdir site '_' num2str(year) '_' num2str(currentday-1,'%03.0f') '_' 'data2.mat'],'FA10min');
            save([procdir site '_' num2str(year) '_' num2str(currentday-1,'%03.0f') '_' 'flux.mat'],'fluxFA');
   
            fluxFA.ubar = ones(nw,1)*nan;
            fluxFA.wbar= ones(nw,1)*nan;
            fluxFA.tsvbar= ones(nw,1)*nan;
            fluxFA.uu= ones(nw,1)*nan;
            fluxFA.vv= ones(nw,1)*nan;
            fluxFA.ww= ones(nw,1)*nan;
            fluxFA.uw= ones(nw,1)*nan;
            fluxFA.vw= ones(nw,1)*nan;
            fluxFA.u3= ones(nw,1)*nan;
            fluxFA.v3= ones(nw,1)*nan;
            fluxFA.w3= ones(nw,1)*nan;
            fluxFA.ustar= ones(nw,1)*nan;
            fluxFA.L= ones(nw,1)*nan;
            fluxFA.WD = ones(nw,1)*nan;
            fluxFA.tair = ones(nw,1)*nan;
            fluxFA.tvair  = ones(nw,1)*nan;
            fluxFA.pbar  = ones(nw,1)*nan;
            fluxFA.humbar = ones(nw,1)*nan;
            fluxFA.Pvaporbar = ones(nw,1)*nan;
            fluxFA.Pvaporsatbar = ones(nw,1)*nan;
            fluxFA.r2 = ones(nw,1)*nan;
            fluxFA.lagq  = ones(nw,1)*nan;
            fluxFA.lagc  = ones(nw,1)*nan;
            fluxFA.tbar  = ones(nw,1)*nan;
            fluxFA.cbar  = ones(nw,1)*nan;
            fluxFA.qbar  = ones(nw,1)*nan;
            fluxFA.tt  = ones(nw,1)*nan;
            fluxFA.qq  = ones(nw,1)*nan;
            fluxFA.cc  = ones(nw,1)*nan;
            fluxFA.wq  = ones(nw,1)*nan;
            fluxFA.wc  = ones(nw,1)*nan;
            fluxFA.wt2  = ones(nw,1)*nan;
            fluxFA.wc2  = ones(nw,1)*nan;
            fluxFA.wq2 = ones(nw,1)*nan;
            fluxFA.t3 = ones(nw,1)*nan;
            fluxFA.q3 = ones(nw,1)*nan;
            fluxFA.c3 = ones(nw,1)*nan;
            fluxFA.dtdt = ones(nw,1)*nan;
            fluxFA.dqdt = ones(nw,1)*nan;
            fluxFA.dcdt = ones(nw,1)*nan;
            fluxFA.clag = ones(nw,1)*nan;
            fluxFA.NetRadbar = ones(nw,1)*nan; 
            fluxFA.PPFDbar = ones(nw,1)*nan; 
            fluxFA.SWRadbar = ones(nw,1)*nan;  
            fluxFA.rhobar = ones(nw,1)*nan; 
            fluxFA.rhovbar = ones(nw,1)*nan;
            fluxFA.wqWPL = ones(nw,1)*nan;
            fluxFA.wcWPL = ones(nw,1)*nan;
            fluxFA.spechum = ones(nw,1)*nan;
            fluxFA.LE = ones(nw,1)*nan;
            fluxFA.wtr= ones(nw,1)*nan; 
            fluxFA.Htr = ones(nw,1)*nan;
            fluxFA.WSbar = ones(nw,1)*nan;
            fluxFA.WC = ones(nw,1)*nan;
            fluxFA.WQ= ones(nw,1)*nan; 
            fluxFA.LPM = ones(nw,1)*nan; 
            fluxFA.RATc = ones(nw,1)*nan;
            fluxFA.RATq = ones(nw,1)*nan;
  
            if floor(data(12000,3))~=h && floor(data(12000,3))~=23
            h = floor(data(12000,3));
            end
            
                  slow = dir([datadirectory char(slowfileA(rowID,g))]); %Checks to see if m46/m34 file slot is present or empty
                  
                  if ~isempty(slow)
                      if any(g==(4:3:length(slowfileA(1,:))-2))
                          [data2A] = readfaset_slow([datadirectory char(slowfileA(rowID,g))],year);   %%%%%%%%%%%%% CHANGE READDAYFILE %%%%%%%%%%%%% 
                      elseif any(g==(5:3:length(slowfileA(1,:))-1))
                          [data2B] = readfaset_slow([datadirectory char(slowfileA(rowID,g))],year);   %%%%%%%%%%%%% CHANGE READDAYFILE %%%%%%%%%%%%%
                      elseif any(g==(6:3:length(slowfileA(1,:))))
                          [data2C] = readfaset_slow([datadirectory char(slowfileA(rowID,g))],year);   %%%%%%%%%%%%% CHANGE READDAYFILE %%%%%%%%%%%%%
                      end
                  else
                      if any(g==(4:3:length(slowfileA(1,:))-2))
                          [data2A] = [year*ones(144,1) currentday*ones(144,1) nan*ones(size(data2(:,3:end)))];
                      elseif any(g==(5:3:length(slowfileA(1,:))-1))
                          [data2B] = [year*ones(144,1) currentday*ones(144,1) nan*ones(size(data2(:,3:end)))];
                      elseif any(g==(6:3:length(slowfileA(1,:))))
                          [data2C] = [year*ones(144,1) currentday*ones(144,1) nan*ones(size(data2(:,3:end)))];
                      end
                  end
                  
                  rows2 = min([size(data2A,2) size(data2B,2) size(data2C,2)]);
                  [data2] = findholesSLOW([data2A(:,1:rows2); data2B(:,1:rows2); data2C(:,1:rows2)],year,currentday);
                        
                    FA10min.year = data2(:,1);
                    FA10min.doy = data2(:,2);
                    FA10min.dectime = data2(:,3);            
                    FA10min.Tm = De_spike_short (144,data2(:,4),4,6,-25,40); %[Celcius]
                    FA10min.Hum = De_spike_short (144,data2(:,5),4,6,0,100); %[ % ]
                    FA10min.PPFD = data2(:,6); %[micro-mol/(m^2*s)]                    
                    FA10min.lpm = De_spike_short (144,data2(:,7),4,6,5,12); %[lpm]
                    FA10min.P = De_spike_short (144,data2(:,8)*100,4,6,90000,120000); %[Pa]
                    FA10min.WS = data2(:,9); %[m/s]
                    FA10min.WD = data2(:,10);
                    FA10min.SWin = data2(:,11); FA10min.SWin(FA10min.SWin==7999)=nan;%[W/m^2] 
                    FA10min.SWout = data2(:,12); FA10min.SWout(FA10min.SWout==7999)=nan;%[W/m^2]
                    FA10min.NR01_T = data2(:,16); FA10min.NR01_T(FA10min.NR01_T==7999)=nan;
                    FA10min.LWin = data2(:,13); FA10min.LWin(FA10min.LWin==7999)=nan;%[W/m^2]
                    FA10min.LWin = FA10min.LWin +(567e-10)*(NR01_T+273.15)^4;
                    FA10min.LWout = data2(:,14); FA10min.LWout(FA10min.LWout==7999)=nan;%[W/m^2]
                    FA10min.LWout = FA10min.LWout +(567e-10)*(NR01_T+273.15)^4;
                    FA10min.NetRad = (FA10min.SWin + FA10min.LWin) - (FA10min.SWout + FA10min.LWout); %[W/m^2]                     

                    Pvaporsat = 611.2*exp(17.67*FA10min.Tm./(243.51+FA10min.Tm)); % Saturated Vapor Pressure
                    Pvapor = FA10min.Hum.*Pvaporsat./100; % Vapor Pressure
                    rho = (FA10min.P./286.9./(FA10min.Tm+273.15)).*(1+0.62197.*(Pvapor./(FA10min.P-Pvapor)))./(1+0.62197.*(Pvapor./(FA10min.P-Pvapor))*461.5/286.9); %air density [kg/m3]
                    rhov = 0.62197*(Pvapor./(FA10min.P-Pvapor)); %vapor density [kg/m3]
                        
                  
    end
    
    [ts c]=size(data);
    
    if ts < 0.6*ntsf
        disp(['there ARE ONLY ' num2str(ts) ' timesteps. cannot work like that'])
        data(ts+1:nstf,:) = nan; %if 2nd half of file is too short, program will crash. If that is the case, add nan to allow run.
    end      
    
    DespikeFA.TIME = data(:,1:3);
    DespikeFA.U = De_spike3 (ts,data(:,4),1200,6,-25,25,'CSAT3',data(:,8)); %[m/s]
    DespikeFA.V = De_spike3 (ts,data(:,5),1200,6,-25,25,'CSAT3',data(:,8)); %[m/s]
    DespikeFA.W = De_spike3 (ts,data(:,6),1200,6,-10,10,'CSAT3',data(:,8)); %[m/s]
    DespikeFA.T = De_spike3 (ts,data(:,7),3000,4,-25,55,'CSAT3',data(:,8)); %[Celsius] 
    DespikeFA.C = De_spike3 (ts,data(:,9),1200,6,340,460,'NONE',ones(size(data(:,9))),1.5,1); %[umol/mol]
    DespikeFA.Q = De_spike3 (ts,data(:,10),1200,6,0,30,'NONE',ones(size(data(:,9))),1.5); %[mmol/mol]   
    DespikeFA.IRGA_T = De_spike3 (ts,data(:,11),1200,6,-25,40); %[Celsius]
    DespikeFA.IRGA_P = De_spike3 (ts,data(:,12),1200,6,50000,120000); %[Pascals]    
    
    r = DespikeFA.Q(1:ts)/1000*18/28.97; % mixing ratio
    Tv = DespikeFA.T(1:ts).*(1+0.61*r); % Virtual Temp 
  
    Nmin30 = floor(ntsf/frec/1800); %Number of half hours per file
    
   
    for b=1:Nmin30
    lastuse=min(b*nn,ts);
        use=((b-1)*nn+1:lastuse)';
        currenthour = h;
        good=find(isnan(DespikeFA.T(use))==0 & isnan(DespikeFA.W(use))==0);        
        
        if b==Nmin30 && h==23
            fluxFA.Pvaporbar((currenthour*Nmin30)+b) = nanmean(Pvapor(1+((currenthour*Nmin30)+b-1)*3:end));
            fluxFA.Pvaporsatbar((currenthour*Nmin30)+b) = nanmean(Pvaporsat(1+((currenthour*Nmin30)+b-1)*3:end));
            fluxFA.pbar((currenthour*Nmin30)+b) = nanmean(FA10min.P(1+((currenthour*Nmin30)+b-1)*3:end));
            fluxFA.r2((currenthour*Nmin30)+b) = 0.622*Pvapor((currenthour*Nmin30)+b)/(fluxFA.pbar((currenthour*Nmin30)+b)-Pvapor((currenthour*Nmin30)+b));
            fluxFA.tair((currenthour*Nmin30)+b) = nanmean(FA10min.Tm(1+((currenthour*Nmin30)+b-1)*3:end));
            fluxFA.tvair((currenthour*Nmin30)+b) = fluxFA.tair((currenthour*Nmin30)+b)*(1+0.61*fluxFA.r2((currenthour*Nmin30)+b));
            fluxFA.humbar((currenthour*Nmin30)+b) = nanmean(FA10min.Hum(1+((currenthour*Nmin30)+b-1)*3:end));
            fluxFA.NetRadbar((currenthour*Nmin30)+b) = nanmean(FA10min.NetRad(1+((currenthour*Nmin30)+b-1)*3:end));
            fluxFA.PPFDbar((currenthour*Nmin30)+b) = nanmean(FA10min.PPFD(1+((currenthour*Nmin30)+b-1)*3:end));            
            fluxFA.rhobar((currenthour*Nmin30)+b) = nanmean(rho(1+((currenthour*Nmin30)+b-1)*3:end));
            fluxFA.rhovbar((currenthour*Nmin30)+b) = nanmean(rhov(1+((currenthour*Nmin30)+b-1)*3:end));
            fluxFA.WSbar((currenthour*Nmin30)+b) = nanmean(FA10min.WS(1+((currenthour*Nmin30)+b-1)*3:end));
            fluxFA.LPM((currenthour*Nmin30)+b) = nanmean(FA10min.lpm(1+((currenthour*Nmin30)+b-1)*3:end));
            
                else fluxFA.Pvaporbar((currenthour*Nmin30)+b) = nanmean(Pvapor(1+((currenthour*Nmin30)+b-1)*3:3+((currenthour*Nmin30)+b-1)*3));
                     fluxFA.Pvaporsatbar((currenthour*Nmin30)+b) = nanmean(Pvaporsat(1+((currenthour*Nmin30)+b-1)*3:3+((currenthour*Nmin30)+b-1)*3));
                     fluxFA.pbar((currenthour*Nmin30)+b) = nanmean(FA10min.P(1+((currenthour*Nmin30)+b-1)*3:3+((currenthour*Nmin30)+b-1)*3));
                     fluxFA.r2((currenthour*Nmin30)+b) = 0.622*Pvapor((currenthour*Nmin30)+b)/(fluxFA.pbar((currenthour*Nmin30)+b)-Pvapor((currenthour*Nmin30)+b));
                     fluxFA.tair((currenthour*Nmin30)+b) = nanmean(FA10min.Tm(1+((currenthour*Nmin30)+b-1)*3:3+((currenthour*Nmin30)+b-1)*3));
                     fluxFA.tvair((currenthour*Nmin30)+b) = fluxFA.tair((currenthour*Nmin30)+b)*(1+0.61*fluxFA.r2((currenthour*Nmin30)+b));
                     fluxFA.humbar((currenthour*Nmin30)+b) = nanmean(FA10min.Hum(1+((currenthour*Nmin30)+b-1)*3:3+((currenthour*Nmin30)+b-1)*3));
                     fluxFA.NetRadbar((currenthour*Nmin30)+b) = nanmean(FA10min.NetRad(1+((currenthour*Nmin30)+b-1)*3:3+((currenthour*Nmin30)+b-1)*3));
                     fluxFA.PPFDbar((currenthour*Nmin30)+b) = nanmean(FA10min.PPFD(1+((currenthour*Nmin30)+b-1)*3:3+((currenthour*Nmin30)+b-1)*3));                    
                     fluxFA.rhobar((currenthour*Nmin30)+b) = nanmean(rho(1+((currenthour*Nmin30)+b-1)*3:3+((currenthour*Nmin30)+b-1)*3));
                     fluxFA.rhovbar((currenthour*Nmin30)+b) = nanmean(rhov(1+((currenthour*Nmin30)+b-1)*3:3+((currenthour*Nmin30)+b-1)*3));
                     fluxFA.WSbar((currenthour*Nmin30)+b) = nanmean(FA10min.WS(1+((currenthour*Nmin30)+b-1)*3:3+((currenthour*Nmin30)+b-1)*3));
                     fluxFA.LPM((currenthour*Nmin30)+b) = nanmean(FA10min.lpm(1+((currenthour*Nmin30)+b-1)*3:3+((currenthour*Nmin30)+b-1)*3));
        end
        if length(good)>0.5*nn
        %-- 2D-rotation the coordinate system in the window----------------
          theta=atan2(nanmean(DespikeFA.V(use)),nanmean(DespikeFA.U(use)));
          alfa=atan2(-nanmean(DespikeFA.W(use)),(nanmean(DespikeFA.U(use))^2+nanmean(DespikeFA.V(use))^2)^0.5);
          Az=[cos(theta) sin(theta) 0; -sin(theta) cos(theta) 0; 0 0 1];
          Ze=[cos(alfa) 0 -sin(alfa); 0 1 0; sin(alfa) 0 cos(alfa)];
          R=Ze*Az;
          u=[DespikeFA.U(use) DespikeFA.V(use) DespikeFA.W(use)]*R(1,:)';
          v=[DespikeFA.U(use) DespikeFA.V(use) DespikeFA.W(use)]*R(2,:)';
          w=[DespikeFA.U(use) DespikeFA.V(use) DespikeFA.W(use)]*R(3,:)';

          %--means
          fluxFA.ubar((currenthour*Nmin30)+b)=nanmean(u);
          fluxFA.wbar((currenthour*Nmin30)+b)=nanmean(DespikeFA.W(use));
          fluxFA.tsvbar((currenthour*Nmin30)+b)=nanmean(Tv(use));
          
          
          un=u-fluxFA.ubar((currenthour*Nmin30)+b);

          fluxFA.uu((currenthour*Nmin30)+b)=nanvar(u);
          fluxFA.vv((currenthour*Nmin30)+b)=nanvar(v);
          fluxFA.ww((currenthour*Nmin30)+b)=nanvar(w);
          fluxFA.u3((currenthour*Nmin30)+b)=moment(u,3);
          fluxFA.v3((currenthour*Nmin30)+b)=moment(v,3);
          fluxFA.w3((currenthour*Nmin30)+b)=moment(w,3);
          fluxFA.uw((currenthour*Nmin30)+b)=nanmean(w.*un);
          fluxFA.vw((currenthour*Nmin30)+b)=nanmean(w.*v);
            
                    %other variables
          fluxFA.ustar((currenthour*Nmin30)+b)=(fluxFA.uw((currenthour*Nmin30)+b)^2+fluxFA.vw((currenthour*Nmin30)+b)^2)^0.25;
          fluxFA.wts=nanmean((Tv(use)-fluxFA.tsvbar((currenthour*Nmin30)+b)).*w);
          fluxFA.L((currenthour*Nmin30)+b)=-(fluxFA.ustar((currenthour*Nmin30)+b)^3)/(0.4*9.81*(fluxFA.wts/(fluxFA.tvair((currenthour*Nmin30)+b)+273.15)));
          fluxFA.WD((currenthour*Nmin30)+b)=theta;
          
        elseif length(good)<=0.5*nn
            theta=nan; alfa=nan; Az=nan; Ze=nan; R=nan; u=nan*ones(size(DespikeFA.U(use))); v=nan*ones(size(DespikeFA.V(use))); w=nan*ones(size(DespikeFA.W(use)));

        end
        
         %--SCALARS--------------------------------
    
         %H2O and CO2
         good_q=find(isnan(DespikeFA.Q(use))==0 & isnan(DespikeFA.C(use))==0);
         if length(good_q)>0.6*nn

            % find lag time between Li-7500 and sonic
            fluxFA.lagq((currenthour*Nmin30)+b)=lag_test(w(good_q),DespikeFA.Q(use(good_q)),IRGA);
            fluxFA.lagc((currenthour*Nmin30)+b)=lag_test(w(good_q),DespikeFA.C(use(good_q)),IRGA);
            lagall = nanmin(fluxFA.lagq((currenthour*Nmin30)+b),fluxFA.lagc((currenthour*Nmin30)+b));
            if lagall>0;
                 C1 = De_spike3 (ts,[data(lagall+1:end,9); datapostlag(1:lagall,9)],1200,6,340,460,'NONE',ones(size(data(:,9))),1.5,1);
                 Q1 = De_spike3 (ts,[data(lagall+1:end,10); datapostlag(1:lagall,10)],1200,6,0,30,'NONE',ones(size(data(:,9))),1.5);
            elseif lagall<0
                C1 = De_spike3 (ts,[dataprelag((end+lagall+1):end,9); data(1:(end+lagall),9)],1200,6,340,460,'NONE',ones(size(data(:,9))),1.5,1);
                Q1 = De_spike3 (ts,[dataprelag((end+lagall+1):end,10); data(1:end+lagall,10)],1200,6,0,30,'NONE',ones(size(data(:,9))),1.5);
            else
                C1=DespikeFA.C;
                Q1=DespikeFA.Q;
            end

         %--apply WPL correction on raw data                    
            
            t=DespikeFA.T(use);
            q=Q1(use);
            c=C1(use);
            
            % evaluete statistics
            fluxFA.clag((currenthour*Nmin30)+b)=lagall;
            fluxFA.tbar((currenthour*Nmin30)+b)=nanmean(DespikeFA.T(use));
            fluxFA.qbar((currenthour*Nmin30)+b)=nanmean(Q1(use));
            fluxFA.cbar((currenthour*Nmin30)+b)=nanmean(C1(use));

            qn=Q1(use)-fluxFA.qbar((currenthour*Nmin30)+b);
            cn=C1(use)-fluxFA.cbar((currenthour*Nmin30)+b);

            fluxFA.tt((currenthour*Nmin30)+b)=nanvar(t);
            fluxFA.qq((currenthour*Nmin30)+b)=nanvar(q);
            fluxFA.cc((currenthour*Nmin30)+b)=nanvar(c);

            fluxFA.wq((currenthour*Nmin30)+b)=nanmean(w.*qn); 
            fluxFA.wc((currenthour*Nmin30)+b)=nanmean(w.*cn);

            if ~isnan(fluxFA.tair((currenthour*Nmin30)+b))
                tairE=fluxFA.tair((currenthour*Nmin30)+b);
                fluxFA.spechum((currenthour*Nmin30)+b) = fluxFA.rhovbar((currenthour*Nmin30)+b)/fluxFA.rhobar((currenthour*Nmin30)+b);
            else
                tairE=nanmean(DespikeFA.T(use));
                fluxFA.spechum((currenthour*Nmin30)+b) = nan;
            end  
            
            [fluxFA.wtr((currenthour*Nmin30)+b),fluxFA.Htr((currenthour*Nmin30)+b),Tr] = KaimalGaynor1990_2(fluxFA.Pvaporbar((currenthour*Nmin30)+b),fluxFA.pbar((currenthour*Nmin30)+b),Q1(use),DespikeFA.T(use),fluxFA.Pvaporsatbar((currenthour*Nmin30)+b),w,fluxFA.tair((currenthour*Nmin30)+b),IRGA);
                             
            [fluxFA.wqWPL((currenthour*Nmin30)+b),fluxFA.wcWPL((currenthour*Nmin30)+b),fluxFA.LE((currenthour*Nmin30)+b)]=raw_WPL_correction(fluxFA.wtr((currenthour*Nmin30)+b),fluxFA.wq((currenthour*Nmin30)+b),fluxFA.wc((currenthour*Nmin30)+b),Q1(use),C1(use),DespikeFA.IRGA_P(use),DespikeFA.IRGA_T(use),fluxFA.pbar((currenthour*Nmin30)+b),tairE,fluxFA.spechum((currenthour*Nmin30)+b),IRGA);
            
            [fluxFA.RATc((currenthour*Nmin30)+b)]=cpirgafreq(fluxFA.ubar((currenthour*Nmin30)+b),lso,z_geom/fluxFA.L((currenthour*Nmin30)+b),AVGtime,fluxFA.tair((currenthour*Nmin30)+b),fluxFA.LPM((currenthour*Nmin30)+b),TubeD,1,1,LT,z_geom);
            [fluxFA.RATq((currenthour*Nmin30)+b)]=cpirgafreq(fluxFA.ubar((currenthour*Nmin30)+b),lso,z_geom/fluxFA.L((currenthour*Nmin30)+b),AVGtime,fluxFA.tair((currenthour*Nmin30)+b),fluxFA.LPM((currenthour*Nmin30)+b),TubeD,2,1,LT,z_geom);
            
            fluxFA.WC((currenthour*Nmin30)+b)=fluxFA.wcWPL((currenthour*Nmin30)+b)/fluxFA.RATc((currenthour*Nmin30)+b);
            fluxFA.WQ((currenthour*Nmin30)+b)=fluxFA.wqWPL((currenthour*Nmin30)+b)/fluxFA.RATq((currenthour*Nmin30)+b);
            
            % compute fluxes as J-block average for stationary test
             J=2;
             for j=1:J
                jjj=nn/J*(j-1)+1:nn/J*j;
                fluxFA.wt2((currenthour*Nmin30)+b)=fluxFA.wt2((currenthour*Nmin30)+b)+(nanmean(w(jjj).*t(jjj))-nanmean(w(jjj)).*nanmean(t(jjj)))/J;  
                fluxFA.wq2((currenthour*Nmin30)+b)=fluxFA.wq2((currenthour*Nmin30)+b)+(nanmean(w(jjj).*q(jjj))-nanmean(w(jjj)).*nanmean(q(jjj)))/J; 
                fluxFA.wc2((currenthour*Nmin30)+b)=fluxFA.wc2((currenthour*Nmin30)+b)+(nanmean(w(jjj).*c(jjj))-nanmean(w(jjj)).*nanmean(c(jjj)))/J;  
             end

            fluxFA.t3((currenthour*Nmin30)+b)=moment(t,3);
            fluxFA.q3((currenthour*Nmin30)+b)=moment(q,3);
            fluxFA.c3((currenthour*Nmin30)+b)=moment(c,3);

         %storage terms - dx is averaged on a sub-interval (2D [sec]) at the beginning and the
         %end of the interval window
            D=120/dt;
            if b>1 && b<Nmin30   
                fluxFA.dtdt((currenthour*Nmin30)+b)=(nanmean(DespikeFA.T(use(end)-D:use(end)+D))-nanmean(DespikeFA.T(use(1)-D:use(1)+D)))/DT;
                fluxFA.dqdt((currenthour*Nmin30)+b)=(nanmean(DespikeFA.Q(use(end)-D:use(end)+D))-nanmean(DespikeFA.Q(use(1)-D:use(1)+D)))/DT;
                fluxFA.dcdt((currenthour*Nmin30)+b)=(nanmean(DespikeFA.C(use(end)-D:use(end)+D))-nanmean(DespikeFA.C(use(1)-D:use(1)+D)))/DT;

            elseif b==1
                fluxFA.dtdt((currenthour*Nmin30)+b)=(nanmean(DespikeFA.T(use(end)-D:use(end)+D))-nanmean(DespikeFA.T(use(1):use(1)+D)))/DT;
                fluxFA.dqdt((currenthour*Nmin30)+b)=(nanmean(DespikeFA.Q(use(end)-D:use(end)+D))-nanmean(DespikeFA.Q(use(1):use(1)+D)))/DT;
                fluxFA.dcdt((currenthour*Nmin30)+b)=(nanmean(DespikeFA.C(use(end)-D:use(end)+D))-nanmean(DespikeFA.C(use(1):use(1)+D)))/DT;

            elseif b==Nmin30
                fluxFA.dtdt((currenthour*Nmin30)+b)=(nanmean(DespikeFA.T(use(end)-D:use(end)))-nanmean(DespikeFA.T(use(1)-D:use(1)+D)))/DT;
                fluxFA.dqdt((currenthour*Nmin30)+b)=(nanmean(DespikeFA.Q(use(end)-D:use(end)))-nanmean(DespikeFA.Q(use(1)-D:use(1)+D)))/DT;
                fluxFA.dcdt((currenthour*Nmin30)+b)=(nanmean(DespikeFA.C(use(end)-D:use(end)))-nanmean(DespikeFA.C(use(1)-D:use(1)+D)))/DT;

            end
         end %if good_q
    end %for b
    
    
    
    save([procdir site '_' num2str(year) '_' num2str(currentday,'%03.0f') '_' num2str(h,'%02.0f') '_data.mat'],'data');
    save([procdir site '_' num2str(year) '_' num2str(currentday,'%03.0f') '_' num2str(h,'%02.0f') '_despike.mat'],'DespikeFA');
    
if f==7538 %(length(fastfile(1,:))-nonleap)
            save([procdir site '_' num2str(year) '_' num2str(currentday,'%03.0f') '_' 'data2.mat'],'FA10min');
            save([procdir site '_' num2str(year) '_' num2str(currentday,'%03.0f') '_flux.mat'],'fluxFA');          
end
     if any(f==(3410:24:length(fastfile(1,:))))%(26:24:length(fastfile(1,:)))) %with fastfile list, currentday changes are known, change currentday on those marks
         currentday = currentday+1; disp(currentday);
         g=g+1;
         h=(-1); % turn h into -1 here so that h+1 at the start of a new day = 0
     end
     f=f+1;
     h=h+1;

end %for f