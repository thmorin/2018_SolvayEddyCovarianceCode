% function [data]=get_Fc_coponents(hh, mm, Doy, year, str)
% 
% 
% load fluxW.mat
% load DespikeW
% load Fdata
% load Sdata
% clear all

%%%%%%%%%%%User Defined Variables%%%%%%%%%%%%%%
Fc_Spike_Days=[36.9 82.45 284.38 215.1 206.58 203.82 230.9 231.8 237];
M_Spike_Days=[83.15 93.32 96.1 133.3];

M_Spike_of_interest=1;
Fc_Spike_of_interest=1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

directory='/home/morin.37/poolA/wetland_data/Matlab_Liel/Processed2011/';



%%%%%%%%%%%Run for Carbon spike data%%%%%%%%%%%%%
ApproxDay=Fc_Spike_Days(Fc_Spike_of_interest);

Day=floor(ApproxDay);
HlfHr=round(ApproxDay-Day)*48;
Despikedaylines=24*60*60*10;
W1mindaylines=24*60;

RelevDays=(Day-1):(Day+1);

Fc=nan(48*3,1);
Fc2=nan(Despikedaylines*3,1);
DOY2=nan(Despikedaylines*3,1);
DOY=nan(48*3,1);
RMYW=nan(Despikedaylines*3,1);
CSATW=nan(Despikedaylines*3,1);

numyplots=3;
numxplots=2;

for CD=RelevDays
    load([directory 'DespikeW/OWR_2011_' num2str(CD,'%03.0f') '_DespikeW.mat']);
    Fc2(((CD-RelevDays(1))*Despikedaylines+1):((CD-RelevDays(1))*Despikedaylines+Despikedaylines))=DespikeW.LI7500_C;
    RMYW(((CD-RelevDays(1))*Despikedaylines+1):((CD-RelevDays(1))*Despikedaylines+Despikedaylines))=DespikeW.RMY_W;
    CSATW(((CD-RelevDays(1))*Despikedaylines+1):((CD-RelevDays(1))*Despikedaylines+Despikedaylines))=DespikeW.CSAT_W;
    DOY2(((CD-RelevDays(1))*Despikedaylines+1):((CD-RelevDays(1))*Despikedaylines+Despikedaylines))=CD:1/(Despikedaylines):CD+(Despikedaylines-1)/(Despikedaylines);
    Tair(((CD-RelevDays(1))*Despikedaylines+1):((CD-RelevDays(1))*Despikedaylines+Despikedaylines))=DespikeW.CSAT_Tmp;
    LI7500_P(((CD-RelevDays(1))*Despikedaylines+1):((CD-RelevDays(1))*Despikedaylines+Despikedaylines))=DespikeW.LI7500_P;
    LI7700_P(((CD-RelevDays(1))*Despikedaylines+1):((CD-RelevDays(1))*Despikedaylines+Despikedaylines))=DespikeW.LI7700_P;
    LI7500_Q(((CD-RelevDays(1))*Despikedaylines+1):((CD-RelevDays(1))*Despikedaylines+Despikedaylines))=DespikeW.LI7500_Q;


    load([directory 'flux/OWR_2011_' num2str(CD,'%03.0f') '_flux.mat']);
    Fc(((CD-RelevDays(1))*48+1):((CD-RelevDays(1))*48+48))=fluxW.Fc;
    DOY(((CD-RelevDays(1))*48+1):((CD-RelevDays(1))*48+48))=CD:1/48:CD+47/48;



    load([directory 'W1min/OWR_2011_' num2str(CD,'%03.0f') '_W1min.mat']);
    Tm(((CD-RelevDays(1))*W1mindaylines+1):((CD-RelevDays(1))*W1mindaylines+W1mindaylines))=W1min.Tm;
    DOY3(((CD-RelevDays(1))*W1mindaylines+1):((CD-RelevDays(1))*W1mindaylines+W1mindaylines))=CD:1/(W1mindaylines):CD+(W1mindaylines-1)/(W1mindaylines);
end

figure(6);
subplot(numyplots,numxplots,1)
plotyy(DOY,Fc,DOY2,Fc2);
grid on;
xlabel('Day')
legend('Half hour carbon flux','10 Hz Carbon density');

subplot(numyplots,numxplots,2)
ax=plotyy(DOY,Fc,DOY2,RMYW);
grid on
hold(ax(1),'on');
hold(ax(2),'on');
plot(DOY2,CSATW,'g','Parent',ax(2));
hold(ax(1),'off');
hold(ax(2),'off');
xlabel('Day');
legend('Half hour carbon flux','RMY_W','CSAT_W');

subplot(numyplots,numxplots,3)
ax=plotyy(DOY,Fc,DOY2,Tair);
hold(ax(1),'on');
hold(ax(2),'on');
plot(DOY3,Tm,'g','Parent',ax(2));
xlabel('Day')
hold(ax(1),'off');
hold(ax(2),'off');
legend('Half hour carbon flux','CSAT Tmp','W1min Tm');

subplot(numyplots,numxplots,4)
ax=plotyy(DOY,Fc,DOY2,LI7500_P);
hold(ax(1),'on');
hold(ax(2),'on');
plot(DOY2,LI7700_P,'g','Parent',ax(2));
xlabel('Day')
hold(ax(1),'off');
hold(ax(2),'off');
legend('Half hour carbon flux','LI7500 P','LI7700 P');

subplot(numyplots,numxplots,5)
plotyy(DOY,Fc,DOY2,LI7500_Q);
grid on;
xlabel('Day')
legend('Half hour carbon flux','10 Hz Water Density');
%%%%%%%%%%End carbon spike data analysis%%%%%%%%%%%%%%


%%%%%%%%%%Begin Methane spike data analysis%%%%%%%%%%%
ApproxDay=M_Spike_Days(M_Spike_of_interest);

Day=floor(ApproxDay);
HlfHr=round(ApproxDay-Day)*48;
Despikedaylines=24*60*60*10;
W1mindaylines=24*60;

RelevDays=(Day-1):(Day+1);

M=nan(48*3,1);
LI7700_M=nan(Despikedaylines*3,1);
DOY2=nan(Despikedaylines*3,1);
DOY=nan(48*3,1);
RMYW=nan(Despikedaylines*3,1);
CSATW=nan(Despikedaylines*3,1);

numyplots=1;
numxplots=1;

for CD=RelevDays
    load([directory 'DespikeW/OWR_2011_' num2str(CD,'%03.0f') '_DespikeW.mat']);
    LI7700_M(((CD-RelevDays(1))*Despikedaylines+1):((CD-RelevDays(1))*Despikedaylines+Despikedaylines))=DespikeW.LI7700_M;
    RMYW(((CD-RelevDays(1))*Despikedaylines+1):((CD-RelevDays(1))*Despikedaylines+Despikedaylines))=DespikeW.RMY_W;
    CSATW(((CD-RelevDays(1))*Despikedaylines+1):((CD-RelevDays(1))*Despikedaylines+Despikedaylines))=DespikeW.CSAT_W;
    DOY2(((CD-RelevDays(1))*Despikedaylines+1):((CD-RelevDays(1))*Despikedaylines+Despikedaylines))=CD:1/(Despikedaylines):CD+(Despikedaylines-1)/(Despikedaylines);
    Tair(((CD-RelevDays(1))*Despikedaylines+1):((CD-RelevDays(1))*Despikedaylines+Despikedaylines))=DespikeW.CSAT_Tmp;
    LI7500_P(((CD-RelevDays(1))*Despikedaylines+1):((CD-RelevDays(1))*Despikedaylines+Despikedaylines))=DespikeW.LI7500_P;
    LI7700_P(((CD-RelevDays(1))*Despikedaylines+1):((CD-RelevDays(1))*Despikedaylines+Despikedaylines))=DespikeW.LI7700_P;
    LI7500_Q(((CD-RelevDays(1))*Despikedaylines+1):((CD-RelevDays(1))*Despikedaylines+Despikedaylines))=DespikeW.LI7500_Q;


    load([directory 'flux/OWR_2011_' num2str(CD,'%03.0f') '_flux.mat']);
    M(((CD-RelevDays(1))*48+1):((CD-RelevDays(1))*48+48))=fluxW.M;
    DOY(((CD-RelevDays(1))*48+1):((CD-RelevDays(1))*48+48))=CD:1/48:CD+47/48;



    load([directory 'W1min/OWR_2011_' num2str(CD,'%03.0f') '_W1min.mat']);
    Tm(((CD-RelevDays(1))*W1mindaylines+1):((CD-RelevDays(1))*W1mindaylines+W1mindaylines))=W1min.Tm;
    DOY3(((CD-RelevDays(1))*W1mindaylines+1):((CD-RelevDays(1))*W1mindaylines+W1mindaylines))=CD:1/(W1mindaylines):CD+(W1mindaylines-1)/(W1mindaylines);
end

figure(6);
subplot(numyplots,numxplots,1)
plotyy(DOY,M,DOY2,LI7700_M);
grid on;
xlabel('Day')
legend('Half hour methane flux','10 Hz Methane density');
%%%%%%%%%%Run for methane spike analysis%%%%%%%%%%%%%%