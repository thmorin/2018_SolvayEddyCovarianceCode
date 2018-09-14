clear all
year=2011;

if year==2008 || year==2012
    numdays=366;
elseif year==2010
    numdays=316;
else
    numdays=365;
end

ustar=[];Lo=[];sv=[];windir=[];
for d=1:142
    
    load(['/home/maurer.189/FASET/' num2str(year) 'Proc/Faset_' num2str(year) '_' num2str(d,'%03.0f') '_flux.mat'])
    
    ustar=[ustar;fluxFA.ustar];
    Lo=[Lo;fluxFA.L];
    sv=[sv;sqrt(fluxFA.ww)];
    windir=[windir;fluxFA.WD];
   
    clear fluxFA
end

addpath('/old-drive/Processdata/gapfill/')

[x ustar]=De_spike_gf(ustar,(numdays)*48,6,48);
[x Lo]=De_spike_gf(Lo,(numdays)*48,6,48);
[x sv]=De_spike_gf(sv,(numdays)*48,6,48);
[x windir]=De_spike_gf(windir,(numdays)*48,6,48);

zo=0.86; zH=34; d=25.70; plotYN=0;

load('FASET.mat')

FOOTSUM=[];
for runs=1:length(ustar)
    
    [foot,xx,yy,footRotate,Fxprime,Fyprime,footsum]=Hsieh2DRotate(ustar(runs),Lo(runs),sv(runs),zo,zH, d, windir(runs), plotYN, FX, FY, F);
    
    FOOTSUM=[FOOTSUM;footsum];
    
    clear foot xx yy footRotate Fxprime Fyprime footsum
    
    if any(250:250:15000==runs)
        disp(runs)
    end
end

% figure(101)
% pcolor(XX,YY,FF);shading('interp')
% hold on
% contour(Fxprime,Fyprime,F,[0 0],'LineColor',[1 1 1])
% axis([-200 500 -500 500])

