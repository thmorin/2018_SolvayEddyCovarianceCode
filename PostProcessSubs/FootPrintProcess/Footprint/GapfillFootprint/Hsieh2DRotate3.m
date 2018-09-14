%Hsieh2D.m :     Computes and plot the 2-D footprint based on the
%                semi-analitycal model Hsieh and Katul (2000) and 
%                extended 2-d by Detto et al. (2006)
%
%        ustar   = friction velocity [m/s]
%
%        Lo      = Obukhov stability length [m]
%            
%        sv      = fluctuation of lateral wind [m/s] (approx. 2 times ustar)
%                       
%        zo      = momentum roughness length [m]
%
%        zH      = height of the measurement [m]
%          
%        d - displacement height
%
%       windir - wind direction, from north, degrees
%       
%        FX - meshgrid, x coordinates, of mask map
%
%        FY - meshgrid, y coordinates, of mask map
%
%        F - mask map on FX,FY coordinates
%
%        usage:[foot]=Hsieh2D(ustar,Lo,sv,zo,zm)

% Dummy inputs:  clear all; close all;ustar = 0.45;Lo = -18.14;sv = 0.91;zo = 0.2;d= 2/3*1.7;zH =3.7; windir = 250; plotYN=0;
%load('FASET.mat')
%FX=(-1000:1000);FY=(-1000:1000);F=zeros(2001,2001);F(1:1002,500:1500)=1;
%load('../MERImap.mat');FX=MERImap.FX;FY=MERImap.FY;F=MERImap.mask;

function [foot,xx,yy,footRotate,Fxprime,Fyprime,footsum]=Hsieh2DRotate(ustar,Lo,sv,zo,zH, d, windir, plotYN, FX, FY, F);

zm=zH-d;
k=0.4;
Lx=100*zm;
%wd = windir/180*pi;

zu=zm*(log(zm/zo)-1+zo/zm);       
P=[0.59 1 1.33];                  
D=[0.28 0.97 2.44];
Mu=[100 500 2000];
stab=zu/Lo;
thresh=0.04;

PI = PatchIndex(F);

if stab<-thresh
    ii=1;
elseif abs(stab)<thresh
    ii=2;
elseif stab>thresh
    ii=3;
end
D1=D(ii);
P1=P(ii);
Mu1=Mu(ii);

% min_x=(Mu1/100)*zo;
% max_x=Mu1*zm;
bin=max(0.5,floor(Lx/500)); % SRG comment
x=[eps:bin:Lx];   %SRG comment

c1=(-1/k/k)*(D1*zu^P1*abs(Lo)^(1-P1))./(x);
Fc=exp(c1);
Fp=-(c1./x).*Fc;
Xp=(1/2/k/k)*(D1*zu^P1*abs(Lo)^(1-P1));
F2H=(D1/0.105/k/k)*(zm^(-1)*abs(Lo)^(1-P1)*zu^(P1));
Xm=min([F2H*zm Lx]);

nn=floor((Xm+1)/bin)-1;
sy=zo*0.3*(sv/ustar).*(x./zo).^0.85;
b=floor(zo*(0.3*(sv/ustar).*(Xm./zo).^0.85)/1.5);
y=(-2*b:bin:2*b);

foot=nan(nn,length(y));
for i=1:nn
    foot(i,:)=Fp(i)*normpdf(y,0,sy(i));
end

% Unrotated footprint
foot=foot';

[xx,yy] = meshgrid(x(1:nn),y);
if plotYN==1
    figure(1)
    pcolor(xx,yy,foot);shading('interp')
    xlabel('distance from the tower [m]')
    ylabel('lateral spread [m]')
end

% % rotate the wind direction to where it is going to rather than where it is
% % coming from
% if winddir <= 180
%     winddir2 = winddir + 180;
% elseif winddir > 180
%     winddir2 = winddir - 180;
% end

% coordinate rotation
% [Fxprime,Fyprime] = rotateToWind(FX,FY,windir,F,1);
[FX2, FY2] = meshgrid(FX,FY);
[Fxprime,Fyprime] = rotateToWind(FX2,FY2,windir);

Rl= length(FY);
Cl=length(FX);

if plotYN==1
    figure(2)
    pcolor(FX2,FY2(Rl:-1:1,:),F);shading('interp')
    xlabel('distance upwind [m]')
    ylabel('distance crosswind [m]')
    title('Area of interest')
end

footRotate = interp2(xx,yy,foot,Fxprime,Fyprime);
% plot the rotated footprint

if plotYN==1
    figure(3)
    pcolor(FX2,FY2(Rl:-1:1,:),footRotate);shading('interp') % this doesn't look right to me...stretched, but maybe it is how it is being plotted
end

%binsize=((max(max(Fxprime))-min(min(Fxprime)))/(Cl-1))*((max(max(Fyprime))-min(min(Fyprime)))/(Rl-1));
binnorm = (nansum(nansum(foot))*(x(2)-x(1))*(y(2)-y(1)))/nansum(nansum(footRotate));

footsum=zeros(size(PI));
for i=1:length(PI)
    footsum(i) = nansum(nansum(footRotate.*(F==PI(i))))*binnorm;
end


end




