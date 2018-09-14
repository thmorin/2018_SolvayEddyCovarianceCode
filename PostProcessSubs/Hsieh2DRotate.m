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
%       windir - wind direction, from north, degrees. Where wind is coming
%       FROM, not going to
%       
%        FX - meshgrid, x coordinates, of mask map
%
%        FY - meshgrid, y coordinates, of mask map
%
%        F - mask map on FX,FY coordinates
%
%        usage:[foot]=Hsieh2D(ustar,Lo,sv,zo,zm)

% Dummy inputs:  clear all; close all;ustar = 1;Lo = -5;sv = 2;zo = 0.1*8;d= 0.76*8;zH =19.6; windir =90; plotYN=1;
%load('FASET.mat')
%FX=(-1000:1000);FY=(-1000:1000);F=zeros(2001,2001);F(1:1002,500:1500)=1;
%load('../MERImap.mat');FX=MERImap.FX;FY=MERImap.FY;F=MERImap.mask;
 %load('ORWRP_Map.mat'); F = WetlandMap; FX=Fx1(1,:)';FY = Fy1(:,1);
 
function [foot,xx,yy,footRotate,Fxprime,Fyprime,footsum]=Hsieh2DRotate(ustar,Lo,sv,zo,zH, d, windir, plotYN, FX, FY, F)

zm=zH-d;
k=0.4;
Lx=100*zm;
%wd = windir/180*pi;

zu=zm*(log(zm/zo)-1+zo/zm);       
P=[0.59 1 1.33];                  
D=[0.28 0.97 2.44];
% Mu=[100 500 2000];
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
% Mu1=Mu(ii);

% min_x=(Mu1/100)*zo;
% max_x=Mu1*zm;
bin=max(0.5,floor(Lx/500)); % SRG comment
x=eps:bin:Lx;   %SRG comment

c1=(-1/k/k)*(D1*zu^P1*abs(Lo)^(1-P1))./(x);
Fc=exp(c1);
Fp=-(c1./x).*Fc;
% Xp=(1/2/k/k)*(D1*zu^P1*abs(Lo)^(1-P1));
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

if ~all(y==0)
[xx,yy] = meshgrid(x(1:nn),y);
% if plotYN==1
% %     figure(1)
%     pcolor(xx,yy,foot);shading('flat')
%     xlabel('distance from the tower [m]')
%     ylabel('lateral spread [m]')
% end

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
% Cl=length(FX);

if plotYN==1
figure;%     figure(2)
    pcolor(FX2,FY2(Rl:-1:1,:),F);shading('interp'); hold on
    xlabel('distance upwind [m]')
    ylabel('distance crosswind [m]')
    title('Area of interest')
end

footRotate = interp2(xx,yy,foot,Fxprime,Fyprime);
% plot the rotated footprint

if plotYN==1
figure(gcf);%     figure(2)
contour(FX2,FY2(Rl:-1:1,:),footRotate,'LevelStep',0.0002,'Fill','off');shading('interp') 

%Testing: contour(FX2,FY2(Rl:-1:1,:),footRotate,'LevelStep',0.0002,'Fill','off');shading('interp') 
% text(0,0,['Wind direction=' num2str(windir)]);
% text(0,0,['Lo=' num2str(Lo)]);
% text(0,0,['sv=' num2str(sv)]);
% text(0,0,['ustar=' num2str(ustar)]);
%     figure(3)
%     pcolor(FX2,FY2(Rl:-1:1,:),footRotate);shading('interp');colorbar 
hold off;
end

%binsize=((max(max(Fxprime))-min(min(Fxprime)))/(Cl-1))*((max(max(Fyprime))-min(min(Fyprime)))/(Rl-1));
binnorm = (nansum(nansum(foot))*(x(2)-x(1))*(y(2)-y(1)))/nansum(nansum(footRotate));

footsum=zeros(size(PI));
for i=1:length(PI)
    footsum(i) = nansum(nansum(footRotate.*(F==PI(i))))*binnorm;
    
end
else
    footsum=zeros(size(PI));
    foot=nan;
    xx=nan;
    yy=nan;
    footRotate=nan;
    Fxprime=nan;
    Fyprime=nan;
    disp('Warning: footprint could not run. Will not return data for this half hour');
end

end




