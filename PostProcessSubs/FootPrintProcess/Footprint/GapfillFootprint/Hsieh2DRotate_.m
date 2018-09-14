function [foot,xx,yy,footRotate,Fxprime,Fyprime,footsum]=Hsieh2DRotate_(ustar,Lo,sv,zo,zH, d, windir, plotYN, FX, FY, F, DOY, HalfHour, Year)

%Hsieh2D.m :     Computes and plot the 2-D footprint based on the
%                semi-analitycal model Hsieh and Katul (2000) and
%                extended 2-d by Detto et al. (2006)
%
%       ustar   = friction velocity [m/s]
%
%       Lo      = Obukhov stability length [m]
%
%       sv      = fluctuation of lateral wind [m/s] (approx. 2 times ustar)
%
%       zo      = momentum roughness length [m]
%
%       zH      = height of the measurement [m]
%
%       d - displacement height
%
%       windir - wind direction, from north, degrees
%
%       plotYN - 0=don't plot, 1=plot.
%
%       FX - meshgrid, x coordinates, of mask map
%
%       FY - meshgrid, y coordinates, of mask map
%
%       F - mask map on FX,FY coordinates
%
%        usage:[foot]=Hsieh2D(ustar,Lo,sv,zo,zm)
%
% Dummy inputs:  clear all; close all;ustar = 0.45;Lo = -18.14;sv = 0.91;zo = 0.2;d= 2/3*1.7;zH =3.7; windir = 250; plotYN=0;
%load('FASET.mat')
%FX=(-1000:1000);FY=(-1000:1000);F=zeros(2001,2001);F(1:1002,500:1500)=1;
%load('../MERImap.mat');FX=MERImap.FX;FY=MERImap.FY;F=MERImap.mask;

zm=zH-d; % Hight of measurment minus displacment hight
k=0.4;
Lx=100*zm;

zu=zm*(log(zm/zo)-1+zo/zm); %Don't understand the science here
P=[0.59 1 1.33];
D=[0.28 0.97 2.44];
%Mu=[100 500 2000];
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
%Mu1=Mu(ii);

% min_x=(Mu1/100)*zo;
% max_x=Mu1*zm;
bin=max(0.5,floor(Lx/500)); % SRG comment
x=(eps:bin:Lx);   % SRG comment

c1=((-1/k/k).*(D1*zu^P1*abs(Lo)^(1-P1)))./(x); %Inside the brackets part of Equation17 Hsieh 2000
Fc=exp(c1);
Fp=-(c1./x).*Fc; %Equation 17 in Hsieh 2000
%Xp=(1/2/k/k)*(D1*zu^P1*abs(Lo)^(1-P1));
F2H=(D1/0.105/k/k)*(zm^(-1)*abs(Lo)^(1-P1)*zu^(P1)); %Fetch to height ratio
Xm=min([F2H*zm Lx]); %Caps it at either 100 times the tower above the displacement height or as calculated above

nn=floor((Xm+1)/bin)-1;
sy=zo*0.3*(sv/ustar).*(x./zo).^0.85; %Cross wind component - Detto 2006, EQN B4
b=floor(zo*(0.3*(sv/ustar).*(Xm./zo).^0.85)/1.5); %Used for lateral spacing
if b~=0
    
    y=(-2*b:bin:2*b);
    
    foot=nan(nn,length(y));
    for i=1:nn
        foot(i,:)=Fp(i)*normpdf(y,0,sy(i));
    end
    
    % Unrotated footprint
    foot=foot';
    
    [xx,yy] = meshgrid(x(1:nn),y);
    
    % % rotate the wind direction to where it is going to rather than where it is
    % % coming from
    % if winddir <= 180
    %     winddir2 = winddir + 180;
    % elseif winddir > 180
    %     winddir2 = winddir - 180;
    % end

    
    [Fxprime,Fyprime] = rotateToWind(FX,FY,windir);
    
    Rl= length(FY(1,:));
    %Cl=length(FX(:,1));
    
    footRotate = interp2(xx,yy,foot,Fxprime,Fyprime,'nearest');
    
    % Plot the map and contour
    if plotYN==1
        figure1=figure('XVisual',...
            '0x23 (TrueColor, depth 24, RGB mask 0xff0000 0xff00 0x00ff)',...
            'NumberTitle','off');
        axes1 = axes('Parent',figure1,'YDir','reverse','Layer','top','FontSize',14);
%         xlim(axes1,[-100 100]);
%         ylim(axes1,[-100 100]);
        box(axes1,'on');
        hold(axes1,'all');
        image(FX(1,:),FY(:,1),F,'Parent',axes1,'CDataMapping','scaled')
        xlabel('Distance Upwind [m]','FontSize',16)
        ylabel('Distance Crosswind [m]')
        [ month , day ] = DayOY2Date( DOY , Year );
        v = [Year, month, day, floor(HalfHour/2), mod(HalfHour,2)*30,0];
        title(['     Area of interest - ORWRP    ',sprintf('\n')]...
            ,'FontSize',20,'Editing','on',...
            'BackgroundColor',[0.831372559070587 0.815686285495758 0.7843137383461])
        % Create textbox
        annotation(figure1,'textbox',...
        [0.647531100478468 0.93573844419391 0.174239234449761 0.0496054114994355],...
        'String',{datestr(v, 'mmm-dd-yyyy HH:MM' )},...
        'FontSize',18,...
        'FitBoxToText','off',...
        'LineStyle','none',...
        'BackgroundColor',[0.831372559070587 0.815686285495758 0.7843137383461]);
        hold on
        % Plotting the footprint cotour on top of the map.
        contour(FX,FY(Rl:-1:1,:),footRotate,'LineWidth',1)% ,'LineColor',[1 1 1]); 
        disp(windir)
    end
    
    %binsize=((max(max(Fxprime))-min(min(Fxprime)))/(Cl-1))*((max(max(Fyprime))-min(min(Fyprime)))/(Rl-1));
    binnorm = (nansum(nansum(foot))*(x(2)-x(1))*(y(2)-y(1)))/nansum(nansum(footRotate));
    
    footsum=zeros(1,length(PI));
    for i=1:length(PI)
        footsum(i) = nansum(nansum(footRotate.*(F==PI(i))))*binnorm;
    end
else
    foot = nan;
    xx = nan;
    yy = nan;
    footRotate = nan;
    Fyprime = nan;
    Fxprime = nan;
    footsum = nan(1,length(PI));
end

end
