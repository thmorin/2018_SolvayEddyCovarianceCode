function [fluxW, u, v, w, alfa, theta] = WindData( U, V, W ,Sonic, MagneticCorr, Azimuth)

% function [fluxW, u, v, w, alfa, theta] = WindData( SonicTmp, U, V, W, NFavg )
% Calculates the Wind related data. 
% Outputs :
% Wind Flux data, rotated wind u, v, w with angles of rotation alfa and theta.
% Inputs are : 
% Sonic Temperature
% Sonic U,V,W
% NFavg - Number of measurments that are being avareged


% Good_Data = find( ~isnan( SonicTmp ) & ~isnan( U ) & ~isnan( V ) & ~isnan( W ) ); % Sonic
% 
%     if length(Good_Data) > 0.5 * NFavg   % NFavd = Number of fast values that are meaned into one value per window
if nargin<4
    Sonic='';
    MagneticCorr=0;
    Azimuth=0;
    disp('Warning: Magnetic correction and azimuth not specified. Defaulting to 0');
end

MagneticCorr=MagneticCorr*pi/180;
Azimuth=Azimuth*pi/180;

        % -- 3D-rotation the coordinate system in the window----------------
        [ u, v, w, alfa, theta ] = RotateWind3D( U, V, W, Sonic);   % Sonic 

        theta=mod(theta+MagneticCorr+Azimuth,2*pi);
        
        % -- Means of Wind
        fluxW.ubar = nanmean( u );      % Sonic wbar = 0; vbar=0
        fluxW.wbar = nanmean( W );      % Sonic Mean unrotated W
        un = u - fluxW.ubar   ;         % Sonic Pertubations of u from mean u
        fluxW.uu = nanvar( u );         % Sonic variance of u
        fluxW.vv = nanvar( v );         % Sonic variance of v 
        fluxW.ww = nanvar( w );         % Sonic variance of w 
        fluxW.u3 = moment( u, 3 );      % Sonic moment of u
        fluxW.v3 = moment( v, 3 );      % Sonic moment of v
        fluxW.w3 = moment( w, 3 );      % Sonic moment of w
        fluxW.uw = nanmean( w.*un );    % Sonic u'w'
        fluxW.vw = nanmean( w.*v );     % Sonic v'w'
        fluxW.U_bar = nanmean( U );
        fluxW.V_bar = nanmean( V );
        fluxW.UW = nanmean((U-fluxW.U_bar).*(W-fluxW.wbar));
        fluxW.VW = nanmean((V-fluxW.V_bar).*(W-fluxW.wbar));
        % -- Horizontal rotation angle. add sonic azimuth in radians, make sure you count the wind dirrection the same.
        fluxW.WD  = theta;  % Sonic 
        fluxW.WDz = alfa;   % Sonic Vertical rotation angle

        % -- U*
        fluxW.ustar=( fluxW.uw^2 + fluxW.vw^2 )^0.25;  % Sonic



%     elseif length( Good_Data ) <= 0.5*NFavg
%         theta = nan; alfa = nan; 
%         u = nan( NFavg,1 );  % Sonic - Bad Data
%         v = nan( NFavg,1 );  % Sonic - Bad Data
%         w = nan( NFavg,1 );  % Sonic - Bad Data
%         fluxW.ubar = nan;    % Sonic wbar = 0; vbar=0
%         fluxW.wbar = nan;    % Sonic Mean unrotated W
%         fluxW.uu = nan;      % Sonic variance of u
%         fluxW.vv = nan;      % Sonic variance of v 
%         fluxW.ww = nan;      % Sonic variance of w 
%         fluxW.u3 = nan;      % Sonic moment of u
%         fluxW.v3 = nan;      % Sonic moment of v
%         fluxW.w3 = nan;      % Sonic moment of w
%         fluxW.uw = nan;      % Sonic u'w'
%         fluxW.vw = nan;      % Sonic v'w'
% 
%         % -- Horizontal rotation angle. add sonic azimuth in radians, make sure you count the wind dirrection the same.
%         fluxW.WD  = nan;     % Sonic 
%         fluxW.WDz = nan;     % Sonic Vertical rotation angle
% 
%         % -- U*
%         fluxW.ustar = nan;   % Sonic
%     end
end

