function [Ts]= Schotanus(T,U,V,W)

% Schotanus 1983 and Liu et al 2000, revisited by Detto
% This corrects the sonnic temp for geometric errors due to speed of sound

Ts=T+(3/4*U.^2+3/4*V.^2+1/2*W.^2)/402.798;  

end
