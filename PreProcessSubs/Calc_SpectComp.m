function [A, B, C] = Calc_SpectComp( T, P, x_v )


% This function is calculating the coefficiants that compensate the CH4 flux 
% that are caused due to spectroscopic light absorption measurements for 
% correlated fluctuations in temperature and water vapor density.
% Inputs:
% x_v = From LI7500. mol fraction of water vapor to dry air [ Unitless ]
% T = From CSAT Sonic. Temperature (Tr from kymal and gaynor) [ Deg C ]
% P = LI7700 pressure  [ Pa ]
% 
% Outputs:
% A = Dimensionless. Corrects methane density for spectroscopic A
% effects of T, P and water vapor.
% 
% B = 1 + (1 − 1.46 xv )αv Pe (κ Pe)/κ
% Dimensionless. Corrects the latent heat B flux term for spectroscopic effects.         
%   
%                    
% C = 1 + (1 − xv) T κT/κ + xv ( B − 1) . 
% Dimensionless. Corrects the sensible heat term for spectroscopic effects.
%   
% 
% Mole fraction of water vapor, dimensionless. xv = ρv ρ = ρv RT P . From the LI-7500, T ,xv and P .
% Water vapor broadening factor, dimensionless. α v = a v − 1 = 0.46 . a v
% Specific to the αv i s the foreign gas broadening coefficient for water vapor on methane,
%  LI-7700. See page 6-6. a v = 1.46 .
%          ∂κ
% κT κT =    . Rate of change in κ w ith T at constant Pe = Pe , K-1. Given in Table 5-4.
%         ∂T
%   
%             ∂κ
%      κ Pe =       . Rate of change in κ w ith Pe at constant T = T , kPa-1.
% κ Pe                                                                        Given in Table 5-5.
%             ∂Pe T
%
%               -2   -1 correct methane flux computed
% Fc = mg CH 4 m    h
%
%
%                                                   0         1             2
% Q = a0 × T 2 + a1 × T + a2               %               
                                           % a     0        -1.3×10-7      3.7×10-5
% R = b0 × T 2 + b1 × T + b2               %       
                                           % b    4.0×10-8   1.1×10-5     2.18×10-3
% S = c0 × T 2 + c1 × T + c2               %                             
                                           % c   2.0×10-6     9.8×10-4    0.378
                                           
                                           
    Pe = P*(1 + 0.46*x_v);    % Pe = equivalent pressure  ( Pe = P 1 + 0.46 x_v )                                       


    a0 = 0 ;    a1 = -1.3E-7 ; a2 = 3.7E-5 ;
    b0 = 4E-8 ; b1 = 1.1E-5  ; b2 = 2.18E-3 ;
    c0 = 2E-6 ; c1 = 9.8E-4  ; c2 = 0.378 ;

    Q = a0 * T.^2 + a1 * T + a2;     % T in [deg C]
    R = b0 * T.^2 + b1 * T + b2;
    S = c0 * T.^2 + c1 * T + c2;

    A = Q .* Pe.^2 + R .* Pe + S;

    a=-8.2E-6;
    b=4.3E-3;
    c=-1.7E-4;
    d=0.03;

    %        κ Pe
    %             = ( a × T + b ) × Pe + ( c × T + d )
    % α v Pe
    %         κ

    B = 1+(1-1.46*x_v).*( a*T +b ).*Pe + c*T +d ;
    % 
    %       (                   )      (                    )
    %   κT
    %      = a × T 2 + b × T + c × Pe + d × T 2 + e × T + f   5-5
    % T
    %    κ

    a=-4.0E-8;
    b=1.55E-5;
    c=-7.0E-3;
    d=-4.7E-6;
    e=3.0E-3;
    f=0.927;

    C = 1+(1-1.45*x_v).*(( a* T.^2 + b*T + c ).* Pe + d* T.^2 + e*T + f )+x_v.*(B-1);

end


