function [ dCdt,Stor ] = CarbonStorageFlux( Cs, s_d, dt)
% [ dCdt ] = CarbonStorageFlux( Cs, theta)
% Written by Ashley Matheny 06/2014
% Modified for carbon storage by Tim Morin, 07/2014
%
% Uses vertical profile of carbon concentration measurements and their
% heights to calculate the change in carbon storage over time using the 
% trapezoidal rule.
%
% Cs   = m X n matrix of carbon concentration values at different heights
%        m=time, n=heights), 1:n in ascending order
% s_d  = 1 X n vector of heights detailing the location of each carbon
%        concentration sensor on the tower, 1:n in ascending order
% dCdt = m X 1 vector of overall carbon storage change per time step input

%%
%%----------- Constants ----------------------------------------------------
if nargin<3
	dt = 30*60;                                                        % [s] seconds per half hour time step
end
dx=s_d-[0 s_d(1:end-1)];                                                   % [m] distance b/w each sensor

%%
%----------- Preallocate---------------------------------------------------
[time, depth] = size(Cs);
dCdt=nan(time,1);

%%
%------------- Calculate storage in entire soil column --------------------
Stor = [Cs(:,1) (Cs(:,2:depth)+Cs(:,1:depth-1))/2]*dx';                    % [ppm] Are these units correct?

%%
%----------- Change in C with respect to time at each layer ---------------
dCdt(2:time,:) = (Stor(2:time,:)-Stor(1:time-1,:))/dt;                     % [ppm/s] T2-T1
dCdt(1,:)=0;                                                               % [ppm/s] 0 for first half hour

end
