function [u ,v, w ,alfa, theta,Az,Ze,R ] = RotateWind3D(U, V, W, Sonic, cont)

          if strcmp(Sonic,'CSAT')
              V=-V;
          elseif strcmp(Sonic,'CSAT3B')
              V=-V;    
          elseif strcmp(Sonic,'RMY')
              X=U;
              U=V;
              V=X;
          else
              disp('Warning: No sonic specified. Running with default coordinate system');
          end

%-- 3D-rotation the coordinate system in the window----------------
          theta=atan2(nanmean(V),nanmean(U));
          alfa=atan2(-nanmean(W),(nanmean(U)^2+nanmean(V)^2)^0.5);
          Az=[cos(theta) sin(theta) 0; -sin(theta) cos(theta) 0; 0 0 1];
          Ze=[cos(alfa) 0 -sin(alfa); 0 1 0; sin(alfa) 0 cos(alfa)];
          R=Ze*Az;
          u=[U V W]*R(1,:)';
          v=[U V W]*R(2,:)';
          w=[U V W]*R(3,:)';
          if nargin == 5
              if cont == 1
                  disp(['sumU = ' num2str(nansum(u))]);
                  disp(['sumV = ' num2str(nansum(v))]);
                  disp(['sumW = ' num2str(nansum(w))]);
              end
          end
end
