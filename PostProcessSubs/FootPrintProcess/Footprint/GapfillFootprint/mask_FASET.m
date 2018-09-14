 FASET=csvread('../FASET/FASET_Grid_outline.csv');
 
 
 Fy = (min(FASET(:,2)):25:max(FASET(:,2)));
 Fy2000 = -1000:1:1000;
 Fx2000 = Fy2000;
 [FX, FY] = meshgrid(Fx2000,Fy2000);
 F = zeros(size(FX));
 
 
 for yi = Fy(1:end-1)
     currenty=find(FASET(:,2)==yi);
     fromx = min(FASET(currenty,1));
     tox=max(FASET(currenty,1));
     F(1001+yi:1026+yi,1001+fromx:1001+tox)=1;
 end
 
 save('FASET.mat','FX','FY','F')
 pcolor(FX,FY,F);shading('interp')