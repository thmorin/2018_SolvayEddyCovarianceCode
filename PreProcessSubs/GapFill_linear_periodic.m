
function [wc2,wc3] = GapFill_linear_periodic(wc,season)

%wc - data to be gapfilled
%season - length of repeat period, in number of observations of wc

NN=length(wc);
%n=season*m;                           %# of points

WC=wc;
wc2=wc;


% Gap-filling with 2-D (griddata) look-up table (along time and day)
nday=floor(NN/season);
if nday*season<NN
    NN=nday*season;
    wc2=wc2(1:NN);
end
run=find(~isnan(wc2));

XI(1,:)=1:nday;
YI(:,1)=1:season;
X=repmat(XI,season,1);
X=reshape(X,1,nday*season);
Y=repmat(YI,nday,1);
Z=wc2;

method='linear';
ZI = griddata(X(run),Y(run),Z(run),XI,YI,method);
wc3=reshape(ZI,NN,1);
wc3(run)=wc2(run);


%estrapol using nearest
run=find(~isnan(wc3));
gap=find(isnan(wc3));
if ~isempty(gap)
    Z=wc3;
    method='nearest';
    ZI = griddata(X(run),Y(run),Z(run),XI,YI,method);
    Z=reshape(ZI,NN,1);
    wc3(gap)=Z(gap);
end

figure(1);clf
plot(wc3,'k');hold all
plot(wc2)
plot(gap,WC(gap),'k+')
pause(.1)



end