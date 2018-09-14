function [wc2,wc3] = gapfill(wc,m,season) %Changed to have maximum and minimum thresholds

NN=length(wc);
n=season*m;                           %# of points
%R=floor(NN/n);                         % # of realizzations
WC=wc;
wc2=wc;

% Gap-filling with 2-D (griddata) look-up table (along time and day)
run=find(~isnan(wc2));
%gap=find(isnan(wc2));

nday=NN/season;
XI(1,:)=1:nday;
YI(:,1)=1:season;
X=repmat(XI,season,1);
X=reshape(X,1,nday*48);
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

% figure(1);clf
% plot(wc3,'k');hold all
% plot(wc2)
% plot(gap,WC(gap),'k+')
% pause(.1)


%nested function
%standarizzation of time series
function [ys,stnd] = standard(y,period)

n=length(y);
y1(:,1)=y;

Y=reshape(y1,period,n/period);
mY=nanmedian(Y');
sY=nanstd(Y');

stnd=repmat(mY',n/period,1);
ys=(y-stnd);%./repmat(sY',n/period,1);

% use=find(isnan(ys)==1);
% ys(use)=0;

end

%find gaps<=M and fill with spline iterpolant

function y_gf=gap_filler(x,y,M)

if isnan(y(1))==1;
    y(1)=y(2);
end

y_gf=y;
use=find(isnan(y)==0);
gap_M=[];

for i=1:length(y)-M-1
    if isnan(y(i))==1 & sum(isnan(y(i:i+M+1)))<M & isnan(y(i-1))==0
       b=find(isnan(y(i:end))==0);
       gap_M=[gap_M i:i+b(1)-2];
       i=i+b(1)-1;
    end
end
  

y_gf(gap_M)=interp1(x(use),y(use),x(gap_M),'spline');

end

%normilize the variable x with zero mean and unit vatiance

function xn=red(x)


     if nanstd(x)==0
       xn=x*0;
     else
       xn=(x-nanmean(x))/nanstd(x);
     end
end

end