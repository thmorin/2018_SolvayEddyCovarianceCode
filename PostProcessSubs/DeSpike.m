function [fspkVector, FlagVec]  = DeSpike (inVector,sizewindow,STDMax,MinVec,MaxVec,TrimPercent,inst,diagno,RSSI, RSSIper) % Modifying to use trimmean

% Identifies and eliminates spikes found within inVector
%
% inVector    = input spiked column vector
%
% sizewindow    = Size of the window around the obsevation
%
% STDMax   = Times variance to identificate a spike
%
% MinVec, MaxVec  = min and max allowed values for inVector
%
% TrimPercent    = TrimPercent/2 highest and lowest points not to be used
% for mean or STD
%
% inst = instrument (string). (optional input)
% diagno = diagnostics of instrument. (optional input)
%
% Output :
% fspkVector  = copy of inVector with spikes replaced with nan
min_length=0;
percentoverlap=0.5; %never change this
% stdOfVec=nan(floor(length(inVector)/sizewindow),1);
% for k=1:floor(length(inVector)/sizewindow)
%     stdOfVec(k)=nanstd((inVector(sizewindow*(k-1)+1:sizewindow*k,1)));
% end
% [~,X]=hist(stdOfVec);
% minSTD=X(2);
%killpercent=0.85;
%TO FIX ON THIS FUNCTION:
%1. Overlap can never be overwritten. This is a problem if you use an overlap
%greater than 50%. I believe there is also a problem if the overlap is less
%than 50%. ex: first pass flags a value. third (so same flag series) does
%not want to flag that value. there is no functionality to deflag it.

spkVector=false(size(inVector));

opts.instrument = 'NONE';
opts.diagnostic = ones(size(inVector));

%  NARGIN Number of function input arguments.
%    Inside the body of a user-defined function, NARGIN returns
%    the number of input arguments that were used to call the function.

if nargin==6       % No Instrument diagnostics
    instrumentDiagnostic = ones( size( inVector ));
elseif nargin==8
    opts = parseargs( opts , 'instrument', inst , 'diagnostic' , diagno );
    instrumentDiagnostic = diagnosticParameters( opts.instrument , opts.diagnostic );
    min_length=5;
elseif nargin==10
    opts = parseargs( opts , 'instrument', inst , 'diagnostic' , diagno );
    instrumentDiagnostic = diagnosticParameters( opts.instrument , opts.diagnostic, RSSI, RSSIper );
    min_length=5;
end

tmpspkVector1=false(size(inVector));
tmpspkVector2=false(size(inVector));

spkVector( (inVector > MaxVec) | (inVector < MinVec) | isnan(inVector) | ~instrumentDiagnostic)=1;

% FlagVec=zeros(length(inVector)+sizewindow,1);
FlagVec=zeros(size(inVector));
inVec=inVector;
% inVecmean=nanmean(inVec);
% inVec=[inVecmean*ones(sizewindow,1) ; inVec];
inVec(spkVector) = nan;

startpt=1;
CW=0;
endpt=sizewindow;
tmpspkVector2(1:floor(endpt*(1-percentoverlap)))=1;
nwindows = round(length(inVec)/sizewindow)+1;
trimmedMean = nan(nwindows*2,1);
trimmedSTD = nan(nwindows*2,1);
trimmedSTDMin = nan(nwindows*2,1);

while endpt<length(inVec)
    CW=CW+1;
    %If last window, extend window to include remainder of vector
    if (endpt+floor(sizewindow*(1-percentoverlap)))>length(inVec)
        endpt=length(inVec);
    end
    inVec_=inVec(startpt:endpt);
    [trimmedMean(CW),~] = trimmean2(inVec_,TrimPercent);
    inVec_=nandetrend(inVec_);
    %trimmedSTD=max(minSTD,trimmedSTD);
    [~, trimmedSTD(CW)] = trimmean2(inVec_,TrimPercent);
     startpt = startpt+floor(sizewindow*(1-percentoverlap));
    endpt = startpt + sizewindow-1;
    
end


if length(trimmedMean) < 11
    trimmedSTDMin(:) = prctile(trimmedSTD(~isnan(trimmedSTD)),25);
else

endwin = 10;
startwin = 1;
while endwin < length(trimmedMean)
    win10 = trimmedSTD(startwin:endwin);
    trimmedSTDMin(startwin:endwin) = prctile(win10(~isnan(win10)),25);
    startwin = startwin+20;
    endwin = min(startwin+19,length(trimmedMean));
end
    win10 = trimmedSTD(startwin:endwin);
    trimmedSTDMin(startwin:endwin) = prctile(win10(~isnan(win10)),25);
end
startpt=1;
CW=0;
endpt=sizewindow;
while endpt<length(inVec)
    CW=CW+1;
     if (endpt+floor(sizewindow*(1-percentoverlap)))>length(inVec)
        endpt=length(inVec);
     end
     inVec_=inVec(startpt:endpt);
     trimmedSTDA= max(trimmedSTD(CW),trimmedSTDMin(CW));
     
    if mod(CW,2)~=0
       
        tmpspkVector1(startpt:endpt) = abs(inVec_-trimmedMean(CW))>(trimmedSTDA*STDMax);
        
        if isnan(trimmedSTDA)
            tmpspkVector1(startpt:endpt)=1;
        end
    else
        tmpspkVector2(startpt:endpt) = abs(inVec_-trimmedMean(CW))>(trimmedSTDA*STDMax);
        
        if isnan(trimmedSTDA)
            tmpspkVector2(startpt:endpt)=1;
        end
        %Returns flagged if in vector is all NaN's (or eliminated)
    end
    startpt = startpt+floor(sizewindow*(1-percentoverlap));
    endpt = startpt + sizewindow-1;
end
%Sets startpt back one iteration to set complementary matrix end to 1's
startpt=startpt-floor(sizewindow*(1-percentoverlap));
%endpt=startpt+sizewindow-1;

if mod(CW,2)==0
    tmpspkVector1((startpt+sizewindow*percentoverlap):end)=1;
else
    tmpspkVector2((startpt+sizewindow*percentoverlap):end)=1;
end


spkVector = or(and(tmpspkVector1,tmpspkVector2),spkVector);

dsig=diff(spkVector);
startIndex=find(dsig<0)+1;
endIndex=find(dsig>0);

if ~spkVector(1)
    startIndex=[1; startIndex];
end

if ~spkVector(end)
    endIndex=[endIndex; length(inVec)];
end

duration=endIndex-startIndex+1;
relevWindows=find(duration<=min_length);

for l=1:length(relevWindows)
    k=relevWindows(l);
    spkVector(startIndex(k):endIndex(k))=1;
end

% QA doesn't account for the nanWinSlaughterer AND needs to be added
FlagVec(and(tmpspkVector1,tmpspkVector2))=5; % 5= nan by DeSpike
FlagVec(isnan(inVector))=6;                     % 6= nan in Original Vector
FlagVec((inVector > MaxVec))=4;                 % 4= nan by MaxVec
FlagVec((inVector < MinVec))=3;                 % 3= nan by MinVec
FlagVec(instrumentDiagnostic==0)=2;          % 2= nan by Diagnostics
fspkVector=inVector;
fspkVector(spkVector==1)=nan;
end