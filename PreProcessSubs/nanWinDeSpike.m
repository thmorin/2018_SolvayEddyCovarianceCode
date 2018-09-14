function outVec=nanWinDeSpike(inVec,min_nans,trimper,window,max_STD,trash_period)
nanVec=isnan(inVec);
%nanWindow=find(nanVec(min_nans)==ones(min_nans));
dsig=diff(nanVec);
outVec=inVec;
startIndex=find(dsig>0)+1;
endIndex=find(dsig<0);

if isnan(inVec(1))
    startIndex=[1 startIndex];
end

if isnan(inVec(end))
    endIndex=[endIndex 1];
end

duration=endIndex-startIndex+1;
relevWindows=find(duration>min_nans);

for l=1:length(relevWindows)
    k=relevWindows(l);
    
    begintrash=max(1,startIndex(k)-trash_period);
    endtrash=min(length(inVec),endIndex(k)+trash_period);
    
    outVec(begintrash:startIndex(k))=nan;
    outVec(endIndex(k):endtrash)=nan;
    
    beginpt=max(1,startIndex(k)-window);
    endpt=min(length(inVec),endIndex(k)+window);
   
    
    %outVec(beginpt:startIndex(k))=windowDeSpike(inVec(beginpt:startIndex(k)),trimper,max_STD);
    %outVec(endIndex(k):endpt)=windowDeSpike(inVec(endIndex(k):endpt),trimper,max_STD);
end
end



function outVec=windowDeSpike(inVec,trimper,max_STD)
    [trimmedmean trimmedstd]=trimmean2(inVec,trimper);
    outVec=inVec;
    outVec(abs(inVec-trimmedmean)>max_STD*trimmedstd)=nan;
end