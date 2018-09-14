function outVec=nanWinSlaughtererT(inVec,min_nans,window,num_wins,STDcutoff)

nanVec=isnan(inVec);
dsig=diff(nanVec);
outVec=inVec;
startIndex=find(dsig>0)+1;
endIndex=find(dsig<0);
if isnan(inVec(1))
    startIndex=[1;startIndex];
end

if isnan(inVec(end))
    endIndex=[endIndex;length(inVec)];
end

duration=endIndex-startIndex+1;
relevWindows=find(duration>min_nans);

for l=1:length(relevWindows)
    k=relevWindows(l);
    stdcurveleft=nan(num_wins,1);
    stdcurveright=nan(num_wins,1);
    for i=1:num_wins
        beginpt=max(1, startIndex(k)-window*i);
        endpt=min(length(inVec), endIndex(k)+window*i);
        stdcurveleft(i)=nanstd(inVec(beginpt:startIndex(k)));
        stdcurveright(i)=nanstd(inVec(endIndex(k):endpt));
    end

%     figure(66);
%     plot(stdcurveleft);
%
%     hold on
%     plot(stdcurveright,'r');
%     hold off
    CutOffLeft=find(stdcurveleft==max(stdcurveleft));
    b=max(startIndex(k)-CutOffLeft*window,1);
    r=sqrt((sum(~isnan(inVec(b:startIndex(k))))-1)/(CutOffLeft*window-1));
    if max(stdcurveleft)>=STDcutoff*r*nanstd(inVec(max(1,startIndex(k)-CutOffLeft*window-window*10):startIndex(k)-CutOffLeft*window))
        outVec(b:startIndex(k))=nan;
    end

    CutOffRight=find(stdcurveright==max(stdcurveright));
    e=min(CutOffRight*window+endIndex(k),length(inVec));
    r=sqrt((sum(~isnan(inVec(endIndex(k):e)))-1)/(CutOffRight*window-1));
    if max(stdcurveright)>=STDcutoff*r*nanstd(inVec(endIndex(k)+CutOffRight*window:min(length(inVec),endIndex(k)+CutOffRight*window+window*10)))

        outVec(endIndex(k):e)=nan;
    end
    
    

    

end








%     sizeWOnan=sum(~isnan( (inVector(startpt:endpt,1) ) ) );
%     r=sqrt((sizeWOnan-1)/(sizewindow-1));
%     trimmedSTD=trimmedSTD*r;













% figure(33);
% hold on
% plot(inVec);
% plot(outVec,'r');
end