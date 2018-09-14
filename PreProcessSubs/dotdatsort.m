function [fileCell1, filecode1]=dotdatsort(filenamelist_,Day1,Dayend,tsOrmet)

% fileCell1 is the cell output. col 1: the relevant DOYs in folder. col 2:
% the firts filename of that DOY. col 3: the second filename of that DOY.
% Col 3 : the 3rd filename of that DOY etc.
% filecode=0 no.dat for this DOY. filecode=1, one .dat file for that day.
% filecode=2 2 .dat files for that DOY etc.
%
% for ts_data year is (9:12)
%             month is (14:15)
%              day is (17:18)
%
% for met_data year is (10:13)
%             month is (15:16)
%              day is (18:19)
% tsOrmet='ts';


%  filenamelist_=fast;



[numOfile, ~]=size(filenamelist_);
filenamelist=cell(numOfile,1);
for k=1:numOfile
    filenamelist{k,1}=filenamelist_(k).name;
end
    


if nargin==4
    namestr=length(tsOrmet);
end
[~,listFL]=size(filenamelist);
listmat=nan(listFL,4);

for k=1:length(filenamelist)
    if ~strcmp('orig',filenamelist{k})
    if strcmp('ts',filenamelist{k}(1:2))
        listmat(k,2)=str2double(filenamelist{k}(9:12));
        listmat(k,3)=str2double(filenamelist{k}(14:15));
        listmat(k,4)=str2double(filenamelist{k}(17:18));
        listmat(k,1)=Date2DOY( listmat(k,3) , listmat(k,4) , listmat(k,2) );
    elseif strcmp('met',filenamelist{k}(1:3))
        listmat(k,2)=str2double(filenamelist{k}(10:13));
        listmat(k,3)=str2double(filenamelist{k}(15:16));
        listmat(k,4)=str2double(filenamelist{k}(18:19));
        listmat(k,1)=Date2DOY( listmat(k,3) , listmat(k,4) , listmat(k,2) );
    else strcmp(tsOrmet,filenamelist{k}(1:namestr))
        listmat(k,2)=str2double(filenamelist{k}(namestr+1:namestr+3)); % get year from file name
        listmat(k,3)=str2double(filenamelist{k}(15:16)); % get month from file name
        listmat(k,4)=str2double(filenamelist{k}(18:19)); % get day from file name
        listmat(k,1)=Date2DOY( listmat(k,3) , listmat(k,4) , listmat(k,2) );
    end
    end
end

FrstD=listmat(1,1);
LastD=listmat(end,1);
LstHst=hist(listmat(:,1),LastD-FrstD+1);

filecnt=1;
fileCell=cell(LastD-FrstD+1,max(LstHst)+1);
for k=1:LastD-FrstD+1
    if LstHst(k)==0
        fileCell{k,1} = k-1+FrstD;
    elseif LstHst(k)==2 % 2 same doy
        fileCell{k,1} = k-1+FrstD;
        fileCell{k,2} = filenamelist{filecnt};
        fileCell{k,3} = filenamelist{filecnt+1};
        filecnt=filecnt+2;
    elseif LstHst(k)==1 
        fileCell{k,1} = k-1+FrstD;
        fileCell{k,2} = filenamelist{filecnt};
    	filecnt=filecnt+1;
    else
        fileCell{k,1} = k-1+FrstD;
        for i=1:LstHst(k)
            fileCell{k,i+1} = filenamelist{filecnt+i-1};
        end
        filecnt=filecnt+LstHst(k);
    end
end

[r,c]=size(fileCell);
fileCellStat=zeros(r,c);
filecode=nan(r,1);
for k=1:r
    for l=1:c
        fileCellStat(k,l)=isempty(fileCell{k,l});
    end
    filecode(k,1)=find(~fileCellStat(k,:),1,'last')-1;
end

% (find(~fileCellStat(CD,:),1,'last')-1) will give me the number of files
% with the same date that we need to patch together. 
% with this I'm generating met list and ts list and i need to figure for
% each day if I need to use merger of just read the file.dat or just
% prealocate for the day.

Day1=max(Day1,fileCell{1,1});
Dayend=min(Dayend,fileCell{end,1});

for i = 1:r
    if fileCell{i,1}==Day1
        Day1i=i;
    end
    if fileCell{i,1}==Dayend
        Dayendi=i;
    end
end
fileCell1=cell(Dayendi-Day1i+1,c);
%filecode=filecode(Dayendi-Day1i+1,1);
for k=1:Dayendi-Day1i+1
    for l=1:c
        fileCell1{k,l}=fileCell{k+Day1i-1,l};
    end
end

filecode1=filecode(Day1i:Dayendi);

end