function [Data, Header] =merge2datasets(directory,year, DOY, Filenew, FileOld, diffrent_table,type)

% Takes 2 .dat files and give a F/SData table.
% If it's due to different structure than diffrent_table='1' if it's the
% same structure than diffrent_table='0' and the MergDayStruct_year_DOY is
% generated here. The output will go to make standart and there the missing
% lines will fill with nan and the timestamps will be ordered.
%
% Fcolnew=[1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22]';
% Fcolold=[1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 nan nan nan nan 22]';
% the way to do this file is to put nan in a column that should be
% eliminated and the number of the column the variable should go to in
% Fcolnew. in this example columns 18,19,20 from the early data should be
% eliminated.
% Fcolnew=[1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22]';
% Fcolold=[1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 22]';
% in this example old column 18 will go into column 22 in the new Data
% matrix.


[DataNew, Header.F1, Header.F2, Header.F3, Header.F4 ] = T0A2Matrix( [directory Filenew] );
[DataOld, ~, ~, ~, ~ ] = T0A2Matrix( [directory FileOld] );

if diffrent_table==1
    eval(['MergDayStruct_' num2str(year) '_' num2str(DOY,'%03.0f') type ';']);
else
    [~, c]=size(DataNew);
    Fcolnew=1:c;
    Fcolold=Fcolnew;
end

[FDlastrowNew,FDlastcolNew]=size(DataNew);
[FDlastrowOld,FDlastcolOld]=size(DataOld);
Data=nan(FDlastrowNew+FDlastrowOld,max(FDlastcolNew,FDlastcolOld));
Data(1:FDlastrowNew,1:FDlastcolNew)=DataNew;

for k=1:min(length(Fcolnew),length(Fcolold))
    if ~isnan(Fcolold(k))
        Data(FDlastrowNew+1:FDlastrowNew+FDlastrowOld,k)=DataOld(:,k);
    end
end

Data=Data(:,~isnan(Fcolold));

end
