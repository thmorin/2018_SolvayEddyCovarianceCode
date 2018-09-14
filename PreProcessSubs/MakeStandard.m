function [StandardData] = MakeStandard( RawData , year , currentday , MeasFreq ,HrsAFile)

% Takes as input the data from T0A2Matrix , year , currentday , frequency
% (MeasFreq) and Number of hours per file (HrsAFile)
% (for 1 measurment per min MeasFreq is 1/60), and returns a matrix with
% ALL the timestamps of the day.

ntimes = HrsAFile*3600*MeasFreq;               % Number of time stamps per day.

Hr0 = floor(RawData(1,3));
Hrend = min(Hr0+HrsAFile,24);
Hrstart = max(0, Hrend-HrsAFile);
timestep = 1/(3600*MeasFreq);
t32 = int32(RawData(:,3)/timestep);
t0 = find(t32==0);
if ~isempty(t0)
    if RawData(t0(end),2) == currentday;
        t32 = t32 +1;
        TimeVector = [currentday*ones(ntimes,1), linspace(Hrstart,Hrend-timestep,ntimes)'];
    else
        t32(t32==0)=24/timestep;
        TimeVector = [currentday*ones(ntimes,1), linspace(Hrstart+timestep,Hrend,ntimes)'];
    end
else
    TimeVector = [currentday*ones(ntimes,1), linspace(Hrstart+timestep,Hrend,ntimes)'];
end
if Hrend ==24
    TimeVector(end,1) = TimeVector(end,1)+1;
    TimeVector(end,2) = 0;
end
    
nc = length(RawData(1,:));
StandardData = nan(ntimes,nc);
StandardData( t32, 4:nc ) = RawData(:,4:nc);
StandardData( 1:ntimes, 1:3 ) = [ year*ones(ntimes,1) , TimeVector]; %creates time columns 

end