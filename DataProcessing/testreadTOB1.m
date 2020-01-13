% clear;close;clc
% metDatOut=[];
% for i=0:22
%     [metData,metHeaders,Var_names]=readTOB1(['C:\Users\thmorin\Desktop\testDrive\thmremote_met_data_2019_07_02_' ...
%         num2str(i*100,'%04d') '.dat']);
%     metDatOut=[metDatOut;metData];
% end
% 
% tsDatOut=[];
% for i=0:21
%    [tsData,tsHeaders,Var_names]=readTOB1(['C:\Users\veron\Documents\Research\Solvay\DataProcessing\SolvayRawData\thmremote_ts_data_2019_10_23_' ...
%        num2str(i*100+59,'%04d') '.dat']);
%    tsDatOut=[tsDatOut;tsData];
% end

function [data,Header,Var_names,year,month,day,hour,minut,secon]=testreadTOB1(name)
fileID=fopen(name);
for i=1:5
    Header{i}=fgetl(fileID);
end
Header{2}(Header{2}=='"')=[];
Var_names=regexp(Header{2},',','split');
cols=length(Var_names);
C=fread(fileID,[cols inf],'float32')';

fclose(fileID);

fileID=fopen(name);
for i=1:5
    Header{i}=fgetl(fileID);
end
Header{2}(Header{2}=='"')=[];
Var_names=regexp(Header{2},',','split');
cols=length(Var_names);
B =fread(fileID,[cols inf],'uint32')';

fclose(fileID);

fileID=fopen(name);
for i=1:5
    H{i}=fgetl(fileID);
end
D=fread(fileID,[cols inf],'uint32')';
fclose(fileID);

zip=zeros(size(D,1),1);

sec=round(D(:,1)+D(:,2)/10^9,1);

T_Stamp=datenum(1990,1,1,0,0,0)+datenum(zip,zip,zip,zip,zip,sec);
[year,month,day,hour,minut,secon]=datevec(T_Stamp);
secon=round(secon,1);
data=[T_Stamp D(:,3) C(:,4:12) B(:,13) C(:,14:end)];
end