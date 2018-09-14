%function [ConvertedT0A, FirstLineInFile, SecondLineInFile, ThirdLineInFile, ForthLineInFile ] = T0A2Matrix( FileToRead )
function [ConvertedT0A, FirstLineInFile, SecondLineInFile] = T0A2Matrix2lH( FileToRead )

% function [ConvertedT0A] = T0A2Matrix(FileToRead,year)
% This function reads a TOA file recorded by a Campbell data logger
% (FileToRead).
%
%"TOA5","CR1000","CR1000","26208","CR1000.Std.16","CPU:wetland_program_Gil_May_6_2010.CR1","39760","met_data"
% "TIMESTAMP","RECORD","t_hmp_Avg","rh_hmp_Avg","e_hmp_Avg","batt_volt_Avg","panel_temp_Avg","SR01Up_Avg","SR01Dn_Avg","IR01Up_Avg","IR01Dn_Avg","TC_Avg","TCC_Avg"
% "TS","RN","C","percent","kPa","V","C","W/m2","W/m2","W/m2","W/m2","Deg C","Deg C"
% "","","Avg","Avg","Avg","Avg","Avg","Avg","Avg","Avg","Avg","Avg","Avg"
% "2010-07-14 10:54:00",74391,26.58174,70.32069,2.445851,14.41931,25.72323,766.2126,120.6859,-9.825724,-17.41301,1.127823,32.86842
% "2010-07-14 10:55:00",74392,26.38912,70.09409,2.410321,14.38056,25.79929,776.91,121.7151,-9.825727,-17.58474,1.12777,32.85489
% "2010-07-14 10:56:00",74393,26.46277,71.17045,2.457966,14.42023,25.87535,778.0906,121.9379,-9.36185,-17.72212,1.127937,32.89807
%
% The function coverts the time to a matlab date number (datenum),
% This function calls the function IsLeapYear.
% Output:
% ConvertedT0A = data as a matrix
% FirstLineInFile, SecondLineInFile, ThirdLineInFile, ForthLineInFile = file header.
% lab = '%d,%d,%d,%d,%f,%f,%f,%f,%f,%f,%f,%f,%f'.

addpath( genpath( char(FileToRead(1))))
fid = fopen( char( FileToRead(2)) ,'r');    % FOPEN(FILENAME,PERMISSION) 
                                % opens the file FILENAME in the mode specified by PERMISSION:
                                % 'r'  open file for reading
                                
FirstLineInFile = fgetl(fid) ;   % function fgetline (file) reads the next line until break
%Header(1) = textscan(FirstLineInFile, '%s') ;
SecondLineInFile = fgetl(fid) ;  % The second Line in the file gives the number of columns
%Header(2) = SecondLineInFile ;                                %
NOC = regexp(SecondLineInFile, '[^",]+','match') ;
NumberOfCol = length(NOC) + 5 ;
%ThirdLineInFile = fgetl(fid) ;
%Header(3) = ThirdLineInFile ;
%ForthLineInFile = fgetl(fid);
%Header(4) = ForthLineInFile ;

%read the rest of the file as char and close the file
%[A, count] = fread(fileID, sizeA-to the end, precision)
[C, NumElFile] = fread(fid, inf, 'char=>char');
fclose(fid);

%replace all quotes (") with space
m = C == '"';
C(m) = ' ';

%Find dashes (-) that occur directly after a number (in the date),
% colons (:) and commas (,) and replace them with spaces
m = C == ',';
C(m) = ' ';

%replace all 'NAN' with 7777
% m1 = C == 777 ;
% m = C == 'N';
% C(m) = '7';
% m = C == 'A';
% C(m) = '7';
% 
% Was777 = size(m1);

%Only numbers remain, read them all in the correct order of coulmns

%Generatr the Format of the matrix
frmt = '%d-%d-%d %d:%d:%f';
for i = 7:NumberOfCol
    frmt = [frmt ' %f'];
end

%Read formatted data from string A = sscanf(str, format, sizeA)
a = sscanf(C,frmt,[NumberOfCol,inf]);

%convert first six columns to matlab date number

year = a(1,:) ;
[day , ~ ] = IsLeapYear (year); % Added function day is a vector of julian days of year 365/366
DecimalTime = a(4,:)+(a(5,:)*60 + a(6,:))/3600;
ConvertedT0A = [a(1,:)' day(a(2,:))'+a(3,:)'-1 DecimalTime' a(8:end,:)'];

%replace fill values (777) with nan
%ConvertedT0A(ConvertedT0A==777) = nan;

%some stats
DataLength = length(ConvertedT0A);
fprintf('Read %d bytes for %d Number Of Columns: %g bytes/row \n',NumElFile,DataLength,NumElFile/DataLength);
fprintf('Number of original 777s that where replaced by nan is %d \n',Was777);

end