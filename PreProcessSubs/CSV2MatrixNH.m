function [ConvertedT0A] = CSV2MatrixNH( FileToRead )

% function [ConvertedT0A] = CSV2Matrix(FileToRead)
% This function reads dayily data files in TOA5 format made by Baler software.
%
%  Example of gold file in CSV format : 
% 2002,4,14,0,0,0,0,0.11,-0.93,0.6,20.82,2.992,2.246,104,30,-2.571,-47.7,-0.015,0.002,16.91,76.1,17.63,99.6,0,0.148,3.261,0.071,0.108,0.643
% 2002,4,14,0,0,0,100,0.02,-0.97,0.63,20.87,2.993,2.251,104,30,-2.571,-47.7,-0.015,0.002,16.91,76.1,17.63,99.6,0,0.148,3.261,0.071,0.108,0.643
% 2002,4,14,0,0,0,200,-0.15,-1.07,0.57,20.94,2.992,2.256,104,30,-2.571,-47.7,-0.015,0.002,16.91,76.1,17.63,99.6,0,0.148,3.261,0.071,0.108,0.643
%   1  2  3 4 5 6  7    8     9    10   11    12    13    14 15    16     17     18   19    20   21    22  23   24 25     26     27     28    29 %
% The function coverts the date and timestamp to a matlab date number
% (datenum), [Date-> Julian date], [14:38:18.9 -> 14.638583333333335]
%                                   14+(38*60+18.9)/3600
% This function calls the function IsLeapYear.
% Output:
% ConvertedT0A = Data as numbers array.
% FirstLineInFile, SecondLineInFile, ThirdLineInFile, ForthLineInFile = file header.
% the data format is = '%d,%d,%d,%d,%f,%f,%f,%f,%f,%f,%f,%f,%f'.

% addpath( genpath( FileToRead));
fid = fopen( (FileToRead ) ,'r');    % FOPEN(FILENAME,PERMISSION) 
                                % opens the file FILENAME in the mode specified by PERMISSION:
                                % 'r'  open file for reading
                                
% FirstLineInFile = fgetl(fid) ;   % function fgetline (file) reads the next line until end-of-line

% SecondLineInFile = fgetl(fid) ;  % The second Line in the file gives the number of columns
% NOC = regexp(SecondLineInFile, '[^",]+','match') ;
NumberOfCol = 29 ;

% ThirdLineInFile = fgetl(fid) ;
% 
% ForthLineInFile = fgetl(fid);


% Read the rest of the file char by char and close the file
[C, NumElFile] = fread(fid, inf, 'char=>char');
fclose(fid);

% Quotes (") are replaced with spaces
% m = C == '"';
% C(m) = ' ';

% Commas (,) are replaced with spaces
m = C == ',';
C(m) = ' ';

% Generate the Format of the matrix %#ok<AGROW> supresses the remark of ineficiency

frmt = '%f %f %f %f %f %f';
for i = 7:NumberOfCol
    frmt = [frmt ' %f'];    %#ok<AGROW>
end

% Read formatted data from string A = sscanf(str, format, sizeA)
a = sscanf(C,frmt,[NumberOfCol,inf]);

% Convert first six columns to matlab date number
year = a(1,:) ;
[day , ~ ] = IsLeapYear (year); % Added function day is a vector of julian days of year 365/366
DecimalTime = a(4,:)+(a(5,:)*60 + a(6,:)+ a(7,:)/1000)/3600 ;
ConvertedT0A = [a(1,:)' day(a(2,:))'+a(3,:)'-1 DecimalTime' a(8:end,:)'];

% Prints out length
DataLength = length(ConvertedT0A);
fprintf('Read ( %d ) charachters for ( %d ) Number Of Rows\n',NumElFile,DataLength);

end