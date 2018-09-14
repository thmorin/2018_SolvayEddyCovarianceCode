function [ConvertedT0A, FirstLineInFile, SecondLineInFile, ThirdLineInFile, ForthLineInFile ] = T0A2Matrix( FileToRead )

% function [ConvertedT0A] = T0A2Matrix(FileToRead)
% This function reads dayily data files in TOA5 format made by Baler software.
%
%  Example of ts file in TOA5 format : 
% "TOA5","CR3000Olen","CR3000","3257","CR3000.Std.09.02","CPU:2011March19OlenWetland.CR3","25678","ts_data"
% "TIMESTAMP","RECORD","UxRMY","UyRMY","UzRMY","sosRMY","diag_RMY","UxCSAT","UyCSAT","UzCSAT","TsCSAT","diag_CSAT","LI7500_CO2","LI7500_H2O","press_LI7500","diag_LI7500","LI7700Diag","LI7700_CH4_dens","LI7700_Press","RSSI","checksum_datalogger"
% "TS","RN","m/s","m/s","m/s","m/s","unitless","m/s","m/s","m/s","Deg C","unitless","mmol/m^3","mmol/m^3","kPa","unitless","","mmol/m^3","kPa","%",""
% "","","Smp","Smp","Smp","Smp","Smp","Smp","Smp","Smp","Smp","Smp","Smp","Smp","Smp","Smp","Smp","Smp","Smp","Smp","Smp"
% "2011-03-23 00:00:00",2737275,"NAN","NAN","NAN","NAN","NAN",1.8625,-2.10325,0.09625001,8.961609,53,17.03966,351.2792,98.13717,248,142,0.0623619,97.9405,27.7657,65
% "2011-03-23 00:00:00.1",2737276,"NAN","NAN","NAN","NAN","NAN",1.983,-2.20475,0.118,8.968323,54,17.04043,351.1445,98.11115,248,142,0.0622026,97.9416,27.752,16
% "2011-03-23 00:00:00.2",2737277,"NAN","NAN","NAN","NAN","NAN",1.71875,-2.159,0.18925,8.944824,55,17.03819,351.2799,98.11115,248,142,0.0624314,97.9387,27.7902,15
% 
%  Example of met file in TOA5 format :
% "TOA5","CR3000Olen","CR3000","3257","CR3000.Std.09.02","CPU:2011March14OlenWetland.CR3","24962","met_data"
% "TIMESTAMP","RECORD","t_hmp_Avg","rh_hmp_Avg","e_hmp_Avg","batt_volt_Avg","panel_temp_Avg","SR01Up_Avg","SR01Dn_Avg","IR01Up_Avg","IR01Dn_Avg","TC_Avg","CO2at7m_Avg","CO2at8m_Avg","CO2at9m_Avg"
% "TS","RN","Deg C","%","kPa","V","Deg C","W/m2","W/m2","W/m2","W/m2","Deg C","ppm","ppm","ppm"
% "","","Avg","Avg","Avg","Avg","Avg","Avg","Avg","Avg","Avg","Avg","Avg","Avg","Avg"
% "2011-03-14 14:43:00",139,9.081827,39.71505,0.4581417,13.28092,14.36889,443.7727,74.31447,-64.76003,20.32286,10.72237,17.75363,17.45917,17.18596
% "2011-03-14 14:44:00",140,8.86272,40.51102,0.4604701,13.21117,14.37775,453.9491,73.57978,-64.74571,19.6335,10.65454,17.74778,17.45457,17.18303
% "2011-03-14 14:45:00",141,8.797,40.32148,0.4563026,13.27981,14.38661,459.9071,74.34117,-63.14458,20.63294,10.6178,17.76464,17.4699,17.1967
% 
%
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
                                
FirstLineInFile = fgetl(fid) ;   % function fgetline (file) reads the next line until end-of-line

SecondLineInFile = fgetl(fid) ;  % The second Line in the file gives the number of columns
NOC = regexp(SecondLineInFile, '[^",]+','match') ;
NumberOfCol = length(NOC) + 5 ;

ThirdLineInFile = fgetl(fid) ;

ForthLineInFile = fgetl(fid);


% Read the rest of the file char by char and close the file
[C, NumElFile] = fread(fid, inf, 'char=>char');
fclose(fid);

% Quotes (") are replaced with spaces
m = C == '"';
C(m) = ' ';

% Commas (,) are replaced with spaces
m = C == ',';
C(m) = ' ';

% Generate the Format of the matrix %#ok<AGROW> supresses the remark of ineficiency

frmt = '%d-%d-%d %d:%d:%f';
for i = 7:NumberOfCol
    frmt = [frmt ' %f'];    %#ok<AGROW>
end

% Read formatted data from string A = sscanf(str, format, sizeA)
a = sscanf(C,frmt,[NumberOfCol,inf]);

% Convert first six columns to matlab date number
year = a(1,:) ;
[day , ~ ] = IsLeapYear (year); % Added function day is a vector of julian days of year 365/366
DecimalTime = a(4,:)+(a(5,:)*60 + a(6,:))/3600;
ConvertedT0A = [a(1,:)' day(a(2,:))'+a(3,:)'-1 DecimalTime' a(8:end,:)'];

% Prints out length
DataLength = length(ConvertedT0A);
fprintf('Read ( %d ) characters for ( %d ) Number Of Rows\n',NumElFile,DataLength);

end