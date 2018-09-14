function errorcodeLI7700()
%    Diagnostics Header           Meaning
%                         Integer
% 1  BOXCONNECTED         1       LI-7550 Attached
% 2  BADAUXTC3            2       Bad thermocouple values
% 3  BADAUXTC2            4       Bad thermocouple values
% 4  BADAUXTC1            8       Bad thermocouple values
% 5  MOTORFAILURE         16      Mirror cleaner motor failure
% 6  CALIBRATING          32      Calibration in process
% 7  BOTTOMHEATERON       64      Lower mirror heater on
% 8  TOPHEATERON          128     Upper mirror heater on
% 9  PUMPON               256     Pump motor running
% 10 MOTORSPINNING        512     Mirror spin motor on
% 11 BLOCKTEMPUNREGULATED 1024    Block temperature unregulated
% 12 LASERTEMPUNREGULATED 2048    Laser cooler unregulated
% 13 BADTEMP              4096    Optical path thermocouple failure
% 14 REFUNLOCKED          8192    Reference methane signal not locked
% 15 NOSIGNAL             16384   No laser signal detected
% 16 NOTREADY             32768   Instrument is not ready

Err = input ('Insert error code  ');

ErrCode = dec2bin(Err,16);


if ErrCode(1)=='1' 
    sprintf('%s','16 NOTREADY             32768   Instrument is not ready')  % 1
end
if ErrCode(2) == '1' 
    sprintf('%s','15 NOSIGNAL             16384   No laser signal detected') % 2
end
if ErrCode(3) == '1' 
    sprintf('%s','14 REFUNLOCKED          8192    Reference methane signal not locked') % 3
end
if  ErrCode(4) == '1' 
    sprintf('%s','13 BADTEMP              4096    Optical path thermocouple failure') % 4
end
if  ErrCode(5) == '1' 
    sprintf('%s','12 LASERTEMPUNREGULATED 2048    Laser cooler unregulated') % 5
end
if  ErrCode(6) == '1'
    sprintf('%s','11 BLOCKTEMPUNREGULATED 1024    Block temperature unregulated') % 6
end
if  ErrCode(7) == '1'
    sprintf('%s','10 MOTORSPINNING        512     Mirror spin motor on') % 7
end
if  ErrCode(8) == '1'
    sprintf('%s','9  PUMPON               256     Pump motor running') % 8
end
if  ErrCode(9) == '1'
    sprintf('%s','8  TOPHEATERON          128     Upper mirror heater on') % 9
end
if  ErrCode(10) == '1'
    sprintf('%s','7  BOTTOMHEATERON       64      Lower mirror heater on') % 10
end
if  ErrCode(11) == '1'
    sprintf('%s','6  CALIBRATING          32      Calibration in process') % 11
end
if  ErrCode(12) == '1'
    sprintf('%s','5  MOTORFAILURE         16      Mirror cleaner motor failure') % 12
     % 8 13
     % 4 14
     % 2 15
     % 1 16
end

end