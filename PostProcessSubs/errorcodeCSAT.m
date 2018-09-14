function errorcodeCSAT()

% The cell diagnostic value is a 1 byte unsigned integer (value between 0 and 255) 
% with the following bit map:
% bit 7   bit 6    bit 5 bit 4 bit 3                bit 2             bit 1                 bit 0
% Chopper Detector PLL   Sync  <------------------------ ----AGC / 6.25 ------------------------>
% 1=ok    1=ok     1=ok  1=ok


    ErrCode = input ('Insert error code  ');

    ErrCodeBin =dec2bin(ErrCode,16);
    
    if  (ErrCode)>= 128
        sprintf('%s','Chopper ok')
        ErrCode=ErrCode-128;
    else sprintf('%s','Chopper NOT ok')
    end
    if  (ErrCode)>=64 
        sprintf('%s','Detector ok')
        ErrCode=ErrCode-64;
    else sprintf('%s','Detector NOT ok')
    end
    if  (ErrCode)>=32 
        sprintf('%s','PLL (Phase Lock Loop) ok')
        ErrCode=ErrCode-32;
    else sprintf('%s','PLL (Phase Lock Loop) NOT ok')    
    end
    if  (ErrCode)>=16 
        sprintf('%s','Sync ok')
        ErrCode=ErrCode-16;
    else sprintf('%s','Sync NOT ok')
    end
    ErrCode = ErrCode*6.25;
    sprintf('AGC is %d . (The smaller the better. Typical values: 50 to 100)',ErrCode)
end