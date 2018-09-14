function Dspk = DeSpikeLimits()
    %ID200 CO2 sensors
    Dspk.CO21STD = 4;
    Dspk.CO21mn = 50;
    Dspk.CO21mx = 700;
    Dspk.CO21tpr = 6;

    Dspk.CO2STD = 4;
    Dspk.CO2mn = 180;
    Dspk.CO2mx = 700;
    Dspk.CO2tpr = 4;

    Dspk.TSTD = 6;
    Dspk.Tmn = [-22,-22,-22,-10,0,0,5,10,5,-7,-10,-15]-5;
    Dspk.Tmx = [19,19,30,33,34,47,50,55,48,32,23,20]+10;
    Dspk.Ttpr = 4;

    % Humidity
    Dspk.HSTD = 6;
    Dspk.Hmn = 0;
    Dspk.Hmx = 100;
    Dspk.Htpr = 4;

    % Vapor Pressure
    Dspk.VPSTD = 6;
    Dspk.VPmn = 0;
    Dspk.VPmx = 10;
    Dspk.VPtpr = 4;

    % Pressure 
    Dspk.PSTD = 6;
    Dspk.Pmn = 50;
    Dspk.Pmx = 150;
    Dspk.Ptpr = 4;

    % Radiation
    Dspk.SWiSTD = 6; % Short Wave in
    Dspk.SWimn = -100;
    Dspk.SWimx = 1500;
    Dspk.SWitpr = 4;

    Dspk.SWoSTD = 6; % Short Wave out
    Dspk.SWomn = -150;
    Dspk.SWomx = 750;
    Dspk.SWotpr = 4;

    Dspk.LWiSTD = 6; % Long wave in
    Dspk.LWimn = -300;
    Dspk.LWmx = 800;
    Dspk.LWitpr = 4;

    Dspk.LWoSTD = 6; % Long wave out
    Dspk.LWomn = -300;
    Dspk.LWomx = 800;
    Dspk.LWotpr = 4;

    Dspk.TPSTD = 6; % Total PAR
    Dspk.TPmn = -200;
    Dspk.TPmx = 3000;
    Dspk.TPtpr = 4;

    Dspk.DPSTD = 6;% Diffuse PAR
    Dspk.DPmn = -200;
    Dspk.DPmx = 1600;
    Dspk.DPtpr = 4;

    % Soil Temperature
    Dspk.STSTD = 4;
    Dspk.STmn = [-35,-35,-25,-20,-10,-10,-10,-10,-15,-20,-25,-35];
    Dspk.STmx = [30,30,30,40,50,65,65,45,45,35,35,30];
    Dspk.STtpr = 4;

    % Vertical wind
    Dspk.VWSTD = 6;
    Dspk.VWmn = -30;
    Dspk.VWmx = 30;
    Dspk.VWtpr = 0;

    % Horizontal wind
    Dspk.HWSTD = 6;
    Dspk.HWmn = -20;
    Dspk.HWmx = 20;
    Dspk.HWtpr = 4;

    % Carbon density
    Dspk.CSTD = 6;
    Dspk.Cmn = 0;
    Dspk.Cmx = 40;
    Dspk.Ctpr = 4;

    % Water density
    Dspk.QSTD = 6;
    Dspk.Qmn = 0;
    Dspk.Qmx = [1000, 1300, 1450, 1200];
    Dspk.Qtpr = 4;

    % Methane density
    Dspk.MSTD = 6;
    Dspk.Mmn = 0;
    Dspk.Mmx = 0.3;
    Dspk.Mtpr = 0;
end
    
