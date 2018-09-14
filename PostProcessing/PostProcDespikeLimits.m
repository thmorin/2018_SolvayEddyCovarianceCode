function [Dspk]=PostProcDespikeLimits(SeasonType)
Dspk=struct;
Dspk.STD=6;
Dspk.Interval=672;
Dspk.trimp=5;
switch SeasonType
    case 'Summer'
    Dspk.SWinbar.min=-10;
    Dspk.SWinbar.max=1200;
    
    Dspk.SWoutbar.min=-10;
    Dspk.SWoutbar.max=400;
    
    Dspk.LWinbar.min=100;
    Dspk.LWinbar.max=700;

    Dspk.LWoutbar.min=100;
    Dspk.LWoutbar.max=700;
    
    Dspk.pressure.min=0;
    Dspk.pressure.max=102000;
    
    Dspk.tair.min=-20;
    Dspk.tair.max=40;
   
    Dspk.net_sw.min=-20;
    Dspk.net_sw.max=1000;
      
    Dspk.net_lw.min=-200;
    Dspk.net_lw.max=40;
   
    Dspk.wind_speed_1.min=0;
    Dspk.wind_speed_1.max=10;
      
    Dspk.wind_speed_2.min=0;
    Dspk.wind_speed_2.max=10;
     
    Dspk.p_vaporSat_bar.min=0;
    Dspk.p_vaporSat_bar.max=8000;
     
    Dspk.ustar_1.min=0;
    Dspk.ustar_1.max=1.5;
     
    Dspk.ustar_2.min=0;
    Dspk.ustar_2.max=1.5;
     
    Dspk.wts.min=-.15;
    Dspk.wts.max=.4;
    
    Dspk.rho_cp.min=1100;
    Dspk.rho_cp.max=1400;
     
    Dspk.H.min=-150;
    Dspk.H.max=500;
    
    Dspk.LE.min=-150;
    Dspk.LE.max=1300;
     
    Dspk.RNET.min=-150;
    Dspk.RNET.max=1000;
     
    Dspk.p_vapor_bar.min=0;
    Dspk.p_vapor_bar.max=3500;
    
    Dspk.rH.min=0;
    Dspk.rH.max=100;
    
    Dspk.Resp.Interval=336;
    Dspk.Resp.min=0;
    Dspk.Resp.max=30;
    
    Dspk.GPP.Interval=336;
    Dspk.GPP.min=-55;
    Dspk.GPP.max=0;
    
    Dspk.L.min=-70000;
    Dspk.L.max=70000;
    
    Dspk.ww.min=0;
    Dspk.ww.max=2;
    
    Dspk.vv.min=-2;
    Dspk.vv.max=2;
    
    Dspk.SoilTemp.min=-50;
    Dspk.SoilTemp.max=100;
    Dspk.SoilTemp.Interval=100;
    
    Dspk.Methane.min=-0.2;
    Dspk.Methane.max=0.6;
    case 'Winter'
    Dspk.SWinbar.min=-10;
    Dspk.SWinbar.max=1200;
    
    Dspk.SWoutbar.min=-10;
    Dspk.SWoutbar.max=400;
    
    Dspk.LWinbar.min=100;
    Dspk.LWinbar.max=700;

    Dspk.LWoutbar.min=100;
    Dspk.LWoutbar.max=700;
    
    Dspk.pressure.min=0;
    Dspk.pressure.max=102000;
    
    Dspk.tair.min=-20;
    Dspk.tair.max=40;
   
    Dspk.net_sw.min=-20;
    Dspk.net_sw.max=1000;
      
    Dspk.net_lw.min=-200;
    Dspk.net_lw.max=40;
   
    Dspk.wind_speed_1.min=0;
    Dspk.wind_speed_1.max=10;
      
    Dspk.wind_speed_2.min=0;
    Dspk.wind_speed_2.max=10;
     
    Dspk.p_vaporSat_bar.min=0;
    Dspk.p_vaporSat_bar.max=8000;
     
    Dspk.ustar_1.min=0;
    Dspk.ustar_1.max=1.5;
     
    Dspk.ustar_2.min=0;
    Dspk.ustar_2.max=1.5;
     
    Dspk.wts.min=-.15;
    Dspk.wts.max=.4;
    
    Dspk.rho_cp.min=1100;
    Dspk.rho_cp.max=1400;
     
    Dspk.H.min=-150;
    Dspk.H.max=500;
    
    Dspk.LE.min=-150;
    Dspk.LE.max=1000;
     
    Dspk.RNET.min=-150;
    Dspk.RNET.max=1000;
     
    Dspk.p_vapor_bar.min=0;
    Dspk.p_vapor_bar.max=3500;
    
    Dspk.rH.min=0;
    Dspk.rH.max=100;
    
    Dspk.L.min=-50000;
    Dspk.L.max=40000;
    
    Dspk.ww.min=0;
    Dspk.ww.max=2;
    
    Dspk.vv.min=-2;
    Dspk.vv.max=2;
    
    Dspk.SoilTemp.min=-50;
    Dspk.SoilTemp.max=100;
    Dspk.SoilTemp.Interval=100;
    
    Dspk.Resp.Interval=336;
    Dspk.Resp.min=0;
    Dspk.Resp.max=20;
    
    Dspk.GPP.Interval=336;
    Dspk.GPP.min=-30;
    Dspk.GPP.max=0;
    
    Dspk.Methane.min=-0.1;
    Dspk.Methane.max=0.15;
    otherwise
    error('Season input to PostProcDespikeLimits must be either ''Summer'' or ''Winter''');
end
