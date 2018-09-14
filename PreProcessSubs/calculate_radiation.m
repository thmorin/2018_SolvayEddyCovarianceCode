function [NetRs, NetRl, Albedo, UpTot, DownTot, NetRad, IRUp, IRDown, SR01Up, SR01Dn, QAstruct ] = calculate_radiation(SR01Up,SR01Dn,IR01Up,IR01Dn, TCC,Dspk )

%NR01TK is despiked outside in main process.
% Up - sensor facing up, radiation to earth (from space)
% Down - sensor facing down, radiation from earth (to space)

NR01TK=TCC+273.15;
IRUp=IR01Up+(567E-10).*NR01TK.^4;
IRDown=IR01Dn+(567E-10).*NR01TK.^4;

[SR01Up, SR01UpQA] = DeSpike(SR01Up,10,Dspk.SWiSTD,Dspk.SWimn,Dspk.SWimx,Dspk.SWitpr );  % NR01 SR01Up [ W/m^2 ]
QAstruct.SR01UpQA=SR01UpQA;
[SR01Dn, SR01DnQA] = DeSpike(SR01Dn,10,Dspk.SWoSTD,-50,550,Dspk.SWotpr );    % NR01 SR01Dn [ W/m^2 ]
QAstruct.SR01DnQA=SR01DnQA;
[IRUp, IRUpQA] = DeSpike(IRUp,10,6,-200,600,Dspk.LWitpr );       % NR01 IR01Up [ W/m^2 ]
QAstruct.IRUpQA=IRUpQA;
[IRDown, IRDownQA] = DeSpike(IRDown,10,6,-200,600,Dspk.LWotpr );   % NR01 IR01Dn [ W/m^2 ]
QAstruct.IRDownQA=IRDownQA;

NetRs = SR01Up - SR01Dn ;
NetRl = IRUp - IRDown ;
Albedo = SR01Dn./SR01Up ;
Albedo(Albedo<0 | Albedo>1) = nan;
UpTot = SR01Up + IRUp ;
DownTot = SR01Dn + IRDown ;
NetRad = UpTot - DownTot ;

end