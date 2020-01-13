function instrumentDiagnostic = diagnosticParameters(instrument , diagnostic ,LI7700RSSI, RSSIper)
    % 1 = good
    % 0 = bad

    if strcmp('NONE', instrument) 
        instrumentDiagnostic = ones(size(diagnostic));
    elseif strcmp('CSAT3', instrument)
        instrumentDiagnostic = ((diagnostic>=0 & diagnostic<=63) | (diagnostic>=4032 & diagnostic<=4095)); %logical variable =1 when diagnostic within good bounds, =0 otherwise
    elseif strcmp('CSAT3B',instrument)
        instrumentDiagnostic = diagnostic<=64; %VERONICA FOLLOW UP ON THIS
    elseif strcmp('RMYOUNG', instrument)
        instrumentDiagnostic = (diagnostic==0); %logical variable =1 when diagnostic within good bounds, =0 otherwise

    elseif strcmp('LI7500', instrument)
        diagnostic(isnan(diagnostic))=99;
        D=dec2bin(diagnostic,8);
        Working=( (sum((D(:,1:4)=='0'),2)==0) & (bin2dec(num2str(D(:,5:8))) < 14) );
        instrumentDiagnostic = (Working==1); %logical variable =1 when diagnostic within good bounds, =0 otherwise
    elseif strcmp('LI7500A',instrument) %<--is ours actually a LI7500 or an LI7500A?
        instrumentDiagnostic=diagnostic>250; %VERONICA PLEAE CHECK
    elseif strcmp('LI7700', instrument)
        diagnostic( diagnostic>(2^16-1) | diagnostic<0 | isnan(diagnostic) ,1) = 2^16-1;
        D=dec2bin(diagnostic,16);
        instrumentDiagnostic = ~(D(:,1) == '1' | D(:,3) == '1' | D(:,5) == '1');
        if nargin==4
            instrumentDiagnostic(LI7700RSSI<RSSIper | isnan(LI7700RSSI) ) = 0;
        end
        
    elseif strcmp('EC150', instrument)
        if nargin==2
            instrumentDiagnostic = (diagnostic==0);        
        else
            instrumentDiagnostic = (diagnostic==0 & RSSIper>.05);
        end 
        
        
        % All the codes are at errorcodeLI7700.m
        % D(1) = 32768   Instrument is not ready
        % D(3) = 8192    Reference methane signal not locked
        % D(5) = 2048    Laser cooler unregulated
    end
end