function [ustk,ustfilt]=GPF3(Fc,spikeflag,Ta,ustar,daynight) 

%Updated June 2013 by T. Morin
%Updated May 2010 by K. Maurer and G. Bohrer

Fcgpf=ones(size(Fc));
Fcgpf(spikeflag==0|Fc==0|Fc<-100|ustar<0|isnan(Fc)|Fc<0)=0;

ustk = zeros(1,6);

Tp=Ta(Fcgpf==1 & daynight==0 & ~isnan(Ta) & ~isnan(ustar));
ustarp=ustar(Fcgpf==1 & daynight==0 & ~isnan(Ta) & ~isnan(ustar));
Fcp=Fc(Fcgpf==1 & daynight==0 & ~isnan(Ta) & ~isnan(ustar));

ll=length(Tp);
Tgaps=sort(Tp);
Tc=Tgaps(1);
for yy=1:5
    Tc(yy+1)=Tgaps(yy*(floor(ll/6)));
end
Tc(7)=Tgaps(ll);
%Tc=min(Tp):(max(Tp)-min(Tp))/6:max(Tp);

uo = zeros(1,6);

if ~isempty(Tc)
    mm=nan(6,1);

    for i=1:6
        Tclass=Tp(Tp>=Tc(i) & Tp<=Tc(i+1));
        uclass=ustarp(Tp>=Tc(i) & Tp<=Tc(i+1));
        nanTclass = Tclass(isnan(Tclass)==0 & isnan(uclass)==0);
        nanuclass = uclass(isnan(Tclass)==0 & isnan(uclass)==0);
        xx=corrcoef(nanTclass,nanuclass);
        if isnan(xx(1))==0 && length(Tclass)>2
            x=abs(xx(1,2));
        else
            x=1;
        end
        if x<.4
            mm(i)=length(uclass);
            Ugaps=sort(uclass);
            uc=Ugaps(1);
            for zz=1:19 %Should be zz=1:18 by paper?
                uc(zz+1)=Ugaps(zz*(floor(mm(i)/20)));
            end
            uc(21)=Ugaps(mm(i)); %Should be uc(20) by paper?
            %uc=(min(uclass):(max(uclass)-min(uclass))/20:max(uclass));
            cnt=1;
            RmRx=0;
            while RmRx<.95 && cnt<20 %Should be 0.95 by paper?
                Rm=nanmean(Fcp(ustarp>=uc(cnt)& ustarp<uc(cnt+1)));
                Rx=nanmean(Fcp(ustarp>=uc(cnt+1)));
                RmRx=Rm/Rx;
                uox=nanmean(ustarp(ustarp>=uc(cnt)& ustarp<uc(cnt+1)));
                cnt=cnt+1;
                %disp(['For class: ' num2str(cnt) ' RmRx was: ' num2str(RmRx), 'and uox: ' num2str(uox)]);
            end
            uo(i)=uox;
        else
            uo(i)=0;
        end
    end

    ustk(1,:)= uo;
end
ustfilt = nanmedian(ustk(1,:),2);

end 


