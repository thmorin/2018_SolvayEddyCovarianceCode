function NEE_Season(NEE,DECDAY,ust)
day1=1;
day2=floor(length(NEE)/48);
NEE(ust<0.25)=nan;
nee=DeSpike(NEE,672,5,-25,15,10);
[~, nee1]=gapfill(nee,(day2-day1+1),48);
a=1; b=ones(1,7*48)/7/48;
NEE_Ave=filtfilt(b,a,nee1);

figure(2);clf
plot(DECDAY,NEE_Ave);
hold on
grid on
j=0;
plot(DECDAY,j,'r');
xlabel('Day of year');
ylabel('NEE (umol/m2/s)');
xlim([min(DECDAY) max(DECDAY)]);
hold off

end
