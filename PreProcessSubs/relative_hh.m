function ind=relative_hh(Strt_date_num,date_num)
% add 1 to ind and that's the index in the year vector
[y_s,m_s,d_s,h_s,mi_s,~] = datevec(Strt_date_num);
[y_c,m_c,d_c,h_c,mi_c,~] = datevec(date_num);


ind=(datenum(y_c,m_c,d_c)-datenum(y_s,m_s,d_s))*48 +(h_c-h_s)*2+round((mi_c-mi_s)./30);
end
% xdate=(datenum(2011,04,01)):1/48:(datenum(2011,04,02)+1);
% xdatem=floor((xdate-Strt_date_num)*48);
% X=strt_11_04-(strt_11_04-end_11_04)/2;
% [~,c]=size(CH4_c_11_04);
% Xc=reshape(CH4_c_11_04,2,c/2);
% Y=nanmean(Xc);
% plot (xdate,...
%     umol_2_grC*meth_gen_gf(xdatem)./(wer_per(xdatem)),...
%     'LineWidth',2,'Color', [0.847058832645416 0.160784319043159 0]);hold on;
% %xlim([datenum(2011,03,19,12,0,0)  datenum(2013,03,20,11,30,0)]);
% datetick('x','mm/dd HH:MM','keepticks')
% % Create labels and title
% xlabel({''},'Interpreter','Tex','FontSize',18);
% ylabel({'CH_4 flux [ grCH4-C m^-^2 hr^-^1] '},'Interpreter','Tex','FontSize',18);
% title({'Gap-filed Methane vs. chamber measurement'},'Interpreter','Tex','FontSize',20);
% % Open water 1,3,5,7
% bar1=bar(strt_11_04(1)-(strt_11_04(1)-end_11_04(1))/2,...
%     Y(1),(strt_11_04(1)-end_11_04(1)));
% baseline1 = get(bar1,'BaseLine');
% set(baseline1,'LineWidth',2,'LineStyle','--',...
%     'Color',[0.313725501298904 0.313725501298904 0.313725501298904]);
% bar(strt_11_04(3)-(strt_11_04(3)-end_11_04(3))/2,...
%     Y(3),(strt_11_04(3)-end_11_04(3)));
% bar(strt_11_04(5)-(strt_11_04(5)-end_11_04(5))/2,...
%     Y(5),(strt_11_04(5)-end_11_04(5)));
% bar(strt_11_04(7)-(strt_11_04(7)-end_11_04(7))/2,...
%     Y(7),(strt_11_04(7)-end_11_04(7)));
% % IM/Veg 2,4,6,8
% bar(strt_11_04(2)-(strt_11_04(2)-end_11_04(2))/2,...
%     Y(2),(strt_11_04(2)-end_11_04(2)),'FaceColor',[0 0.498039215803146 0]);
% bar(strt_11_04(4)-(strt_11_04(4)-end_11_04(4))/2,...
%     Y(4),(strt_11_04(4)-end_11_04(4)),'FaceColor',[0 0.498039215803146 0]);
% bar(strt_11_04(6)-(strt_11_04(6)-end_11_04(6))/2,...
%     Y(6),(strt_11_04(6)-end_11_04(6)),'FaceColor',[0 0.498039215803146 0]);
% bar(strt_11_04(8)-(strt_11_04(8)-end_11_04(8))/2,...
%     Y(8),(strt_11_04(8)-end_11_04(8)),'FaceColor',[0 0.498039215803146 0]);
% % Upland 9
% bar(strt_11_04(9)-(strt_11_04(9)-end_11_04(9))/2,...
%     Y(9),(strt_11_04(9)-end_11_04(9)),'FaceColor',[0.682352960109711 0.466666668653488 0]);
% 
% minXc=min(Xc);maxXc=max(Xc);
% errorbar(X,Y,Y-minXc,maxXc-Y,'LineStyle','none','Color',[0 0 1]);