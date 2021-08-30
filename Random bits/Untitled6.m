Fighandle=figure;
set(Fighandle, 'Position', [10, 10, 1000, 500]);
delay=10;
ha=tight_subplot(2,5);
axes(ha(1));
a=area(Speed_flow(1,49-delay:176),'FaceColor','m');a.FaceAlpha = 0.3;a.LineStyle='None';a.BaseLine.LineStyle='None';set(gca,'visible','off');hold on;
plot(Int_data_temp_mean(49+5:176)-min(Int_data_temp_mean(49:176)),'k','LineWidth',4);axis([0 220 -0.5 3.5]);
axes(ha(2));
a=area(Speed_flow(1,1100-delay:1272),'FaceColor','m');a.FaceAlpha = 0.3;a.LineStyle='None';a.BaseLine.LineStyle='None';set(gca,'visible','off');hold on;
plot(Int_data_temp_mean(1100+5:1272)-min(Int_data_temp_mean(1100+5:1272)),'k','LineWidth',4);axis([0 220 -0.5 3.5]);set(gca,'visible','off');
axes(ha(3));
a=area(Speed_flow(1,1000-delay:1100),'FaceColor','m');a.FaceAlpha = 0.3;a.LineStyle='None';a.BaseLine.LineStyle='None';set(gca,'visible','off');hold on;
plot(Int_data_temp_mean(1000+5:1120)-min(Int_data_temp_mean(1000+5:1120)),'k','LineWidth',4);axis([0 220 -0.5 3.5]);set(gca,'visible','off');
axes(ha(4));
a=area(Speed_flow(1,249-delay:440),'FaceColor','m');a.FaceAlpha = 0.3;a.LineStyle='None';a.BaseLine.LineStyle='None';set(gca,'visible','off');hold on;
plot(Int_data_temp_mean(249+5:440)-min(Int_data_temp_mean(249+5:440)),'k','LineWidth',4);axis([0 220 -0.5 3.5]);set(gca,'visible','off');
axes(ha(5));
a=area(Speed_flow(1,549-delay:800),'FaceColor','m');a.FaceAlpha = 0.3;a.LineStyle='None';a.BaseLine.LineStyle='None';set(gca,'visible','off');hold on;
plot(Int_data_temp_mean(549+5:800)-min(Int_data_temp_mean(549+5:800)),'k','LineWidth',4);axis([0 220 -0.5 3.5]);set(gca,'visible','off');

axes(ha(6));i=1;bar(nanmean(int_quant_perFish(:,i),1),'FaceColor',[0.5 0.5 0.5]);axis([0.5 2 0 275]);hold on;er = errorbar(nanmean(int_quant_perFish(:,i),1),nanstd(int_quant_perFish(:,i),1,1));er.Color = [0 0 0];er.LineWidth = 2;hold off;set(gca,'box','off');set(gca,'color','none');set(gca,'xtick',[])
axes(ha(7));i=5;bar(nanmean(int_quant_perFish(:,i),1),'FaceColor',[0.5 0.5 0.5]);axis([0.5 2 0 275]);hold on;er = errorbar(nanmean(int_quant_perFish(:,i),1),nanstd(int_quant_perFish(:,i),1,1));er.Color = [0 0 0];er.LineWidth = 2;hold off;set(gca,'visible','off');
axes(ha(8));i=4;bar(nanmean(int_quant_perFish(:,i),1),'FaceColor',[0.5 0.5 0.5]);axis([0.5 2 0 275]);hold on;er = errorbar(nanmean(int_quant_perFish(:,i),1),nanstd(int_quant_perFish(:,i),1,1));er.Color = [0 0 0];er.LineWidth = 2;hold off;set(gca,'visible','off');
axes(ha(9));i=2;bar(nanmean(int_quant_perFish(:,i),1),'FaceColor',[0.5 0.5 0.5]);axis([0.5 2 0 275]);hold on;er = errorbar(nanmean(int_quant_perFish(:,i),1),nanstd(int_quant_perFish(:,i),1,1));er.Color = [0 0 0];er.LineWidth = 2;hold off;set(gca,'visible','off');
axes(ha(10));i=3;bar(nanmean(int_quant_perFish(:,i),1),'FaceColor',[0.5 0.5 0.5]);axis([0.5 2 0 275]);hold on;er = errorbar(nanmean(int_quant_perFish(:,i),1),nanstd(int_quant_perFish(:,i),1,1));er.Color = [0 0 0];er.LineWidth = 2;hold off;set(gca,'visible','off');
