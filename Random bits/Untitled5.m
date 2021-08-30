Int_data_temp=ZS2(Test_Int_idx,:);
[Int_idxKmeans Int_Cmap]=kmeans(Int_data_temp,2,'Options',options,'Distance','correlation','Replicates',3,'MaxIter',1000,'Display','final');
figure;plot(Int_Cmap');hold on;plot(Speed_flow(1,:)/10-0.1)

%% ble
int_quant=nan(length(Fish_list),7);
int_toplot=nan(length(Fish_list),size(ZS2,2));
Int_fish_temp=idx_Fish(Test_Int_idx);
good_cluster=2;
delay=-5;
Nb_neurons=5;
for fish_nb=1:length(Fish_list)
    idx_fish_temp=find(Int_fish_temp==Fish_list(fish_nb));        
    idx_fish_temp=intersect(find(Int_idxKmeans==good_cluster),idx_fish_temp);
    if length(idx_fish_temp)>Nb_neurons                
        Int_data_temp_mean=mean(Int_data_temp(idx_fish_temp,:),1);
        int_toplot(fish_nb,:)=Int_data_temp_mean;
        int_quant(fish_nb,1)=trapz(Int_data_temp_mean(49-delay:176)-min(Int_data_temp_mean(49:176)));
        int_quant(fish_nb,2)=trapz(Int_data_temp_mean(249-delay:440)-min(Int_data_temp_mean(249:440)));
        int_quant(fish_nb,3)=trapz(Int_data_temp_mean(549-delay:800)-min(Int_data_temp_mean(549:800)));
        int_quant(fish_nb,4)=trapz(Int_data_temp_mean(1000-delay:1120)-min(Int_data_temp_mean(1000:1120)));
        int_quant(fish_nb,5)=trapz(Int_data_temp_mean(1120-delay:1272)-min(Int_data_temp_mean(1120:1272)));
        int_quant(fish_nb,6)=trapz(Int_data_temp_mean(1825-delay:2005)-min(Int_data_temp_mean(1825:2005)));
        int_quant(fish_nb,7)=trapz(Int_data_temp_mean(2080-delay:2245)-min(Int_data_temp_mean(2080:2245)));
        int_quant(fish_nb,8)=length(idx_fish_temp>Nb_neurons);
        int_quant(fish_nb,9)=trapz(Int_data_temp_mean(1473-delay:1575)-min(Int_data_temp_mean(1473:1575)));
    end
end

int_quant_perFish=nan(size(int_quant,1)/3,9);
for i=0:5
    int_quant_perFish(i+1,:)=nanmean(int_quant(1+i*3:3+i*3,1:9),1);
end
int_quant_perFish(:,4)=(int_quant_perFish(:,4)+int_quant_perFish(:,9))/2;

Int_data_temp_mean=mean(Int_data_temp(Int_idxKmeans==good_cluster,:),1);
Fighandle=figure;
set(Fighandle, 'Position', [10, 10, 1000, 500]);
delay=10;
ha=tight_subplot(2,7);
axes(ha(1));
a=area(Speed_flow(1,49-delay:176),'FaceColor','m');a.FaceAlpha = 0.3;a.LineStyle='None';a.BaseLine.LineStyle='None';set(gca,'visible','off');hold on;
plot(Int_data_temp_mean(49+5:176)-min(Int_data_temp_mean(49:176)),'k','LineWidth',4);axis([0 220 -0.5 3.5]);
axes(ha(5));
a=area(Speed_flow(1,1100-delay:1272),'FaceColor','m');a.FaceAlpha = 0.3;a.LineStyle='None';a.BaseLine.LineStyle='None';set(gca,'visible','off');hold on;
plot(Int_data_temp_mean(1100+5:1272)-min(Int_data_temp_mean(1100+5:1272)),'k','LineWidth',4);axis([0 220 -0.5 3.5]);set(gca,'visible','off');
axes(ha(2));
a=area(Speed_flow(1,249-delay:440),'FaceColor','m');a.FaceAlpha = 0.3;a.LineStyle='None';a.BaseLine.LineStyle='None';set(gca,'visible','off');hold on;
plot(Int_data_temp_mean(249+5:440)-min(Int_data_temp_mean(249+5:440)),'k','LineWidth',4);axis([0 220 -0.5 3.5]);set(gca,'visible','off');
axes(ha(3));
a=area(Speed_flow(1,549-delay:800),'FaceColor','m');a.FaceAlpha = 0.3;a.LineStyle='None';a.BaseLine.LineStyle='None';set(gca,'visible','off');hold on;
plot(Int_data_temp_mean(549+5:800)-min(Int_data_temp_mean(549+5:800)),'k','LineWidth',4);axis([0 220 -0.5 3.5]);set(gca,'visible','off');
axes(ha(4));
a=area(Speed_flow(1,1008-delay:1100),'FaceColor','m');a.FaceAlpha = 0.3;a.LineStyle='None';a.BaseLine.LineStyle='None';set(gca,'visible','off');hold on;
plot(mean([Int_data_temp_mean(1013+5:1115)-min(Int_data_temp_mean(1013+5:1115));Int_data_temp_mean(1473+5:1575)-min(Int_data_temp_mean(1473+5:1575))],1),'k','LineWidth',4);axis([0 220 -0.5 3.5]);set(gca,'visible','off');
axes(ha(6));
a=area(Speed_flow(1,1825-delay:2005),'FaceColor','m');a.FaceAlpha = 0.3;a.LineStyle='None';a.BaseLine.LineStyle='None';set(gca,'visible','off');hold on;
plot(Int_data_temp_mean(1825+5:2005)-min(Int_data_temp_mean(1825+5:2005)),'k','LineWidth',4);axis([0 220 -0.5 3.5]);set(gca,'visible','off');
axes(ha(7));
a=area(Speed_flow(1,2080-delay:2245),'FaceColor','m');a.FaceAlpha = 0.3;a.LineStyle='None';a.BaseLine.LineStyle='None';set(gca,'visible','off');hold on;
plot(Int_data_temp_mean(2080+5:2245)-min(Int_data_temp_mean(2080+5:2245)),'k','LineWidth',4);axis([0 220 -0.5 3.5]);set(gca,'visible','off');

axes(ha(8));i=1;bar(nanmean(int_quant_perFish(:,i),1),'FaceColor',[0.5 0.5 0.5]);axis([0.5 2 0 275]);hold on;er = errorbar(nanmean(int_quant_perFish(:,i),1),nanstd(int_quant_perFish(:,i),1,1));er.Color = [0 0 0];er.LineWidth = 2;hold off;set(gca,'box','off');set(gca,'color','none');set(gca,'xtick',[])
axes(ha(9));i=2;bar(nanmean(int_quant_perFish(:,i),1),'FaceColor',[0.5 0.5 0.5]);axis([0.5 2 0 275]);hold on;er = errorbar(nanmean(int_quant_perFish(:,i),1),nanstd(int_quant_perFish(:,i),1,1));er.Color = [0 0 0];er.LineWidth = 2;hold off;set(gca,'visible','off');
axes(ha(10));i=3;bar(nanmean(int_quant_perFish(:,i),1),'FaceColor',[0.5 0.5 0.5]);axis([0.5 2 0 275]);hold on;er = errorbar(nanmean(int_quant_perFish(:,i),1),nanstd(int_quant_perFish(:,i),1,1));er.Color = [0 0 0];er.LineWidth = 2;hold off;set(gca,'visible','off');
axes(ha(11));i=4;bar(nanmean(int_quant_perFish(:,i),1),'FaceColor',[0.5 0.5 0.5]);axis([0.5 2 0 275]);hold on;er = errorbar(nanmean(int_quant_perFish(:,i),1),nanstd(int_quant_perFish(:,i),1,1));er.Color = [0 0 0];er.LineWidth = 2;hold off;set(gca,'visible','off');
axes(ha(12));i=5;bar(nanmean(int_quant_perFish(:,i),1),'FaceColor',[0.5 0.5 0.5]);axis([0.5 2 0 275]);hold on;er = errorbar(nanmean(int_quant_perFish(:,i),1),nanstd(int_quant_perFish(:,i),1,1));er.Color = [0 0 0];er.LineWidth = 2;hold off;set(gca,'visible','off');
axes(ha(13));i=6;bar(nanmean(int_quant_perFish(:,i),1),'FaceColor',[0.5 0.5 0.5]);axis([0.5 2 0 275]);hold on;er = errorbar(nanmean(int_quant_perFish(:,i),1),nanstd(int_quant_perFish(:,i),1,1));er.Color = [0 0 0];er.LineWidth = 2;hold off;set(gca,'visible','off');
axes(ha(14));i=7;bar(nanmean(int_quant_perFish(:,i),1),'FaceColor',[0.5 0.5 0.5]);axis([0.5 2 0 275]);hold on;er = errorbar(nanmean(int_quant_perFish(:,i),1),nanstd(int_quant_perFish(:,i),1,1));er.Color = [0 0 0];er.LineWidth = 2;hold off;set(gca,'visible','off');

print(Fighandle,strcat('D:\Dropbox\Papers\Flow paper\Revision\Test_integrator.png'),'-dpng','-r0');


PrismTemp=int_quant_perFish(:,[1 4 5 7 6 2 3]);