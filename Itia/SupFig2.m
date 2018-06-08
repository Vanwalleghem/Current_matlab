%Graph figure with thick lines and no SD
counter=1;xplot=3;yplot=1;
x = linspace(0.25,246/4,246);
for region_nb=9
    Fighandle=figure;
    set(Fighandle, 'Position', [100, 100, 600, 250]);
    counter=1;
    regionName=ItiaList{region_nb};
    idx_temp=LinReg.(regionName).KmeansIdx_select_AVG;
    GoodBet_temp=LinReg.(regionName).GoodBetas_AVG_final;
    if region_nb==5
        counter=counter+1;
    end
    for j=1:length(GoodBet_temp)
        idx_temp2=find(idx_temp==GoodBet_temp(j));
        for k=0:5
            meanToPlot=mean(LinReg.(regionName).ZS_AVG(idx_temp2,4+(k*41):40+(k*41)),1);
            hold on;
            plot(x(4+(k*40):40+(k*40)),meanToPlot-mean(meanToPlot(1:4)),'color',colors{region_nb}(j,:),'LineWidth',3);axis([0 60 -3 3]);            
            if (j==1 & k<5)
                rectangle('EdgeColor','none','FaceColor',[0.25, 0.25, 0.25, 0.4-k*0.06],'Position',[x(9+(k*40)) -3 1 7]);
            end
        end
        hold on;
    end
    hold off;
    
    print(Fighandle,regionName,'-dsvg','-r0');
    close all;
end


figure;
counter=1;
x = linspace(0.25,760/4,760);
order=[5 9 1 12 6 8 3 2 10 11];xplot=length(order);yplot=3;counter2=1;counter3=1;
ha = tight_subplot(xplot,yplot,[.01 .01],[.01 .01],[.01 .01]);
MeanClusPerBrain={};
STDClusPerBrain={};
for i=order
    regionName=ItiaList{i};
    if i>10
        idx_temp=LinReg.(regionName).KmeansIdx_merge;
        GoodBet_temp=LinReg.(regionName).GoodBeta_merge;
    else
        idx_temp=LinReg.(regionName).KmeansIdx_select_AVG;
        GoodBet_temp=LinReg.(regionName).GoodBetas_AVG_final;
    end
    counter=counter2;
    if i==5
        counter=counter+1;
    end
    for j=1:length(GoodBet_temp)
        idx_temp2=find(idx_temp==GoodBet_temp(j));
        %subplot(xplot,yplot,counter+(j-1));
        axes(ha(counter+(j-1)));
        meanToPlot=mean(LinReg.(regionName).ZS_rsq(idx_temp2,:),1);
        MeanClusPerBrain{counter3,j}=meanToPlot;
        std_95=std(LinReg.(regionName).ZS_rsq(idx_temp2,:),1,1);%confidence_intervals(LinReg.(regionName).ZS_rsq(idx_temp2,:),95);
        %std(LinReg.(regionName).ZS_rsq(idx_temp2,:),1,1);
        STDClusPerBrain{counter3,j}=std_95;
        H=shadedErrorBar(x(:),meanToPlot-mean(meanToPlot(1:4)), std_95);axis([0 190 -5 5]);
        H.mainLine.Color=colors{i}(j,:);
        H.patch.FaceColor=colors{i}(j,:)/2;
        H.edge(1).Color=colors{i}(j,:)/2;
        H.edge(2).Color=colors{i}(j,:)/2;
        set(gca,'XTickLabel',[]);set(gca,'YTickLabel',[]);
    end
    counter2=counter2+yplot;
    counter3=counter3+1;
end
clearvars i j k counter counter2 H idx_temp GoodBet_temp std_95 meanToPlot


idxStart=40;
Max_response=zeros(10,18,3);
Max_response_STD=zeros(10,18,3);
Min_response=zeros(10,18,3);
Min_response_STD=zeros(10,18,3);
for j=1:size(MeanClusPerBrain,1)
    for k=1:3
        if ~isempty(MeanClusPerBrain{j,k}) 
            for i=1:18
                [Max_response(j,i,k), idx]=max(MeanClusPerBrain{j,k}((idxStart+(i-1)*40):(idxStart+(i-1)*40)+10));
                Max_response_STD(j,i,k)=STDClusPerBrain{j,k}(idxStart+(i-1)*40)+idx-1;
                [Min_response(j,i,k), idx]=min(MeanClusPerBrain{j,k}((idxStart+(i-1)*40):(idxStart+(i-1)*40)+10));
                Min_response_STD(j,i,k)=STDClusPerBrain{j,k}(idxStart+(i-1)*40)+idx-1;
            end
        end
    end
end


PrismTemp=zeros(9,10);
for i=1:10
    for j=1:3
        PrismTemp((j*3)-2,i)=Max_response(i;
        PrismTemp((j*3)-1,i)=Max_response_STD;
        PrismTemp((j*3),i)=13;
    end
end