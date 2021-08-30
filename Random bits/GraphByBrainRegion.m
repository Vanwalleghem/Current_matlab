RegionList_select={'Thalamus','Cerebellum','Semicircularis','Telencephalon','Tectum','Tegmentum','Habenula','Pretectum','MON','Hindbrain','pLLG'};

Node_selectRegion=zeros(1,length(Node_ID));
for i=1:length(Node_ID)
    Min_dist=[inf,0];
    for j=1:length(RegionList_select)
        regionName=RegionList_select{j};        
        if strcmp(regionName,'Telencephalon')
            Mask=Zbrain_Masks{294,3};
        elseif strcmp(regionName,'Hindbrain')
            Hindbrain_Mask=Zbrain_Masks{259,3};
            Mask=Zbrain_Masks{131,3};
            IsInEyes_temp=ismember(Hindbrain_Mask,Mask,'rows');IsInEyes_temp=find(IsInEyes_temp==1);%remove cerebellum
            Hindbrain_Mask(IsInEyes_temp,:)=[];
            Mask=Zbrain_Masks{295,3};
            IsInEyes_temp=ismember(Hindbrain_Mask,Mask,'rows');IsInEyes_temp=find(IsInEyes_temp==1);%remove MON
            Hindbrain_Mask(IsInEyes_temp,:)=[];
            Mask=Hindbrain_Mask;
            clearvars Hindbrain_Mask IsInEyes_temp;
        elseif strcmp(regionName,'pLLG')
            Mask=Zbrain_Masks{90,3};
            Mask=vertcat(Mask,Zbrain_Masks{93,3});
        elseif strcmp(regionName,'pLLG')
            Mask=Zbrain_Masks{90,3};
            Mask=vertcat(Mask,Zbrain_Masks{93,3});
        else
            Mask=[];
            IndexC=strfind({Zbrain_Masks{:,2}}, regionName);
            IndexC=find(not(cellfun('isempty', IndexC)));
            for h=IndexC
                if isempty(Mask)
                    Mask=Zbrain_Masks{h,3};
                else
                    Mask=vertcat(Mask,Zbrain_Masks{h,3});
                end
            end
        end
        Mask=unique(Mask,'rows');
        dist_temp=min(pdist2(round(Node_all_temp(i,:)),Mask));
        if dist_temp<Min_dist(1)
            Min_dist(1)=dist_temp;
            Min_dist(2)=j;
        end
    end
    Node_selectRegion(i)=Min_dist(2);
end

spectral=cbrewer('div','Spectral',max(Node_selectRegion));
spectral=cbrewer('qual','Paired',max(Node_selectRegion));
for i=1:length(Node_selectRegion)
    Node_colors(i,:)=spectral(Node_selectRegion(i),:);
end

Fighandle=figure;
graph_temp=graph(Correlation_bin_prop);
set(Fighandle, 'Position', [10,10, 600, 1400]);
plot(Zbrain_brainMask2D(:,2),1400-Zbrain_brainMask2D(:,1),'k');hold on;
scatter(Node_all(:,2),Node_all(:,1),(degree(graph_temp)/5)+30,'k','filled'); %colormap jet;
scatter(Node_all(:,2),Node_all(:,1),(degree(graph_temp)/5)+20,Node_colors,'filled');axis([0 600 0 1400]);camroll(180); %colormap jet;
set(gca,'Visible','off')
print(Fighandle,strcat('D:\Pictures\processed\Flow\Basic\Figure\Graph Theory\','PerBrainRegion.tif'),'-dtiff','-r0');

P_brain=struct();
P_brain.coef=participation_coef(Correlation_bin_prop,Node_selectRegion);
Fighandle=figure;
set(Fighandle, 'Position', [10,10, 600, 1400]);
plot(Zbrain_brainMask2D(:,2),1400-Zbrain_brainMask2D(:,1),'k');hold on;
scatter(Node_all(:,2),Node_all(:,1),90,'k','filled');hold on;
scatter(Node_all(:,2),Node_all(:,1),55,P_brain.coef,'filled');colormap hot;hold off;axis([0 600 0 1400]);camroll(180);caxis([0 1]);
set(gca,'Visible','off')
print(Fighandle,strcat('D:\Pictures\processed\Flow\Basic\Figure\Graph Theory\','Participation_perBrain.tif'),'-dtiff','-r0');

edges=[0:0.1:1];
h=nan(max(Node_selectRegion),length(edges)-1,13);
h_mean=[];
for i=1:max(Node_selectRegion)    
    for fish_nb=1:13
        Corr_temp=squeeze(CorrelationMatrices_minusFish(:,:,fish_nb));
        Corr_temp=1-Corr_temp;
        Corr_temp=weight_conversion(Corr_temp,'autofix');
        Corr_temp=weight_conversion(Corr_temp,'normalize');
        P_temp=participation_coef(threshold_proportional(Corr_temp,0.25),Node_selectRegion);
        h(i,:,fish_nb)=histcounts(P_temp(Node_selectRegion==i),edges,'normalization','probability');
        h_mean(i,fish_nb)=mean(P_temp(Node_selectRegion==i));
    end
end
clearvars temp i j k h h_mean s S edges Q Q0 Q1

edges=[0:0.1:1];
h=nan(max(Node_selectRegion),length(edges)-1,13);
h_mean=[];
for i=1:max(Node_selectRegion)    
    thr_nb=1;
    for Threshold=[0.95:-0.05:0.5]        
        P_temp=participation_coef(threshold_absolute(CorrelationMatrix,Threshold),Node_selectRegion);
        h(i,:,thr_nb)=histcounts(P_temp(Node_selectRegion==i),edges,'normalization','probability');
        h_mean(i,thr_nb)=mean(P_temp(Node_selectRegion==i));
        thr_nb=thr_nb+1;
    end
end


edges=[0:0.1:1];
h=nan(max(Node_selectRegion),length(edges)-1,length([0.95:-0.05:0.5]),13);
h_mean=[];PrismTemp2=[];
for i=1:max(Node_selectRegion)
    thr_nb=1;
    PrismTemp=[];    
    for Threshold=[0.95:-0.05:0.5]
        for fish_nb=1:13        
            Corr_temp=squeeze(CorrelationMatrices_minusFish(:,:,fish_nb));
            Corr_temp=1-Corr_temp;
            Corr_temp=weight_conversion(Corr_temp,'autofix');
            Corr_temp=weight_conversion(Corr_temp,'normalize');
            P_temp=participation_coef(threshold_absolute(Corr_temp,Threshold),Node_selectRegion);
            h(i,:,thr_nb,fish_nb)=histcounts(P_temp(Node_selectRegion==i),edges,'normalization','probability');
            h_mean(i,thr_nb,fish_nb)=mean(P_temp(Node_selectRegion==i));
            PrismTemp=horzcat(PrismTemp,mean(P_temp(Node_selectRegion==i)));
        end
        thr_nb=thr_nb+1;
    end
    PrismTemp2=vertcat(PrismTemp2,PrismTemp);
end