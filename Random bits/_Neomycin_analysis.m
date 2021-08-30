load('CNMF_analysis_flow_Neomycin.mat')

flow=zeros(6,655);
GCaMP6=[-0.104392135015146,1.69644104899772,5.13796058542217,8.27886020152244,10.3756715204800,11.8173714529814,12.2425184714093,10.8571417354877,8.80831829681196,6.91339112244670,5.46959264663869,4.30868766622567,3.42533619066766,2.75378443486879,2.18017250852183,1.72816235135824,1.32732537295463,1.00684435500268,0.730210038304555,0.530242444093118,0.362253250339685,0.227668255288566,0.0869242416152502,0.000718266708050853,-0.0828334873368325]';
%GCaMP6=[0.000256990000000000;0.00850739000000000;0.0654158300000000;0.0784609000000000;0.0764130100000000;0.0665958600000000;0.0579028900000000;0.0467942900000000;0.0232079800000000;0.0144564400000000;0.00695772000000000;0.00526551000000000;0.00299500000000000;0.00198520000000000;0.00128512000000000;0.00134175000000000;0.000403170000000000;0];
back=[57 257 457];
back_off=[106 306 506];
fwd=[157 357 557];
fwd_off=[207 407 607];
flow(1,back(1):back(1)+size(GCaMP6,1)-1)=GCaMP6';
flow(1,back(2):back(2)+size(GCaMP6,1)-1)=GCaMP6';
flow(1,back(3):back(3)+size(GCaMP6,1)-1)=GCaMP6';
flow(2,back_off(1):back_off(1)+size(GCaMP6,1)-1)=GCaMP6';
flow(2,back_off(2):back_off(2)+size(GCaMP6,1)-1)=GCaMP6';
flow(2,back_off(3):back_off(3)+size(GCaMP6,1)-1)=GCaMP6';
flow(3,fwd(1):fwd(1)+size(GCaMP6,1)-1)=GCaMP6';
flow(3,fwd(2):fwd(2)+size(GCaMP6,1)-1)=GCaMP6';
flow(3,fwd(3):fwd(3)+size(GCaMP6,1)-1)=GCaMP6';
flow(4,fwd_off(1):fwd_off(1)+size(GCaMP6,1)-1)=GCaMP6';
flow(4,fwd_off(2):fwd_off(2)+size(GCaMP6,1)-1)=GCaMP6';
flow(4,fwd_off(3):fwd_off(3)+size(GCaMP6,1)-1)=GCaMP6';
flow(5,back(1):back(1)+43)=1;
flow(5,back(2):back(2)+43)=1;
flow(5,back(3):back(3)+43)=1;
flow(6,fwd(1):fwd(1)+43)=1;
flow(6,fwd(2):fwd(2)+43)=1;
flow(6,fwd(3):fwd(3)+43)=1;
clearvars GCaMP6 back back_off fwd fwd_off;

parfor i=1:size(ZS2,1)
    mdl=fitlm(flow',ZS2(i,:));    
    model_flow(i).coef=mdl.Coefficients;        
    model_flow(i).rsquared=mdl.Rsquared.Adjusted;
end
clearvars slow_fast slow_fast_fwd i

idx_rsq_flow=find([model_flow.rsquared]>0.25);
[idxKmeans_rsq_flow Cmap_rsq_flow]=kmeans(ZS2(idx_rsq_flow,:),5,'Options',options,'Distance','cityblock','Replicates',3,'MaxIter',1000,'Display','final');
[Model_rsq_f,GoodBetas_rsq_f]=Test_Regress(Cmap_rsq_flow,flow,idxKmeans_rsq_flow,0.4);


parfor i=1:size(ZS2,1)
    mdl=fitlm(Basic_Clusters',ZS2(i,:));    
    model_basic(i).coef=mdl.Coefficients;        
    model_basic(i).rsquared=mdl.Rsquared.Adjusted;
end
clearvars slow_fast slow_fast_fwd i


idx_rsq=find([model_DF_Thr5.rsquared]>0.25);
[idxKmeans_rsq Cmap_rsq]=kmeans(ZS2(idx_rsq,:),3,'Options',options,'Distance','cityblock','Replicates',3,'MaxIter',1000,'Display','final');
[Model_rsq,GoodBetas_rsq]=Test_Regress(Cmap_rsq,flow,idxKmeans_rsq,0.4);

Fighandle=figure;
set(Fighandle, 'Position', [100, 100, 800, 800]);
imagesc(ZS2(idx_rsq(randperm(length(idx_rsq))),:),[-1 4]);colormap hot
print(Fighandle,strcat('D:\Pictures\processed\Flow\Basic\Neomycin\_ZS_rsq.svg'),'-dsvg','-r0');

Fighandle=figure;
set(Fighandle, 'Position', [100, 100, 1200, 600]);
histogram([model_DF_Thr5.rsquared],'Normalization','probability');hold on;histogram(rsq_Basic,'Normalization','probability');
print(Fighandle,strcat('D:\Pictures\processed\Flow\Basic\Neomycin\_rsq_histogram.svg'),'-dsvg','-r0');

colors_b=colors;

%Make figure of full length + STD + rasterhisto
Fighandle=figure;
set(Fighandle, 'Position', [100, 100, 800, 800]);
set(findall(Fighandle,'type','text'),'fontSize',12,'fontWeight','bold','FontName','Arial')
counter=1;counter2=1;xplot=2;yplot=2;
StimLength=655;
x = linspace(0.2,StimLength/5,StimLength);
ha = tight_subplot(xplot,yplot,[.01 .01],[.01 .01],[.01 .01]);
colors=colors_b([4 7],:);
for i=GoodBetas_rsq    
    idx_temp=find(idxKmeans_rsq==i);
    temp=mean(ZS2(idx_rsq(idx_temp),:),1);
    std_temp=std(ZS2(idx_rsq(idx_temp),:),1,1);
    %subplot(xplot,yplot,(2*counter)-1);
    axes(ha((2*counter)-1));
    H=shadedErrorBar(x, temp, std_temp);axis([0 130 -1 4]);
    H.mainLine.Color=colors(counter,:)/256;
    H.patch.FaceColor=colors(counter,:)/256;
    H.edge(1).Color=colors(counter,:)/256;
    H.edge(2).Color=colors(counter,:)/256;
    H.mainLine.LineWidth=3;    
    set(gca,'XTickLabel',[]);set(gca,'YTickLabel',[]);  
    %subplot(xplot,yplot,(2*counter));
    axes(ha((2*counter)));
    imagesc(ZS2(idx_rsq(idx_temp(randperm(length(idx_temp)))),:),[-0.5 5]);colormap hot
    set(gca,'XTickLabel',[]);set(gca,'YTickLabel',[]);    
    counter=counter+1;
end
print(Fighandle,strcat('D:\Pictures\processed\Flow\Basic\Figure\FullLength_neo_clusters.svg'),'-dsvg','-r0');    

GoodBetas_select=GoodBetas_rsq;
idxKmeans=zeros(size(ZS2,1),1);
idxKmeans(idx_rsq)=idxKmeans_rsq;
GoodClustersData=[];
for i=1:length(GoodBetas_select)
    GoodClustersData(i).ZS=ZS2(idxKmeans==GoodBetas_select(i),:);
    GoodClustersData(i).Mean=mean(GoodClustersData(i).ZS,1);
    GoodClustersData(i).STD=std(GoodClustersData(i).ZS,1,1);
end

for i=1:numel(GoodClustersData)
    corr_temp=zeros(size(GoodClustersData(i).ZS,1),1);
    parfor j=1:size(GoodClustersData(i).ZS,1)
        temp=corrcoef(GoodClustersData(i).Mean, GoodClustersData(i).ZS(j,:));
        corr_temp(j)=temp(1,2);
    end
    GoodClustersData(i).CorrCoef=corr_temp;
end


GoodClusters_goodmembers=[];Threshold=0.5;
idxKmeans_ZS_goodmembers=zeros(1,size(ZS2,1));
for i=1:length(GoodBetas_select)
%GoodClusters_goodmembers(i).Spikes=GoodClustersData(i).Spikes(find(GoodClustersData(i).CorrCoef>=0.5),:);
%GoodClusters_goodmembers(i).ZS=zscore(GoodClustersData(i).DF(find(GoodClustersData(i).CorrCoef>=0.5),:),1,2);
GoodClusters_goodmembers(i).ZS=GoodClustersData(i).ZS(find(GoodClustersData(i).CorrCoef>=Threshold),:);
temp=find(idxKmeans==GoodBetas_select(i));
GoodClusters_goodmembers(i).idx=temp(find(GoodClustersData(i).CorrCoef>=0.5));
GoodClusters_goodmembers(i).mean=mean(GoodClusters_goodmembers(i).ZS,1);
GoodClusters_goodmembers(i).STD=std(GoodClusters_goodmembers(i).ZS,1,1);
idx=find(idxKmeans==GoodBetas_select(i));
idx=idx(find(GoodClustersData(i).CorrCoef>=Threshold));
idxKmeans_ZS_goodmembers(idx)=GoodBetas_select(i);
%GoodClusters_goodmembers(i).Fish=idx_Fish(idx);
end

colors_backup=colors;


colors=colors([4 7],:);

%Make figure of full length + STD + rasterhisto
Fighandle=figure;
set(Fighandle, 'Position', [100, 100, 1200, 600]);
set(findall(Fighandle,'type','text'),'fontSize',12,'fontWeight','bold','FontName','Arial')
counter=1;counter2=1;xplot=2;yplot=2;
StimLength=655;
x = linspace(0.2,StimLength/5,StimLength);
ha = tight_subplot(xplot,yplot,[.01 .01],[.01 .01],[.01 .01]);
for i=GoodBetas_select    
    temp=GoodClusters_goodmembers(counter).mean;
    std_temp=GoodClusters_goodmembers(counter).STD;        
    %subplot(xplot,yplot,(2*counter)-1);
    axes(ha((2*counter)-1));
    H=shadedErrorBar(x, temp, std_temp);axis([0 130 -2 5]);
    H.mainLine.Color=colors(counter,:)/256;
    H.patch.FaceColor=colors(counter,:)/256;
    H.edge(1).Color=colors(counter,:)/256;
    H.edge(2).Color=colors(counter,:)/256;
    H.mainLine.LineWidth=3;    hold on;
    rectangle('FaceColor','m','Position',[11 -2 10 0.25]);rectangle('FaceColor','m','Position',[51 -2 10 0.25]);rectangle('FaceColor','m','Position',[91 -2 10 0.25]);rectangle('FaceColor','g','Position',[31 -2 10 0.25]);rectangle('FaceColor','g','Position',[71 -2 10 0.25]);rectangle('FaceColor','g','Position',[111 -2 10 0.25]);    
    set(gca,'XTickLabel',[]);set(gca,'YTickLabel',[]);  
    %subplot(xplot,yplot,(2*counter));
    axes(ha((2*counter)));
    imagesc(GoodClusters_goodmembers(counter).ZS(randperm(size(GoodClusters_goodmembers(counter).ZS,1)),:),[-0.5 5]);colormap hot
    set(gca,'XTickLabel',[]);set(gca,'YTickLabel',[]);    
    counter=counter+1;
end
print(Fighandle,'D:\Pictures\processed\Flow\Basic\Figure\FullLength_Neo_clusters.svg','-dsvg','-r0');    




back=[57 257 457];
fwd=[157 357 557];
PrismBackTemp=nan(length(GoodBetas_select),length(unique(idx_Fish)));counter=1;
PrismFWDTemp=nan(length(GoodBetas_select),length(unique(idx_Fish)));counter=1;
for i=GoodBetas_select
    idx_temp=find(idxKmeans_ZS_goodmembers==i);
    for fish_nb=1:3
        idx_fish_temp=idx_Fish(idx_temp)==fish_nb;
        ZS_temp=mean(ZS2(idx_temp(idx_fish_temp),:),1);
        sum(idx_fish_temp)
        max_resp=zeros(2,3);        
        for rep=2:3
        max_resp(1,rep)=max(ZS_temp(:,back(rep):back(rep)+20),[],2)-min(ZS_temp(:,back(rep)-5:back(rep)+20),[],2);
        max_resp(2,rep)=max(ZS_temp(:,fwd(rep):fwd(rep)+20),[],2)-min(ZS_temp(:,fwd(rep)-5:fwd(rep)+20),[],2);
        end        
        PrismBackTemp(counter,fish_nb)=round(mean(max_resp(1,:)),3);
        PrismFWDTemp(counter,fish_nb)=round(mean(max_resp(2,:)),3);
    end    
    counter=counter+1;
end

back=[57 257 457];
fwd=[157 357 557];
FishList=unique(basic_idx_Fish.idx_Fish);
PrismBackTemp=nan(length(GoodBetas_select),length(unique(basic_idx_Fish.idx_Fish)));counter=1;
PrismFWDTemp=nan(length(GoodBetas_select),length(unique(basic_idx_Fish.idx_Fish)));counter=1;
for i=[4 7]
    idx_temp=basic_data.GoodClusters_goodmembers(i).idx;
    for fish_nb=1:length(FishList)
        idx_fish_temp=basic_idx_Fish.idx_Fish(idx_temp)==FishList(fish_nb);
        ZS_temp=basic_data.GoodClusters_goodmembers(i).ZS;
        ZS_temp=mean(ZS_temp(idx_fish_temp,:),1);
        max_resp=zeros(2,3);        
        for rep=2:3
        max_resp(1,rep)=max(ZS_temp(:,back(rep):back(rep)+20),[],2)-min(ZS_temp(:,back(rep)-5:back(rep)+20),[],2);
        max_resp(2,rep)=max(ZS_temp(:,fwd(rep):fwd(rep)+20),[],2)-min(ZS_temp(:,fwd(rep)-5:fwd(rep)+20),[],2);
        end        
        PrismBackTemp(counter,fish_nb)=round(mean(max_resp(1,:)),3);
        PrismFWDTemp(counter,fish_nb)=round(mean(max_resp(2,:)),3);
    end    
    counter=counter+1;
end

%histogram of strength of response
counter=1;
mean_resp_strength={};
for i=GoodBetas_select
    idx_temp=find(idxKmeans_ZS_goodmembers==i);
    ZS_temp=ZS2(idx_temp,:);
    max_resp=zeros(size(ZS_temp,1),3,2);
    for rep=2:3
        max_resp(:,rep,1)=max(ZS_temp(:,back(rep):back(rep)+20),[],2)-min(ZS_temp(:,back(rep)-5:back(rep)+20),[],2);
        max_resp(:,rep,2)=max(ZS_temp(:,fwd(rep):fwd(rep)+20),[],2)-min(ZS_temp(:,fwd(rep)-5:fwd(rep)+20),[],2);
    end
    mean_resp_strength{counter,1}=squeeze(mean(max_resp(:,:,1),2));
    mean_resp_strength{counter,2}=squeeze(mean(max_resp(:,:,2),2));
    counter=counter+1;
end

counter=1;
mean_resp_strength_basic={};
for i=[4 7]    
    ZS_temp=basic_data.GoodClusters_goodmembers(i).ZS;
    max_resp=zeros(size(ZS_temp,1),3,2);
    for rep=2:3
        max_resp(:,rep,1)=max(ZS_temp(:,back(rep):back(rep)+20),[],2)-min(ZS_temp(:,back(rep)-5:back(rep)+20),[],2);
        max_resp(:,rep,2)=max(ZS_temp(:,fwd(rep):fwd(rep)+20),[],2)-min(ZS_temp(:,fwd(rep)-5:fwd(rep)+20),[],2);
    end
    mean_resp_strength_basic{counter,1}=squeeze(mean(max_resp(:,:,1),2));
    mean_resp_strength_basic{counter,2}=squeeze(mean(max_resp(:,:,2),2));
    counter=counter+1;
end

%Make figure of full length + STD + rasterhisto
edges=[0:0.4:6];
Fighandle=figure;
set(Fighandle, 'Position', [100, 100, 1200, 600]);
set(findall(Fighandle,'type','text'),'fontSize',12,'fontWeight','bold','FontName','Arial')
counter=1;counter2=1;xplot=2;yplot=2;
StimLength=655;
x = linspace(0.2,StimLength/5,StimLength);
ha = tight_subplot(xplot,yplot,[.01 .01],[.01 .01],[.01 .01]);
axes(ha(1));
histogram(mean_resp_strength{1,1},edges,'Normalization','probability');hold on;histogram(mean_resp_strength_basic{1,1},edges,'Normalization','probability');ylim([0 0.4]);set(gca,'XTickLabel',[]);set(gca,'YTickLabel',[]);
axes(ha(2));histogram(mean_resp_strength{1,2},edges,'Normalization','probability');hold on;histogram(mean_resp_strength_basic{1,2},edges,'Normalization','probability');ylim([0 0.4]);set(gca,'XTickLabel',[]);set(gca,'YTickLabel',[]);
axes(ha(3));histogram(mean_resp_strength{2,1},edges,'Normalization','probability');hold on;histogram(mean_resp_strength_basic{2,1},edges,'Normalization','probability');ylim([0 0.4]);set(gca,'XTickLabel',[]);set(gca,'YTickLabel',[]);
axes(ha(4));histogram(mean_resp_strength{2,2},edges,'Normalization','probability');hold on;histogram(mean_resp_strength_basic{2,2},edges,'Normalization','probability');ylim([0 0.4]);set(gca,'XTickLabel',[]);set(gca,'YTickLabel',[]);
print(Fighandle,'D:\Pictures\processed\Flow\Basic\Neomycin\StrengthHistogram.svg','-dsvg','-r0');    


% %--------------------------------
% load('__BasicFlow_final_lite.mat')
% Neo=load('__BasicFlow_Neo.mat')
% 
% %Make figure of full length + STD + rasterhisto
% Fighandle=figure;
% set(Fighandle, 'Position', [100, 100, 1200, 1200]);
% set(findall(Fighandle,'type','text'),'fontSize',12,'fontWeight','bold','FontName','Arial')
% counter=1;counter2=1;xplot=4;yplot=2;
% StimLength=655;
% x = linspace(0.2,StimLength/5,StimLength);
% ha = tight_subplot(xplot,yplot,[.01 .01],[.01 .01],[.01 .01]);
% for i=[5 7]   
%     temp=GoodClusters_goodmembers(i).mean;
%     std_temp=GoodClusters_goodmembers(i).STD;        
%     %subplot(xplot,yplot,(2*counter)-1);
%     axes(ha((2*counter)-1));
%     H=shadedErrorBar(x, temp, std_temp);axis([0 130 -1 4]);
%     H.mainLine.Color=colors(i,:)/256;
%     H.patch.FaceColor=colors(i,:)/256;
%     H.edge(1).Color=colors(i,:)/256;
%     H.edge(2).Color=colors(i,:)/256;
%     H.mainLine.LineWidth=3;    
%     set(gca,'XTickLabel',[]);set(gca,'YTickLabel',[]);%title(size(GoodClusters_goodmembers(i).ZS,1));
%     %subplot(xplot,yplot,(2*counter));
%     axes(ha((2*counter)));
%     imagesc(GoodClusters_goodmembers(i).ZS(randperm(size(GoodClusters_goodmembers(i).ZS,1)),:),[-1 3]);colormap hot
%     set(gca,'XTickLabel',[]);set(gca,'YTickLabel',[]);    
%     counter=counter+1;
%     temp=Neo.GoodClusters_goodmembers(counter2).mean;
%     std_temp=Neo.GoodClusters_goodmembers(counter2).STD;        
%     %subplot(xplot,yplot,(2*counter)-1);
%     axes(ha((2*counter)-1));
%     H=shadedErrorBar(x, temp, std_temp);axis([0 130 -1 4]);
%     H.mainLine.Color=colors(counter,:)/256;
%     H.patch.FaceColor=colors(counter,:)/256;
%     H.edge(1).Color=colors(counter,:)/256;
%     H.edge(2).Color=colors(counter,:)/256;
%     H.mainLine.LineWidth=3;    %title(size(Neo.GoodClusters_goodmembers(counter2).ZS,1));
%     set(gca,'XTickLabel',[]);set(gca,'YTickLabel',[]);  
%     %subplot(xplot,yplot,(2*counter));
%     axes(ha((2*counter)));
%     imagesc(Neo.GoodClusters_goodmembers(counter2).ZS(randperm(size(Neo.GoodClusters_goodmembers(counter2).ZS,1)),:),[-1 3]);colormap hot
%     set(gca,'XTickLabel',[]);set(gca,'YTickLabel',[]);    
%     counter=counter+1;
%     counter2=counter2+1;
% end
% print(Fighandle,strcat('D:\Pictures\processed\Flow\Basic\Figure\FullLength_neoAndbasic_clusters.svg'),'-dsvg','-r0');    
% 
% 
% 
% %Make figure of full length + STD + rasterhisto
% counter=1;counter2=1;
% for i=[5 7]   
%     size(GoodClusters_goodmembers(i).ZS,1)/length(unique(GoodClusters_goodmembers(i).Fish))
%     idx_temp=Neo.GoodClusters_goodmembers(counter2).idx;
%     size(Neo.GoodClusters_goodmembers(counter2).ZS,1)/length(unique(Neo.idx_Fish(idx_temp)))
%     counter2=counter2+1;
% end

