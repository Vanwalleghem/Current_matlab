load('CNMF_analysis_flow_Neomycin.mat')

parfor i=1:size(ZS2,1)
    mdl=fitlm(Basic_Clusters',ZS2(i,:));    
    model_basic(i).coef=mdl.Coefficients;        
    model_basic(i).rsquared=mdl.Rsquared.Adjusted;
end
clearvars slow_fast slow_fast_fwd i


idx_rsq=find([model_DF_Thr5.rsquared]>0.2);
[idxKmeans_rsq Cmap_rsq]=kmeans(ZS2(idx_rsq,:),5,'Options',options,'Distance','cityblock','Replicates',3,'MaxIter',1000,'Display','final');
[Model_rsq,GoodBetas_rsq]=Test_Regress(Cmap_rsq,flow,idxKmeans_rsq,0.4);

Fighandle=figure;
set(Fighandle, 'Position', [100, 100, 800, 800]);
imagesc(ZS2(idx_rsq(randperm(length(idx_rsq))),:),[-1 4]);colormap hot
print(Fighandle,strcat('D:\Pictures\processed\Flow\Basic\Neomycin\_ZS_rsq.svg'),'-dsvg','-r0');

Fighandle=figure;
set(Fighandle, 'Position', [100, 100, 800, 800]);
histogram([model_DF_Thr5.rsquared],'Normalization','probability');hold on;histogram(rsq_Basic,'Normalization','probability');
print(Fighandle,strcat('D:\Pictures\processed\Flow\Basic\Neomycin\_rsq_histogram.svg'),'-dsvg','-r0');


%Make figure of full length + STD + rasterhisto
Fighandle=figure;
set(Fighandle, 'Position', [100, 100, 800, 800]);
set(findall(Fighandle,'type','text'),'fontSize',12,'fontWeight','bold','FontName','Arial')
counter=1;counter2=1;xplot=2;yplot=2;
StimLength=655;
x = linspace(0.2,StimLength/5,StimLength);
ha = tight_subplot(xplot,yplot,[.01 .01],[.01 .01],[.01 .01]);
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
    imagesc(ZS2(idx_rsq(idx_temp(randperm(length(idx_temp)))),:),[-1 3]);colormap hot
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

%Make figure of full length + STD + rasterhisto
Fighandle=figure;
set(Fighandle, 'Position', [100, 100, 800, 800]);
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
    H=shadedErrorBar(x, temp, std_temp);axis([0 130 -1 4]);
    H.mainLine.Color=colors(counter,:)/256;
    H.patch.FaceColor=colors(counter,:)/256;
    H.edge(1).Color=colors(counter,:)/256;
    H.edge(2).Color=colors(counter,:)/256;
    H.mainLine.LineWidth=3;    
    set(gca,'XTickLabel',[]);set(gca,'YTickLabel',[]);  
    %subplot(xplot,yplot,(2*counter));
    axes(ha((2*counter)));
    imagesc(GoodClusters_goodmembers(counter).ZS(randperm(size(GoodClusters_goodmembers(counter).ZS,1)),:),[-1 3]);colormap hot
    set(gca,'XTickLabel',[]);set(gca,'YTickLabel',[]);    
    counter=counter+1;
end
print(Fighandle,'D:\Pictures\processed\Flow\Basic\Neomycin\FullLength_Neo_clusters.svg','-dsvg','-r0');    

%--------------------------------
load('__BasicFlow_final_lite.mat')
Neo=load('__BasicFlow_Neo.mat')

%Make figure of full length + STD + rasterhisto
Fighandle=figure;
set(Fighandle, 'Position', [100, 100, 1200, 1200]);
set(findall(Fighandle,'type','text'),'fontSize',12,'fontWeight','bold','FontName','Arial')
counter=1;counter2=1;xplot=4;yplot=2;
StimLength=655;
x = linspace(0.2,StimLength/5,StimLength);
ha = tight_subplot(xplot,yplot,[.01 .01],[.01 .01],[.01 .01]);
for i=[5 7]   
    temp=GoodClusters_goodmembers(i).mean;
    std_temp=GoodClusters_goodmembers(i).STD;        
    %subplot(xplot,yplot,(2*counter)-1);
    axes(ha((2*counter)-1));
    H=shadedErrorBar(x, temp, std_temp);axis([0 130 -1 4]);
    H.mainLine.Color=colors(i,:)/256;
    H.patch.FaceColor=colors(i,:)/256;
    H.edge(1).Color=colors(i,:)/256;
    H.edge(2).Color=colors(i,:)/256;
    H.mainLine.LineWidth=3;    
    set(gca,'XTickLabel',[]);set(gca,'YTickLabel',[]);%title(size(GoodClusters_goodmembers(i).ZS,1));
    %subplot(xplot,yplot,(2*counter));
    axes(ha((2*counter)));
    imagesc(GoodClusters_goodmembers(i).ZS(randperm(size(GoodClusters_goodmembers(i).ZS,1)),:),[-1 3]);colormap hot
    set(gca,'XTickLabel',[]);set(gca,'YTickLabel',[]);    
    counter=counter+1;
    temp=Neo.GoodClusters_goodmembers(counter2).mean;
    std_temp=Neo.GoodClusters_goodmembers(counter2).STD;        
    %subplot(xplot,yplot,(2*counter)-1);
    axes(ha((2*counter)-1));
    H=shadedErrorBar(x, temp, std_temp);axis([0 130 -1 4]);
    H.mainLine.Color=colors(counter,:)/256;
    H.patch.FaceColor=colors(counter,:)/256;
    H.edge(1).Color=colors(counter,:)/256;
    H.edge(2).Color=colors(counter,:)/256;
    H.mainLine.LineWidth=3;    %title(size(Neo.GoodClusters_goodmembers(counter2).ZS,1));
    set(gca,'XTickLabel',[]);set(gca,'YTickLabel',[]);  
    %subplot(xplot,yplot,(2*counter));
    axes(ha((2*counter)));
    imagesc(Neo.GoodClusters_goodmembers(counter2).ZS(randperm(size(Neo.GoodClusters_goodmembers(counter2).ZS,1)),:),[-1 3]);colormap hot
    set(gca,'XTickLabel',[]);set(gca,'YTickLabel',[]);    
    counter=counter+1;
    counter2=counter2+1;
end
print(Fighandle,strcat('D:\Pictures\processed\Flow\Basic\Figure\FullLength_neoAndbasic_clusters.svg'),'-dsvg','-r0');    



%Make figure of full length + STD + rasterhisto
counter=1;counter2=1;
for i=[5 7]   
    size(GoodClusters_goodmembers(i).ZS,1)/length(unique(GoodClusters_goodmembers(i).Fish))
    idx_temp=Neo.GoodClusters_goodmembers(counter2).idx;
    size(Neo.GoodClusters_goodmembers(counter2).ZS,1)/length(unique(Neo.idx_Fish(idx_temp)))
    counter2=counter2+1;
end

