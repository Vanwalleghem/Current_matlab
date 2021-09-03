for i=1:length(IdealSpikes_bin)    
    Spikes2=Suite2pSpikes{i,1};
    temp=conv2(1,gausswin(5),Spikes2,'Same');
    temp2=conv2(1,gausswin(5),IdealSpikes_bin{i},'Same');
    Suite2pSpikes{i,2}=pdist2(temp,temp2,'correlation');
    Suite2pSpikesResults(i).Neg_correl_max=max(1-Suite2pSpikes{i,2}(IdealNegIdx{i},:),[],2);
    temp=diag(Suite2pSpikes{i,2});
    Suite2pSpikesResults(i).Neg_correl=1-temp(IdealNegIdx{i});
    Suite2pSpikesResults(i).Pos_correl_max=max(1-Suite2pSpikes{i,2}(setdiff(1:end,IdealNegIdx{i}),:),[],2);
    Suite2pSpikesResults(i).Pos_correl=1-temp(setdiff(1:end,IdealNegIdx{i}));
    Suite2pSpikesResults(i).mean_correl=mean(diag(Suite2pSpikes{i,2}));
end

CascadeSpikesResults=struct();
CascadeFiles=dir('pred*_DF.mat');
for i=1:length(IdealSpikes_bin)    
    CascadeSpikes{i,2}=pdist2(conv2(1,gausswin(5),CascadeSpikes{i,1},'Same'),conv2(1,gausswin(5),IdealSpikes_bin{i}(:,33:end-32),'Same'),'correlation');
    CascadeSpikesResults(i).Neg_correl_max=nanmax(1-CascadeSpikes{i,2}(IdealNegIdx{i},:),[],2);
    temp=diag(CascadeSpikes{i,2});
    CascadeSpikesResults(i).Neg_correl=1-temp(IdealNegIdx{i});
    CascadeSpikesResults(i).Pos_correl_max=nanmax(1-CascadeSpikes{i,2}(setdiff(1:end,IdealNegIdx{i}),:),[],2);
    CascadeSpikesResults(i).Pos_correl=1-temp(setdiff(1:end,IdealNegIdx{i}));
    CascadeSpikesResults(i).mean_correl=nanmean(diag(CascadeSpikes{i,2}));    
end

CaImAnSpikes={};
CaImAnSpikesResults=struct();
for i=1:length(IdealSpikes_bin)
    FluorescentTraces=IdealResponses{i}(1:end-1,:);
    SpikeInfer=zeros(size(FluorescentTraces));
    for ij=1:size(FluorescentTraces,1)
        [~, temp, ~]=deconvolveCa(FluorescentTraces(ij,:)','ar1','foopsi');
        SpikeInfer(ij,:)=temp';
    end
    CaImAnSpikes{i,1}=SpikeInfer;
    CaImAnSpikes{i,2}=pdist2(conv2(1,gausswin(5),SpikeInfer,'Same'),conv2(1,gausswin(5),IdealSpikes_bin{i},'Same'),'correlation');
    CaImAnSpikesResults(i).Neg_correl_max=max(1-CaImAnSpikes{i,2}(IdealNegIdx{i},:),[],2);
    temp=diag(CaImAnSpikes{i,2});
    CaImAnSpikesResults(i).Neg_correl=1-temp(IdealNegIdx{i});
    CaImAnSpikesResults(i).Pos_correl_max=max(1-CaImAnSpikes{i,2}(setdiff(1:end,IdealNegIdx{i}),:),[],2);
    CaImAnSpikesResults(i).Pos_correl=1-temp(setdiff(1:end,IdealNegIdx{i}));
    CaImAnSpikesResults(i).mean_correl=mean(diag(CaImAnSpikes{i,2}));
end

CellSortSpikes={};
CellSortSpikesResults=struct();
for i=1:length(IdealSpikes_bin)
    CellSortSpikes{i}=pdist2(conv2(1,gausswin(5),single(full(CellsortFindspikes(zscore(IdealResponses{i}(1:end-1,:),1,2),2,0.2,2,1))'),'Same'),conv2(1,gausswin(5),IdealSpikes_bin{i},'Same'),'correlation');
    CellSortSpikesResults(i).Neg_correl_max=max(1-CellSortSpikes{i}(IdealNegIdx{i},:),[],2);
    temp=diag(CellSortSpikes{i});
    CellSortSpikesResults(i).Neg_correl=1-temp(IdealNegIdx{i});
    CellSortSpikesResults(i).Pos_correl_max=max(1-CellSortSpikes{i}(setdiff(1:end,IdealNegIdx{i}),:),[],2);
    CellSortSpikesResults(i).Pos_correl=1-temp(setdiff(1:end,IdealNegIdx{i}));
    CellSortSpikesResults(i).mean_correl=mean(diag(CellSortSpikes{i}));
end

MeanCorrel=[];
for i=1:length(IdealSpikes_bin)
    MeanCorrel(i,2)=mean(CellSortSpikesResults(i).Pos_correl);
    MeanCorrel(i,1)=mean(CaImAnSpikesResults(i).Pos_correl);
    MeanCorrel(i,3)=mean(Suite2pSpikesResults(i).Pos_correl);
    MeanCorrel(i,4)=mean(CascadeSpikesResults(i).Pos_correl);
    MeanCorrel(i,7)=mean(CellSortSpikesResults(i).Neg_correl);
    MeanCorrel(i,6)=mean(CaImAnSpikesResults(i).Neg_correl);
    MeanCorrel(i,8)=mean(Suite2pSpikesResults(i).Neg_correl);
    MeanCorrel(i,9)=mean(CascadeSpikesResults(i).Neg_correl);
end



Fighandle=figure;
set(Fighandle, 'Position', [100, 100, 600, 600]);
plot(x,zscore(temp(1,:)),'color',[166 33 255]/255,'LineWidth',2);
hold on;plot(x,zscore(temp2(1,:)));