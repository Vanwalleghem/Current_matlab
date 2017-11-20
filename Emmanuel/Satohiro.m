ROI=reshape(full(ROIs),[size(Correlation_image) size(ROIs,2)]);
centroids=regionprops(BW,'centroid');
parfor i=1:size(ROIs,2)
    BW = imbinarize(ROI(:,:,i));
    centroids=regionprops(BW,'centroid');
    s{i}=centroids.Centroid;
end

temp=zeros(2,size(ROIs,2));
for i=1:size(ROIs,2)
    temp(:,i)=s{i};
end
all_cent=temp;

figure;scatter(temp(1,:),-temp(2,:));
temp=all_cent(:,idx_components+1);
[~,telencephalon]=find(temp(1,:)>110 & temp(1,:)<275 & temp(2,:)<85);
[~,pretectum]=find(temp(1,:)>94 & temp(1,:)<145 & temp(2,:)>116 & temp(2,:)<173);
[~,thalamus]=find(temp(1,:)>155 & temp(1,:)<245 & temp(2,:)>100 & temp(2,:)<200);
[~,hindbrain]=find(temp(1,:)>78 & temp(1,:)<160 & temp(2,:)>249 & temp(2,:)<297);

figure;
imagesc(Correlation_image);hold on;
scatter(temp(1,[pretectum telencephalon thalamus hindbrain]),temp(2,[pretectum telencephalon thalamus hindbrain]));

ModelResults=[];
parfor i=1:length(telencephalon)
    mdl=stepwiselm(Stimuli2',ZS(telencephalon(i),:),'linear','Criterion','adjrsquared','Intercept',false,'Upper','interactions','Verbose',0);
    %mdl=stepwiselm(Stimuli',ZS(i,:),'linear','Criterion','adjrsquared','Upper','interactions','Verbose',0);
    ModelResults(i).rsquared=mdl.Rsquared.Adjusted;
end
rsquare_loom=[ModelResults.rsquared];
idx_rsq=find(rsquare_loom>0.15);
figure;
imagesc(ZS(telencephalon(idx_rsq),:));


[~,idx_100]=datasample(telencephalon,100);

figure;
imagesc(ZS(pretectum,:),[0 4]);colormap hot;

figure;
imagesc(ZS(telencephalon(idx_100),:),[0 4]);colormap hot;

figure;
imagesc(ZS(thalamus,:),[0 4]);colormap hot;

figure;
imagesc(ZS(hindbrain,:),[0 4]);colormap hot;

LoomResponseDataset.preTectum=ZS(pretectum,:);
LoomResponseDataset.telencephalon=ZS(telencephalon(idx_rsq),:);
LoomResponseDataset.thalamus=ZS(thalamus,:);
LoomResponseDataset.hindbrain=ZS(hindbrain,:);
LoomResponseDataset.image=Correlation_image;
LoomResponseDataset.Centroids_preTectum=temp(:,pretectum);
LoomResponseDataset.Centroids_telencephalon=temp(:,telencephalon(idx_rsq));
LoomResponseDataset.Centroids_thalamus=temp(:,thalamus);
LoomResponseDataset.Centroids_hindbrain=temp(:,hindbrain);
LoomResponseDataset.stimuli_timing=Stimuli2;

figure;
imagesc(Correlation_image);hold on;
scatter(temp(1,telencephalon),temp(2,telencephalon));