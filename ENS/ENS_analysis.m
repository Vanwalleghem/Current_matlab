ZS=zscore(DenoisedTraces(idx_components+1,:),1,2);
ZS_std=std(ZS,1,2);
ZS=ZS(ZS_std>0,:);
%ZS=ZS([1:297 299:431 434:size(ZS,1)],:);
ZS=detrend(ZS')';
framerate=2;
x = linspace(0.2,size(ZS,2)/framerate,size(ZS,2));
y = linspace(1,size(ZS,1)/100,size(ZS,1));

CorrMatrix_dist=pdist(ZS,'correlation');
CorrMatrix=squareform(CorrMatrix_dist);

[RdBu]=cbrewer('div','RdBu',101);
[BuGn]=cbrewer('seq','BuGn',51);

Fighandle=figure;
set(Fighandle, 'Position', [100, 100, 1000, 1000]);set(gca,'xtick',[]);set(gca,'ytick',[]);set(gca,'xticklabel',[]);set(gca,'yticklabel',[]);
imagesc(1-CorrMatrix,[-1 1]);colormap(RdBu);

Fighandle=figure;
set(Fighandle, 'Position', [100, 100, 1000, 500]);
imagesc(x,y,ZS,[-0.5 3]);colormap(hot);%colormap(flipud(BuGn));

ROI_temp=ROIs(:,idx_components+1);
ROI_temp=reshape(full(ROI_temp),[size(Correlation_image) length(idx_components)]);

figure;imagesc(squeeze(sum(ROI_temp,3)));

ROI_GPU=gpuArray(ROI_temp);
ROI_centroid=zeros(size(ROI_temp,3),2);

% ROI_GPU=gpuArray(ROI_temp(:,:,1:500));
for i=1:size(ROI_GPU,3)
    temp=squeeze(ROI_GPU(:,:,i));   
    s=regionprops(temp>0,temp,'WeightedCentroid');
    ROI_centroid(i,:)=s.WeightedCentroid;    
end
% ROI_GPU=gpuArray(ROI_temp(:,:,501:1000));
% for i=1:size(ROI_GPU,3)
%     s=regionprops(ROI_GPU(:,:,i)>0,'WeightedCentroid');
%     ROI_centroid(i+500,:)=s.WeightedCentroid;    
% end
% ROI_GPU=gpuArray(ROI_temp(:,:,1001:1500));
% for i=1:size(ROI_GPU,3)
%     s=regionprops(ROI_GPU(:,:,i)>0,'WeightedCentroid');
%     ROI_centroid(i+1000,:)=s.WeightedCentroid;    
% end
% ROI_GPU=gpuArray(ROI_temp(:,:,1501:end));
% for i=1:size(ROI_GPU,3)
%     s=regionprops(ROI_GPU(:,:,i)>0,'WeightedCentroid');
%     ROI_centroid(i+1500,:)=s.WeightedCentroid;    
% end

figure;scatter(ROI_centroid(:,2),ROI_centroid(:,1));axis([0 size(Correlation_image,2) 0 size(Correlation_image,1)]);

ROI_centroid=ROI_centroid(ZS_std>0,:);
Good_rois = ROI_centroid(:,1)>300;
figure;scatter(ROI_centroid(Good_rois,1),ROI_centroid(Good_rois,2));axis([0 size(Correlation_image,1) 0 size(Correlation_image,1)]);

Correlation_ROIs=pdist(ZS(Good_rois,:),'correlation');

edges=[0:0.05:2];
Fighandle=figure;
%set(Fighandle, 'Position', [10,10, 800, 200]);
histogram(Correlation_ROIs,edges,'normalization','probability','EdgeAlpha',0.5,'EdgeColor','k','FaceColor','g');hold on;
ax = gca;
outerpos = ax.OuterPosition;
ti = ax.TightInset; 
left = outerpos(1) + ti(1);
bottom = outerpos(2) + ti(2);
ax_width = outerpos(3) - ti(1) - ti(3);
ax_height = outerpos(4) - ti(2) - ti(4);
ax.Position = [left bottom ax_width ax_height];

Corr_temp=(1-squareform(Correlation_ROIs));
Corr_temp=weight_conversion(Corr_temp,'autofix');

Bu=cbrewer('div','RdBu',500);
Fighandle=figure;
set(Fighandle, 'Position', [10,10, 500, 500]);
imagesc(Corr_temp,[-1 1]);colormap(Bu)
ax = gca;
outerpos = ax.OuterPosition;
ti = ax.TightInset; 
left = outerpos(1) + ti(1);
bottom = outerpos(2) + ti(2);
ax_width = outerpos(3) - ti(1) - ti(3);
ax_height = outerpos(4) - ti(2) - ti(4);
ax.Position = [left bottom ax_width ax_height];

Correlation_matrix=threshold_absolute(Corr_temp,0.75);
Fighandle=figure;
graph_temp=graph(Correlation_matrix,'omitselfloops');
set(Fighandle, 'Position', [10,10, 1000, 1000]);
imagesc(Correlation_image);hold on;
%scatter(ROI_centroid(~test,2),ROI_centroid(~test,1),20,'g','filled');hold on;
LWidths = 1*graph_temp.Edges.Weight/max(graph_temp.Edges.Weight);
LWidths(~isfinite(LWidths))=0.01;
LWidths(LWidths==0)=0.01;
LColors = ones(length(LWidths),3);LColors=LColors-LWidths;
plot(graph_temp,'EdgeColor','k','LineStyle','-','NodeColor','m','LineWidth',LWidths,'MarkerSize',5,'Marker','o','MarkerSize',2,'XData',ROI_centroid(Good_rois,1),'YData',ROI_centroid(Good_rois,2),'NodeLabel',char.empty(100,0));hold off;
%axis([20 size(Correlation_image,1)-20 20 size(Correlation_image,1)-20]);
ax = gca;
outerpos = ax.OuterPosition;
ti = ax.TightInset; 
left = outerpos(1) + ti(1);
bottom = outerpos(2) + ti(2);
ax_width = outerpos(3) - ti(1) - ti(3);
ax_height = outerpos(4) - ti(2) - ti(4);
ax.Position = [left bottom ax_width ax_height];
set(gca,'Visible','off')

Baseline=cell2mat(Baseline);
Baseline=Baseline(idx_components+1);
Traces=DenoisedTraces(idx_components+1,:);
Traces=Traces(Good_rois,:)-repmat(Baseline(Good_rois)',1,size(Traces,2));
X=Traces;
L=40;
tic
Ws = {};
Hs = {};
numfits = 5; %number of fits to compare
for k = 1:10
    display(sprintf('running seqNMF with K = %i',k))
    for ii = 1:numfits
        [Ws{ii,k},Hs{ii,k}] = seqNMF(X,'K',k, 'L', L,'lambda', 0,'maxiter',30,'showplot',0); 
        % note that max iter set low (30iter) for speed in demo (not recommended in practice)
    end
    inds = nchoosek(1:numfits,2);
    for i = 1:size(inds,1) % consider using parfor for larger numfits
            Diss(i,k) = helper.DISSX(Hs{inds(i,1),k},Ws{inds(i,1),k},Hs{inds(i,2),k},Ws{inds(i,2),k});
    end
    
end

%% Plot Diss and choose K with the minimum average diss.
figure;
plot(1:20,Diss,'ko'), hold on
h1 = plot(1:20,median(Diss,1),'k-','linewidth',2);
h2 = plot([3,3],[0,0.5],'r--');
legend([h1 h2], {'median Diss','true K'})
xlabel('K')
ylabel('Diss')

%% Procedure for choosing lambda
nLambdas = 40; % increase if you're patient
K = 8; 
lambdas = sort([logspace(-1,-5,nLambdas)], 'ascend'); 
loadings = [];
regularization = []; 
cost = []; 
for li = 1:length(lambdas)
    [N,T] = size(X);
    [W, H, ~,loadings(li,:),power]= seqNMF(X,'K',K,'L',L,...
        'lambdaL1W', .1, 'lambda', lambdas(li), 'maxiter', 100, 'showPlot', 0); 
    [cost(li),regularization(li),~] = helper.get_seqNMF_cost(X,W,H);
    display(['Testing lambda ' num2str(li) '/' num2str(length(lambdas))])
end
%% plot costs as a function of lambda
figure;
windowSize = 3; 
b = (1/windowSize)*ones(1,windowSize);
a = 1;
Rs = filtfilt(b,a,regularization); 
minRs = prctile(regularization,10); maxRs= prctile(regularization,90);
Rs = (Rs-minRs)/(maxRs-minRs); 
R = (regularization-minRs)/(maxRs-minRs); 
Cs = filtfilt(b,a,cost); 
minCs =  prctile(cost,10); maxCs =  prctile(cost,90); 
Cs = (Cs -minCs)/(maxCs-minCs); 
C = (cost -minCs)/(maxCs-minCs); 

clf; hold on
plot(lambdas,Rs, 'b')
plot(lambdas,Cs,'r')
scatter(lambdas, R, 'b', 'markerfacecolor', 'flat');
scatter(lambdas, C, 'r', 'markerfacecolor', 'flat');
xlabel('Lambda'); ylabel('Cost (au)')
set(legend('Correlation cost', 'Reconstruction cost'), 'Box', 'on')
set(gca, 'xscale', 'log', 'ytick', [], 'color', 'none')
set(gca,'color','none','tickdir','out','ticklength', [0.025, 0.025])

%%

ops.nCall=[3 20];
ZS_ENS=ZS(Good_rois,:);
[isort1, isort2, Sm] = mapTmap(ZS_ENS,ops);

DarkGrey=zeros(100,3);
for i=2:100
    DarkGrey(i,:)=[i/100 i/100 i/100];
end

Fighandle=figure;
set(Fighandle, 'Position', [100, 100, 1275, 495]);set(gca,'xtick',[]);set(gca,'ytick',[]);set(gca,'xticklabel',[]);set(gca,'yticklabel',[]);
imagesc(Sm(:,1:400),[0,3]);colormap(DarkGrey);
imagesc(ZS_ENS(isort1,1:400), [0 3])
print(Fighandle,strcat('D:\Dropbox\Projets\Jobs\Figures\','RastermapENS'),'-dtiff','-r0');

pks_smooth=cell(1,size(Sm,1));
locs_smooth=cell(1,size(Sm,1));
pks_raw=cell(1,size(Sm,1));
locs_raw=cell(1,size(Sm,1));
for ij=1:size(Sm,1)
    [pks_smooth{ij},locs_smooth{ij}] = findpeaks(Sm(ij,:),'MinPeakProminence',1,'MinPeakDistance',10);
    [pks_raw{ij},locs_raw{ij}] = findpeaks(ZS_ENS(isort1(ij),:),'MinPeakProminence',1,'MinPeakDistance',10);
end

Periodicity_smooth=zeros(size(Sm,1),1);
Periodicity_raw=zeros(size(Sm,1),1);
for ij=1:size(Sm,1)
    locs_temp=locs_smooth{ij}(2:end)-locs_smooth{ij}(1:end-1);
    locs_temp=locs_temp/2;
    Periodicity_smooth(ij)=mean(locs_temp);
    locs_temp=locs_raw{ij}(2:end)-locs_raw{ij}(1:end-1);
    locs_temp=locs_temp/2;
    Periodicity_raw(ij)=mean(locs_temp);
end

ROI_centroid_temp=ROI_centroid(Good_rois,:);
Timing_smooth=nan(size(Sm,1),1);
Timing_raw=nan(size(Sm,1),1);
for ij=1:size(Sm,1)
    Timing_smooth(ij)=locs_smooth{ij}(1);
    if locs_raw{ij}
        Timing_raw(ij)=locs_raw{ij}(1);
    end
end
Timing_smooth=Timing_smooth-min(Timing_smooth);
Timing_raw=Timing_raw-min(Timing_raw);

Velocity_smooth=1.5*(ROI_centroid_temp(isort1,1)-min(ROI_centroid_temp(isort1,1)))./Timing_smooth;
alt=pdist(ROI_centroid_temp(isort1,1))./pdist(Timing_smooth);alt=alt*1.5;alt(~isfinite(alt))=[];


%%

ROI_figures=zeros(size(Correlation_image));
ROI_ENS=ROI_temp(:,:,Good_rois);
for i=1:size(ROI_ENS,3)
    temp=squeeze(ROI_ENS(:,:,isort1(i)));
    ROI_figures=ROI_figures+double(temp>mean(mean(temp(temp>0))+1*std(temp(:))))*i;
end
    
Fighandle=figure;
set(Fighandle, 'Position', [100, 100, size(Correlation_image,2), size(Correlation_image,1)]);set(gca,'xtick',[]);set(gca,'ytick',[]);set(gca,'xticklabel',[]);set(gca,'yticklabel',[]);
imagesc(ROI_figures);colormap viridis

ROI_centroid_temp=ROI_centroid(Good_rois,:);
mean_img=imread('MAX_MAX_GV_ENS_SEM_20200408_fish3_SL23_2Hz_range120_step5_exposure20_power80-1.tif');
best_img=imread('HD_fish3_ENS_Gilles.tif');
Bu=cbrewer('div','RdYlBu',length(isort1));
Fighandle=figure;
set(Fighandle, 'Position', [10,10, size(Correlation_image,2)*2, size(Correlation_image,1)*2]);
imagesc(mean_img);colormap Gray; hold on;
scatter(ROI_centroid_temp(isort1,1),ROI_centroid_temp(isort1,2),30,Bu,'o','filled');hold off;
print(Fighandle,strcat('D:\Dropbox\Projets\Jobs\Figures\','ScatterENS'),'-dtiff','-r0');
Fighandle=figure;colormap(Bu);
colorbar
print(Fighandle,strcat('D:\Dropbox\Projets\Jobs\Figures\','colorbarENS'),'-dtiff','-r0');

Fighandle=figure;
set(Fighandle, 'Position', [10,10, size(Correlation_image,2)*2, size(Correlation_image,1)*2]);
imagesc(best_img);colormap Gray; hold on;
scatter(ROI_centroid_temp(isort1,1)*4,ROI_centroid_temp(isort1,2)*4,30,Bu,'o','filled');hold off;

% [iclustup, isort, Vout] = activityMap(X2);
% figure;imagesc(ZS_ENS(isort1,:), [0 3])

%%
Fighandle=figure;
Time_end=400;
set(Fighandle, 'Position', [100, 100, 1275, 495]);set(gca,'xtick',[]);set(gca,'ytick',[]);set(gca,'xticklabel',[]);set(gca,'yticklabel',[]);
ax = gca;
outerpos = ax.OuterPosition;
ti = ax.TightInset; 
left = outerpos(1) + ti(1);
bottom = outerpos(2) + ti(2);
ax_width = outerpos(3) - ti(1) - ti(3);
ax_height = outerpos(4) - ti(2) - ti(4);
ax.Position = [left bottom ax_width ax_height];
imagesc(ZS_ENS(isort1,1:Time_end), [-1,5]);colormap inferno;
print(Fighandle,strcat('D:\Dropbox\Projets\Jobs\Figures\','RasterENSb'),'-dtiff','-r0');

colorbar;
print(Fighandle,strcat('D:\Dropbox\Projets\Jobs\Figures\','ColorbarENSb'),'-dtiff','-r0');