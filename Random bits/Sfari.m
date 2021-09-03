%Fish per Fish for pretty figure
load('D:\Pictures\processed\ENS\MAX_GV_ENS_fish1_2_7DPF_range135_step5_exposure20_power60_output_analysis_matlab.mat')
load('D:\Pictures\processed\ENS\MAX_GV_ENS_fish1_2_7DPF_range135_step5_exposure20_power60_output_correlation.mat')
image_gut=(imread('MAX_GV_ENS_fish1_2_7DPF_range135_step5_exposure20_power60_mean.tiff'));
ZS=zscore(DenoisedTraces(idx_components+1,:),1,2);
ZS_std=std(ZS,1,2);
ZS=ZS(ZS_std>0,:);
ZS=ZS([1:297 299:431 434:size(ZS,1)],:);
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
print(Fighandle,strcat('D:\Pictures\ENS\Figures\','CorrMatrix_4DPF.svg'),'-dsvg','-r0');

Fighandle=figure;
set(Fighandle, 'Position', [100, 100, 1000, 500]);
imagesc(x,y,ZS,[-0.5 3]);colormap(hot);%colormap(flipud(BuGn));
print(Fighandle,strcat('D:\Pictures\ENS\Figures\','ZS_raster_4DPF.svg'),'-dsvg','-r0');

options = statset('UseParallel',1); %parallelize the replicates
[idxKmeans_ZS Cmap_ZS]=kmeans(ZS,5,'Options',options,'Distance','correlation','Replicates',3,'MaxIter',1000,'Display','final');
figure;
imagesc(Cmap_ZS);
figure;
plot(Cmap_ZS');

dims=size(Correlation_image);
Rs=ROIs;
ROI=reshape(full(Rs),dims(1),dims(2),size(Rs,2));counter=1;
ROI_temp=nan(size(ROI,3),2);
for roi_nb=1:size(ROI,3)
    progressbar([],roi_nb/size(ROI,3));
    temp=regionprops(uint16(squeeze(ROI(:,:,roi_nb)))==max(max(uint16(squeeze(ROI(:,:,roi_nb))))),'Centroid');
    temp=temp.Centroid;
    ROI_temp(roi_nb,1:2)=temp;
end
ROI_temp=ROI_temp(idx_components+1,:);
ROI_pool=ROI_temp(ZS_std>0,:);
ROI_pool=ROI_pool([1:297 299:431 434:size(ROI_pool,1)],:);
Dist_matrices=pdist(ROI_pool);
Max_distance_ROIs=max(Dist_matrices);
edges={linspace(0,2,51),linspace(0,500,51)};
CorrVsDist = hist3([CorrMatrix_dist' Dist_matrices'],'Edges',edges);

Fighandle=figure;
set(Fighandle, 'Position', [200,300, 1000, 1000]);
imagesc(CorrVsDist);
ax=gca;
ax.set('XTickLabel',edges{2}(ax.XTick));
ax.YTick=ax.YTick+1;
ax.set('YTickLabel',1-edges{1}(ax.YTick));
print(Fighandle,strcat('D:\Pictures\ENS\Figures\','_CorrelationVsDistance_4DPF.png'),'-dpng','-r0');

Correlation_matrix_binary=CorrMatrix<0.2;
Correlation_matrix_binary=triu(Correlation_matrix_binary,1);
Fighandle=figure;
set(Fighandle, 'Position', [200,300, 1000, 1000]);
imagesc(image_gut);colormap(gray);hold on;
for x_nb=1:size(Correlation_matrix_binary,1)
    for y_nb=1:size(Correlation_matrix_binary,2)
        if Correlation_matrix_binary(x_nb,y_nb)
            hold on;
            plot([ROI_pool(x_nb,1) ROI_pool(y_nb,1)], [ROI_pool(x_nb,2) ROI_pool(y_nb,2)])
        end
    end
end
print(Fighandle,strcat('D:\Pictures\ENS\Figures\','_GraphTheory_3DPF.png'),'-dpng','-r0');

Correlation_matrix_binary=CorrMatrix<0.2;
Graph_temp=graph(Correlation_matrix_binary,'omitselfloops');
Fighandle=figure;
set(Fighandle, 'Position', [200,300, 1000, 1000]);
imagesc(image_gut);colormap gray;hold on;
plot(Graph_temp,'w','XData',ROI_pool(:,1),'YData',ROI_pool(:,2));
Degree_weightGraph = degree(Graph_temp);

Correlation_matrix_graph=1-CorrMatrix;Correlation_matrix_graph(Correlation_matrix_graph<0.8)=0;Correlation_matrix_graph=triu(Correlation_matrix_graph,1);
digraph_temp=digraph(Correlation_matrix_graph,'omitselfloops');
Fighandle=figure;
set(Fighandle, 'Position', [200,300, size(image_gut,2)*2, size(image_gut,1)*2]);
imagesc(image_gut);colormap gray;hold on;
LWidths = 2*digraph_temp.Edges.Weight;
plot(digraph_temp,'w','XData',ROI_pool(:,1),'YData',ROI_pool(:,2),'LineWidth',LWidths);
print(Fighandle,strcat('D:\Pictures\ENS\Figures\','_GraphTheory_weighted.png'),'-dpng','-r0');

temp_peaks={};
temp_locs={};
for j=1:size(ZS,1)
    [temp_peaks{j} temp_locs{j}]= findpeaks(ZS(j,:),'MinPeakDistance',5,'MinPeakHeight',2);
end

GCaMP6=[0,1.69644104899772,5.13796058542217,8.27886020152244,10.3756715204800,12.2425184714093,10.8571417354877,8.80831829681196,6.91339112244670,5.46959264663869,4.30868766622567,3.42533619066766,2.75378443486879,2.18017250852183,1.72816235135824,1.32732537295463,1.00684435500268,0.730210038304555,0.530242444093118,0.362253250339685,0.227668255288566,0.0869242416152502,0.000718266708050853,-0.0828334873368325]';
GCaMP6=GCaMP6(1:16);

Peaks_shape=[];
Peaks_timing={};
for i=1:length(temp_locs)
    Timing_temp=[];
    Rost_temp=[];
    if ~isempty(temp_locs{i})
        for peak_nb=1:length(temp_locs{i})
            idx_peak=temp_locs{i}(peak_nb);
            if idx_peak<(length(ZS(i,:))-10) & idx_peak>5
                temp=ZS(i,idx_peak-5:idx_peak+10)';temp=temp-temp(1);
                Rost_temp=horzcat(Rost_temp,temp/abs((max(temp))-(min(temp))));
                Timing_temp=[Timing_temp idx_peak];
            end
        end
        if ~isempty(Rost_temp)
            test=pdist2(Rost_temp',GCaMP6','correlation');
            Rost_temp=Rost_temp(:,test<0.6);
            Timing_temp=Timing_temp(test<0.6);
            Peaks_timing{i}=Timing_temp;
            Peaks_shape=horzcat(Peaks_shape,Rost_temp);
        end
    end
end

test=cellfun(@length,Peaks_timing);
Gut_ROIs=find(ROI_pool(:,2)>300);
Gut_ROIs=find(ROI_pool(:,2)<120 & ROI_pool(:,2)>80 & ROI_pool(:,1)<275);%3DPF
Gut_ROIs=find(ROI_pool(:,2)<95 & ROI_pool(:,2)>40 & ROI_pool(:,1)<275);%4DPF
test_temp=intersect(Gut_ROIs,find(test>3));

options = statset('UseParallel',1); %parallelize the replicates
[idxKmeans_ZS Cmap_ZS]=kmeans(ZS(test_temp,:),5,'Options',options,'Distance','correlation','Replicates',3,'MaxIter',1000,'Display','final');
figure;
imagesc(Cmap_ZS);
figure;
plot(Cmap_ZS');

[B,index] = sortrows(ZS(test_temp,:)>2,[1:1:size(ZS,2)]);
figure;imagesc(B);

lags_signal=zeros(1,length(index));
for i=1:length(index)
    %[~,~,lags_signal(i)]=alignsignals(ZS(test_temp(i),:),Cmap_ZS(2,:));
    lags_signal(i)=finddelay(ZS(test_temp(i),:),Cmap_ZS(1,:),10);
end

[~,idx_sort]=sort(lags_signal);

Fighandle=figure;
set(Fighandle, 'Position', [200,300, size(image_gut,2)*2, size(image_gut,1)*2]);
imagesc(image_gut);colormap gray;hold on;
[seqTiming]=cbrewer('seq','BuPu',length(idx_sort));
scatter(ROI_pool(test_temp(idx_sort),1),ROI_pool(test_temp(idx_sort),2),16,seqTiming,'filled');
print(Fighandle,strcat('D:\Pictures\ENS\Figures\','DelayedSignal_4DPF.png'),'-dpng','-r0');

Fighandle=figure;
set(Fighandle, 'Position', [200,300, 700, 300]);
imagesc(image_gut);colormap gray;hold on;
[seqTiming]=cbrewer('seq','BuPu',length(index));
scatter(ROI_pool(test_temp(index),1)cl,ROI_pool(test_temp(index),2),16,seqTiming,'filled');

image_gut_welsh=image_gut;
image_gut_welsh(image_gut_welsh>prctile(image_gut(:),99.5))=prctile(image_gut(:),99.5);
CorrMatrix_dist_gut=pdist(ZS(Gut_ROIs,:),'correlation');
CorrMatrix_gut=squareform(CorrMatrix_dist_gut);
Correlation_matrix_graph=1-CorrMatrix_gut;Correlation_matrix_graph(Correlation_matrix_graph<0.8)=0;Correlation_matrix_graph=triu(Correlation_matrix_graph,1);
digraph_temp=digraph(Correlation_matrix_graph,'omitselfloops');
Fighandle=figure;
set(Fighandle, 'Position', [200,300, size(image_gut,2)*2, size(image_gut,1)*2]);
imagesc(image_gut_welsh);colormap gray;hold on;
LWidths = 2*digraph_temp.Edges.Weight;
plot(digraph_temp,'w','XData',ROI_pool(Gut_ROIs,1),'YData',ROI_pool(Gut_ROIs,2),'LineWidth',LWidths);
print(Fighandle,strcat('D:\Pictures\ENS\Figures\','_GraphTheory_weighted_gut_4DPF.png'),'-dpng','-r0');

ROI_full=ROI(:,:,idx_components+1);
ROI_full=ROI_full(:,:,ZS_std>0);
ROI_full=ROI_full(:,:,[1:297 299:431 434:size(ROI_full,3)]);
ROI_full=ROI_full(:,:,test_temp);
ROI_full(ROI_full==0)=0;

selected_ROIs=[3 52 90];
colors={};
colors{1}=cbrewer('seq','Reds',100);
colors{2}=cbrewer('seq','Greens',100);
colors{3}=cbrewer('seq','Blues',100);
colors{1}(1,:)=[0 0 0];colors{2}(1,:)=[0 0 0];colors{3}(1,:)=[0 0 0];

Fighandle=figure;
set(Fighandle, 'Position', [200,300, 1000, 1000]);
imagesc(image_gut);colormap gray;hold on;
print(Fighandle,strcat('D:\Pictures\ENS\Figures\','image_gut.png'),'-dpng','-r0');
Fighandle=figure;
set(Fighandle, 'Position', [200,300, 1000, 1000]);
counter=1;
for ROI_nb=selected_ROIs
    imagesc(squeeze(ROI_full(:,:,ROI_nb)));    
    colormap gray;
	print(Fighandle,strcat('D:\Pictures\ENS\Figures\','image_gut_ROI',num2str(counter),'.png'),'-dpng','-r0');
    counter=counter+1;
end


Fighandle=figure;
set(Fighandle, 'Position', [200,300, 1200, 300]);
colors=[1 0 0 0.5;0 1 0 0.5; 0 0 1 0.5];counter=1;
for ROI_nb=selected_ROIs
    plot(ZS(test_temp(ROI_nb),:),'Color',colors(counter,:),'LineWidth',2.5);hold on;
    counter=counter+1;
end
print(Fighandle,strcat('D:\Pictures\ENS\Figures\','image_gut_ROI_ZS.png'),'-dsvg','-r0');


Gut_ROIs=find(ROI_pool(:,2)<120 & ROI_pool(:,2)>80 & ROI_pool(:,1)<275);

lags_signal=zeros(1,length(index));
for i=1:length(index)
    %[~,~,lags_signal(i)]=alignsignals(ZS(test_temp(i),:),Cmap_ZS(2,:));
    lags_signal(i)=finddelay(ZS(test_temp(i),:),Cmap_ZS(1,:),10);
end

[~,idx_sort]=sort(lags_signal);

Fighandle=figure;
set(Fighandle, 'Position', [200,300, 1000, 1000]);
imagesc(image_gut);colormap gray;hold on;
[seqTiming]=cbrewer('seq','BuPu',length(idx_sort));
scatter(ROI_pool(test_temp(idx_sort),1),ROI_pool(test_temp(idx_sort),2),16,seqTiming,'filled');
print(Fighandle,strcat('D:\Pictures\ENS\Figures\','DelayedSignal.png'),'-dpng','-r0');