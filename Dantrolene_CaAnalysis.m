Traces=DenoisedTraces(idx_components+1,:);
figure;imagesc(zscore(Traces,1,2),[-0.5 4]);colormap hot;
ROI_temp=ROIs(:,idx_components+1);
ROI_temp=reshape(full(ROI_temp),[size(Correlation_image) length(idx_components)]);

figure;imagesc(squeeze(sum(ROI_temp,3)));

ROI_GPU=gpuArray(ROI_temp);
ROI_centroid=zeros(size(ROI_temp,3),2);
for i=1:size(ROI_temp,3)
    s=regionprops(ROI_GPU(:,:,i)>0,'centroid');
    ROI_centroid(i,:)=s.Centroid;    
end
figure;scatter(ROI_centroid(:,1),ROI_centroid(:,2));axis([0 size(Correlation_image,1) 0 size(Correlation_image,1)]);

%test=ROI_centroid(:,1)<50 | ROI_centroid(:,1)>length(Correlation_image(:,1))-50 | ROI_centroid(:,2)<50 | ROI_centroid(:,2)>length(Correlation_image(:,1))-50;
test=ROI_centroid(:,1)<30 | ROI_centroid(:,1)>200 | ROI_centroid(:,2)<30 | ROI_centroid(:,2)>200;
% %Piece 8
% Good_rois = ROI_centroid(:,1)>100 & ROI_centroid(:,1)<500 & ROI_centroid(:,2)>100 & ROI_centroid(:,2)<500;

%Piece F63
Good_rois = ROI_centroid(:,1)>20 & ROI_centroid(:,1)<length(Correlation_image)-20 & ROI_centroid(:,2)>20 & ROI_centroid(:,2)<length(Correlation_image)-20;

%Piece 7
Good_rois = ROI_centroid(:,1)>20 & ROI_centroid(:,1)<length(Correlation_image)-20 & ROI_centroid(:,2)>20 & ROI_centroid(:,2)<length(Correlation_image)-20;

%Piece F92
Good_rois = ROI_centroid(:,1)>20 & ROI_centroid(:,1)<length(Correlation_image)-20 & ROI_centroid(:,2)>20 & ROI_centroid(:,2)<length(Correlation_image)-20;

%Piece 8 - Green is actually red...
Good_rois = ROI_centroid(:,1)>200 & ROI_centroid(:,1)<length(Correlation_image)-50 & ROI_centroid(:,2)>100 & ROI_centroid(:,2)<length(Correlation_image)-30;


figure;scatter(ROI_centroid(Good_rois,1),ROI_centroid(Good_rois,2));axis([0 size(Correlation_image,1) 0 size(Correlation_image,1)]);

Correlation_ROIs={};
Correlation_ROIs{1}=pdist(Traces(Good_rois,20:size(Traces,2)-100),'correlation');
Correlation_ROIs{2}=pdist(Traces(Good_rois,size(Traces,2)-100:end),'correlation');
% for i=1:5
%     Correlation_ROIs{i}=pdist(Traces(~test,20+100*(i-1):20+100*(i)),'correlation');    
% end
% Correlation_ROIs{6}=pdist(Traces(~test,607:707),'correlation');    

edges=[0:0.05:2];
Fighandle=figure;
%set(Fighandle, 'Position', [10,10, 800, 200]);
histogram(Correlation_ROIs{1},edges,'normalization','probability','EdgeAlpha',0.5,'EdgeColor','k','FaceColor','g');hold on;histogram(Correlation_ROIs{2},edges,'normalization','probability','EdgeAlpha',0.5,'EdgeColor','k','FaceColor','m');hold off;
ax = gca;
outerpos = ax.OuterPosition;
ti = ax.TightInset; 
left = outerpos(1) + ti(1);
bottom = outerpos(2) + ti(2);
ax_width = outerpos(3) - ti(1) - ti(3);
ax_height = outerpos(4) - ti(2) - ti(4);
ax.Position = [left bottom ax_width ax_height];
print(Fighandle,strcat('D:\Felicity\Dantrolene\Caiman\Figure\','CorrelationHist3'),'-dtiff','-r0');


ROI_centroid=ROI_centroid(~test,:);

Dant(1).Traces=Traces;
Dant(1).Corr=Correlation_ROIs;

Dant(1).ROI_centroids=ROI_centroid;
% 
% 
% Dant(2).Traces=Traces;
% Dant(2).Corr=Correlation_ROIs;
% Dant(2).ROI_centroids=ROI_centroid(~test,:);


Corr_temp=(1-squareform(Correlation_ROIs{1}));
Corr_temp=weight_conversion(Corr_temp,'autofix');
%Corr_temp=weight_conversion(Corr_temp,'normalize');

Bu=cbrewer('div','RdBu',500);
Fighandle=figure;
set(Fighandle, 'Position', [10,10, 500, 500]);
imagesc(Corr_temp,[-1 1]);colorax = gca;
outerpos = ax.OuterPosition;
ti = ax.TightInset; 
left = outerpos(1) + ti(1);
bottom = outerpos(2) + ti(2);
ax_width = outerpos(3) - ti(1) - ti(3);
ax_height = outerpos(4) - ti(2) - ti(4);
ax.Position = [left bottom ax_width ax_height];map(Bu)

set(gca,'Visible','off')
print(Fighandle,strcat('D:\Felicity\Dantrolene\Caiman\Figure\','CorrelationMatrixPref3_pre100'),'-dtiff','-r0');
Correlation_matrix{1}=threshold_absolute(Corr_temp,0.75);

Corr_temp=(1-squareform(Correlation_ROIs{2}));
Corr_temp=weight_conversion(Corr_temp,'autofix');
%Corr_temp=weight_conversion(Corr_temp,'normalize');
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
set(gca,'Visible','off')
print(Fighandle,strcat('D:\Felicity\Dantrolene\Caiman\Figure\','CorrelationMatrixPref3_post100'),'-dtiff','-r0');
Correlation_matrix{2}=threshold_absolute(Corr_temp,0.75);

Green=zeros(500,3);
for i=1:size(Green,1)
    Green(i,2)=(i/size(Green,1));
end

%Backgd_img=imread('C1-Dant_Piece80000_LongReg_3D_green_mean.tiff');
Fighandle=figure;
graph_temp=graph(Correlation_matrix{1},'omitselfloops');
set(Fighandle, 'Position', [10,10, 1000, 1000]);
imagesc(Correlation_image');colormap(Green); hold on;
%scatter(ROI_centroid(~test,2),ROI_centroid(~test,1),20,'g','filled');hold on;
LWidths = 1*graph_temp.Edges.Weight/max(graph_temp.Edges.Weight);
LWidths(~isfinite(LWidths))=0.01;
LWidths(LWidths==0)=0.01;
LColors = ones(length(LWidths),3);LColors=LColors-LWidths;
plot(graph_temp,'EdgeColor','k','LineStyle','-','NodeColor','m','LineWidth',LWidths,'MarkerSize',5,'Marker','o','MarkerSize',2,'XData',ROI_centroid(:,2),'YData',ROI_centroid(:,1),'NodeLabel',char.empty(100,0));hold off;
axis([20 size(Correlation_image,1)-20 20 size(Correlation_image,1)-20]);
ax = gca;
outerpos = ax.OuterPosition;
ti = ax.TightInset; 
left = outerpos(1) + ti(1);
bottom = outerpos(2) + ti(2);
ax_width = outerpos(3) - ti(1) - ti(3);
ax_height = outerpos(4) - ti(2) - ti(4);
ax.Position = [left bottom ax_width ax_height];
set(gca,'Visible','off')
print(Fighandle,strcat('D:\Felicity\Dantrolene\Caiman\Figure\','Graph3_pre100'),'-dtiff','-r0');

Fighandle=figure;
graph_temp=graph(Correlation_matrix{2},'omitselfloops');
set(Fighandle, 'Position', [10,10, 1000, 1000]);
imagesc(Correlation_image');colormap(Green); hold on;
%scatter(ROI_centroid(~test,2),ROI_centroid(~test,1),20,'g','filled');hold on;
LWidths = 1*graph_temp.Edges.Weight/max(graph_temp.Edges.Weight);
LWidths(~isfinite(LWidths))=0.01;
LWidths(LWidths==0)=0.01;
LColors = ones(length(LWidths),3);LColors=LColors-LWidths;
plot(graph_temp,'EdgeColor','k','LineStyle','-','NodeColor','m','LineWidth',LWidths,'MarkerSize',5,'Marker','o','MarkerSize',2,'XData',ROI_centroid(:,2),'YData',ROI_centroid(:,1),'NodeLabel',char.empty(100,0));hold off;
axis([20 size(Correlation_image,1)-20 20 size(Correlation_image,1)-20]);
ax = gca;
outerpos = ax.OuterPosition;
ti = ax.TightInset; 
left = outerpos(1) + ti(1);
bottom = outerpos(2) + ti(2);
ax_width = outerpos(3) - ti(1) - ti(3);
ax_height = outerpos(4) - ti(2) - ti(4);
ax.Position = [left bottom ax_width ax_height];
set(gca,'Visible','off')
print(Fighandle,strcat('D:\Felicity\Dantrolene\Caiman\Figure\','Graph3_post100'),'-dtiff','-r0');
