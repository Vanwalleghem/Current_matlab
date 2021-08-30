edges=[0:0.05:2];
Fighandle=figure;
set(Fighandle, 'Position', [10,10, 1500, 1500]);
WHPlot(W(indSort,:,:),H, X(indSort,:), 1); title('lambdaOrthoH -> events based')
print(Fighandle,strcat('D:\Dropbox\Papers\Felicity_lactation\events_f63.tif'),'-dtiff','-r0');

WHPlot(W(indSort,:,:),H, X(indSort,:), 0); title('lambdaOrthoH -> events based')
print(Fighandle,strcat('D:\Dropbox\Papers\Felicity_lactation\events_f63b.tif'),'-dtiff','-r0');

Fighandle=figure(2);
set(Fighandle, 'Position', [10,10, 1500, 1500]);
yplot=min([size(H,1); length(unique(hybrid(:,2)))]);
ha=tight_subplot(2,yplot);

for cluster_nb=1:min([size(H,1); length(unique(hybrid(:,2)))])
    axes(ha(0+cluster_nb));
    plot(H(cluster_nb,:),'k', 'linewidth', 2);xlim([0 size(X,2)]);
    axes(ha(yplot+cluster_nb));
    imagesc(squeeze(W(indSort,cluster_nb,:)));colormap hot;
end

ROI_centroid_temp=ROI_centroid(Good_rois,:);
Bu=cbrewer('div','RdYlBu',length(indSort));
Fighandle=figure;
set(Fighandle, 'Position', [10,10, size(Correlation_image,1)*2, size(Correlation_image,2)*2]);
imagesc(Correlation_image');colormap Gray; hold on;
scatter(ROI_centroid_temp(indSort,2),ROI_centroid_temp(indSort,1),25,Bu,'o','filled');hold off;

dist_ROIs=pdist(ROI_centroid_temp(indSort,:));
figure;
imagesc(squareform(dist_ROIs));

temp=squareform(dist_ROIs);

%For F63 :
[~,i1]=max(W(indSort(10),1,:))
[~,i2]=max(W(indSort(185),1,:))
temp(10,185)*0.8286408/((i2-i1)*2.917)

velocity=zeros(185-10,1);
for i=1:185-10
    [~,i1]=max(W(indSort(10),1,:));
    [~,i2]=max(W(indSort(i+10),1,:));
	velocity(i)=temp(10,i+10)*0.8286408/((i2-i1)*2.917);
end

velocity=zeros(185-10,1);
counter=1;
for i=1:length(indSort)
    for j=2:length(indSort)
        if j>i
            [~,i1]=max(W(indSort(i),1,:));
            [~,i2]=max(W(indSort(j),1,:));
            if i2>i1
                velocity(counter)=temp(i,j)*0.8286408/((i2-i1)*2.917);
                counter=counter+1;
            end
        end
    end
end


Fighandle=figure;
set(Fighandle, 'Position', [10,10, 1000, 1000],'visible','off');
WHPlot(W_parts(indSort_parts,:,:),H(:,:), X(indSort_parts,:), 1); title('lambdaOrthoW -> parts based')
print(Fighandle,strcat('D:\Dropbox\Papers\Felicity_lactation\parts_f63.png'),'-dpng','-r0');

ha = tight_subplot(3,2,[.01 .01],[.01 .01],[.01 .01]);

for i=1:6
    axes(ha(i));
    histogram(Correlation_ROIs{i},edges,'normalization','probability','EdgeAlpha',0.5,'EdgeColor','k');
end
ax = gca;
outerpos = ax.OuterPosition;
ti = ax.TightInset; 
left = outerpos(1) + ti(1);
bottom = outerpos(2) + ti(2);
ax_width = outerpos(3) - ti(1) - ti(3);
ax_height = outerpos(4) - ti(2) - ti(4);
ax.Position = [left bottom ax_width ax_height];

Baseline=cell2mat(Baseline);
Traces=DenoisedTraces(idx_components+1,:);
Traces=Traces-repmat(Baseline(idx_components+1)',1,size(Traces,2));

[W,H,cost,loadings,power] = seqNMF(Traces,'K',5,'L',50,'lambda',0.01);
clust=2;
figure;findpeaks(H(clust,:),'MinPeakProminence',0.05,'MinPeakDistance',10);
[pks,locs] = findpeaks(H(clust,:),'MinPeakProminence',0.05,'MinPeakDistance',10);
locs_temp=locs(2:end)-locs(1:end-1);
locs_temp=locs_temp*3.598;
%for f63
locs_temp=locs_temp*2.917;

Select_rois = ROI_centroid(:,1)>20 & ROI_centroid(:,1)<300 & ROI_centroid(:,2)>300 & ROI_centroid(:,2)<size(Correlation_image,1)-20;
figure;scatter(ROI_centroid(Select_rois,1),ROI_centroid(Select_rois,2));axis([0 size(Correlation_image,1) 0 size(Correlation_image,1)]);

Baseline=cell2mat(Baseline);
Baseline=Baseline(idx_components+1);

Traces=DenoisedTraces(idx_components+1,:);
Traces=Traces(Select_rois,:)-repmat(Baseline(Select_rois)',1,size(Traces,2));
X=Traces;


Fighandle=figure(2);
set(Fighandle, 'Position', [10,10, 1500, 1500]);
yplot=min([size(H,1); length(unique(hybrid(:,2)))]);
ha=tight_subplot(2,yplot);
for cluster_nb=1:min([size(H,1); length(unique(hybrid(:,2)))])
    axes(ha(0+cluster_nb));
    plot(H(cluster_nb,:),'k', 'linewidth', 2);xlim([0 size(X,2)]);
    axes(ha(yplot+cluster_nb));
    imagesc(squeeze(W(indSort,cluster_nb,:)));colormap hot;
end

figure;
imagesc(X(indSort(hybrid(:,2)==1),:));
figure;
imagesc(X(indSort(hybrid(:,2)==2),:));

ops.nCall=[3 20];
[isort1, isort2, Sm] = mapTmap(X,ops);
DarkGrey=zeros(100,3);
for i=2:100
    DarkGrey(i,:)=[i/100 i/100 i/100];
end

Fighandle=figure;
set(Fighandle, 'Position', [100, 100, 1275, 495]);set(gca,'xtick',[]);set(gca,'ytick',[]);set(gca,'xticklabel',[]);set(gca,'yticklabel',[]);
imagesc(Sm,[0,3]);colormap(DarkGrey);
imagesc(X(isort1,:))

%velocity_all={};
ROI_centroid_temp=ROI_centroid;
ROI_centroid_temp=ROI_centroid(Good_rois,:);
velocity={}; 
Fighandle=figure;
set(Fighandle, 'Position', [10,10, size(Correlation_image,1)*4, size(Correlation_image,2)*4]);
yplot=min([size(H,1); length(unique(hybrid(:,2)))]);
%yplot=2;
fps=3.598;
%fps=2.917;
ha = tight_subplot(yplot,3,[.03 .03],[.03 .03],[.03 .03]);
cluster_order=[2 1 3 4 5];
x=linspace(0,size(X,2)*fps,size(X,2));
for cluster_nb=1:min([size(H,1); length(unique(hybrid(:,2)))])
    idx_temp=find(hybrid(:,2)==cluster_nb);
    lag_temp=hybrid(idx_temp,1);
    idx_sort=hybrid(idx_temp,3);
    roi_temp=ROI_centroid_temp(idx_sort,:);
    Bu=cbrewer('div','RdYlBu',length(idx_sort));
    dist_ROIs=pdist(roi_temp);
    dist_ROIs=squareform(dist_ROIs);
    temp_vel=[];
    counter=1;
    for i=1:length(idx_temp)
        for j=2:length(idx_temp)
            if j>i
                i1=lag_temp(i);
                i2=lag_temp(j);
                if (i2>i1 && dist_ROIs(i,j)*0.8286408 > 3)
                    %temp_vel(counter)=dist_ROIs(i,j)*0.8286408/((i2-i1)*2.917);
                    if dist_ROIs(i,j)*0.8286408 <200
                        temp_vel(counter)=dist_ROIs(i,j)*0.8286408/((i2-i1)*fps);
                        counter=counter+1;
                    end
                end
            end
        end
    end
    velocity{cluster_nb}=temp_vel;
    axes(ha(1+(cluster_order(cluster_nb)-1)*3));
    %cluster_order(cluster_nb)
    plot(x,H(cluster_nb,:),'k', 'linewidth', 2);xlim([0 size(X,2)*fps]);title(strcat('velocity : ',num2str(mean(temp_vel)),' +/- ',num2str(std(temp_vel))));
    %plot(mean(X(idx_temp,:),1),'k', 'linewidth', 2);xlim([0 size(X,2)]);title(strcat('velocity : ',num2str(mean(temp_vel)),' +/- ',num2str(std(temp_vel))));
    axes(ha(2+(cluster_order(cluster_nb)-1)*3));
    imagesc(squeeze(W(idx_sort,cluster_nb,:)));colormap hot;
    set(gca,'Visible','off')
    axes(ha(3+(cluster_order(cluster_nb)-1)*3));
    imagesc(Correlation_image'); hold on;
    scatter(roi_temp(:,2),roi_temp(:,1),50,Bu,'o','filled');hold off;
    set(gca,'Visible','off')     
end
for cluster_nb=1:min([size(H,1); length(unique(hybrid(:,2)))])
    colormap(ha(3+(cluster_nb-1)*3),gray);
end
print(Fighandle,strcat('D:\Dropbox\Papers\Felicity_lactation\velocity_f94_8.tif'),'-dtiff','-r0');
print(Fighandle,strcat('D:\Dropbox\Papers\Felicity_lactation\velocity_f94_8.eps'),'-depsc','-r0');
print(Fighandle,strcat('D:\Dropbox\Papers\Felicity_lactation\velocity_f94_8.svg'),'-dsvg','-r0');
close(Fighandle);