load('Fish201709191_ERO_11_output_analysis_matlab.mat')

imageSizeY = size(Correlation_image,1);
imageSizeX = size(Correlation_image,2);
[columnsInImage rowsInImage] = meshgrid(1:imageSizeX, 1:imageSizeY);radius =3;
image=zeros(size(Correlation_image));
for roi_nb=1:size(ROIs,2)
    ROI=reshape(full(ROIs(:,roi_nb)),imageSizeY,imageSizeX);%reshape cause of resize
    temp=regionprops(uint16(squeeze(ROI))==max(max(uint16(squeeze(ROI)))),'Centroid');
    temp=temp.Centroid;
    xcoord=uint16(temp(1));
    ycoord=uint16(temp(2));
    if xcoord>1 & ycoord>1        
        image([ycoord-1 ycoord ycoord+1],[xcoord-1 xcoord xcoord])=roi_nb;
    end
end
 imwrite(uint16(image),'ROIs_09191_ERO_11.tif','tif');

 
 BoundingBox=[260 175;70 145];
 image=zeros(size(Correlation_image));
 IndexInsideBB=[];
 GoodROIS=ROIs(:,idx_components+1);
 for roi_nb=1:size(GoodROIS,2)
    ROI=reshape(full(GoodROIS(:,roi_nb)),imageSizeY,imageSizeX);%reshape cause of resize
    temp=regionprops(uint16(squeeze(ROI))==max(max(uint16(squeeze(ROI)))),'Centroid');
    temp=temp.Centroid;
    xcoord=uint16(temp(1));
    ycoord=uint16(temp(2));
    if xcoord>BoundingBox(1,1) & ycoord>BoundingBox(2,1) & xcoord<BoundingBox(1,1)+BoundingBox(1,2) & ycoord<BoundingBox(2,1)+BoundingBox(2,2)
        %temp=ROI/max(ROI(:));
        image=image+ROI;
        IndexInsideBB=[IndexInsideBB roi_nb];
    end
 end
 figure;imagesc(image);
 
 
%correction for 09111_ERO et 09112_ELO
ZS_temp=zscore(DenoisedTraces(idx_components(IndexInsideBB)+1,:)+Noise(idx_components(IndexInsideBB)+1,:),1,2);
ZS_temp=ZS_temp(:,21:end);

Fighandle=figure;
set(Fighandle, 'Position', [100, 100, 1000, 800]);
imagesc(ZS_temp(:,1:700),[-3 3]);
set(gca,'XTickLabel',[]);set(gca,'YTickLabel',[]);

idx_select=[28 58];
idx_select=[48 58];
figure;plot(ZS_temp(idx_select,:)');xlim([0 700]);

colors=[0.14 1 0.14; 0.7 0.4 1];
counter=1;
Fighandle=figure;
set(Fighandle, 'Position', [100, 100, 1000, 400]);
for i=idx_select
    subplot(length(idx_select),1,counter);
    plot(ZS_temp(i,:),'color',colors(counter,:),'LineWidth',3);xlim([0 700]);ylim([-2 6]);set(gca,'XTickLabel',[]);set(gca,'YTickLabel',[]);
    counter=counter+1;    
end

colors=[0.14 1 0.14; 0.7 0.4 1];
counter=1;
Fighandle=figure;
set(Fighandle, 'Position', [100, 100, 1000, 400]);
for i=idx_select
    subplot(length(idx_select),1,counter);
    plot(RAW_Fluorescence(21:721,counter),'color',colors(counter,:),'LineWidth',3);xlim([0 700]);set(gca,'XTickLabel',[]);set(gca,'YTickLabel',[]);
    counter=counter+1;    
end


% Fighandle=figure;
% set(Fighandle, 'Position', [100, 100, 800, 800]);
% for i=1:size(ZS_sample,1)
%     subplot(size(ZS_sample,1),1,i);
%     plot(ZS_sample(i,:));
% end

image=zeros(size(Correlation_image));
 for roi_nb=IndexInsideBB(idx_select)
    ROI=reshape(full(GoodROIS(:,roi_nb)),imageSizeY,imageSizeX);%reshape cause of resize
    temp=regionprops(uint16(squeeze(ROI))==max(max(uint16(squeeze(ROI)))),'Centroid');
    temp=temp.Centroid;
    xcoord=uint16(temp(1));
    ycoord=uint16(temp(2));
    if xcoord>BoundingBox(1,1) & ycoord>BoundingBox(2,1) & xcoord<BoundingBox(1,1)+BoundingBox(1,2) & ycoord<BoundingBox(2,1)+BoundingBox(2,2)
        temp=ROI;%/max(ROI(:));     
        image=image+temp;        
    end
 end
 figure;imagesc(image);
 imwrite(uint16(image),'SelectROIS_09191_ERO_11.tif','tif');
 
imageSizeY = size(Correlation_image,1);
imageSizeX = size(Correlation_image,2);
[columnsInImage rowsInImage] = meshgrid(1:imageSizeX, 1:imageSizeY);radius =3;
image=zeros(size(Correlation_image));
idx_temp=IndexInsideBB(idx_select);
for roi_nb=idx_temp
    ROI=reshape(full(GoodROIS(:,roi_nb)),imageSizeY,imageSizeX);%reshape cause of resize
    temp=regionprops(uint16(squeeze(ROI))==max(max(uint16(squeeze(ROI)))),'Centroid');
    temp=temp.Centroid;
    xcoord=uint16(temp(1));
    ycoord=uint16(temp(2));
    if xcoord>1 & ycoord>1        
        image([ycoord-1 ycoord ycoord+1],[xcoord-1 xcoord xcoord])=roi_nb;
    end
end
 imwrite(uint16(image),'ROIs_09191_ERO_11.tif','tif');
 
 
 