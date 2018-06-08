load('Fish201709131_ERO_10_output_analysis_matlab.mat');
load('Fish201709131_ERO_10_output_correlation.mat');

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
 imwrite(uint16(image),'ROIs_09131_ERO_11.tif','tif');
 
imageSizeY = size(Correlation_image,1);
imageSizeX = size(Correlation_image,2);
[columnsInImage rowsInImage] = meshgrid(1:imageSizeX, 1:imageSizeY);radius =3;
image=zeros(size(Correlation_image));
GoodROIS=ROIs(:,idx_components+1);
for roi_nb=1:size(GoodROIS,2)
    ROI=reshape(full(GoodROIS(:,roi_nb)),imageSizeY,imageSizeX);%reshape cause of resize
    temp=regionprops(uint16(squeeze(ROI))==max(max(uint16(squeeze(ROI)))),'Centroid');
    temp=temp.Centroid;
    xcoord=uint16(temp(1));
    ycoord=uint16(temp(2));
    if xcoord>1 & ycoord>1        
        image([ycoord-1 ycoord ycoord+1],[xcoord-1 xcoord xcoord])=roi_nb;
    end
end
 imwrite(uint16(image),'GoodROIS_09131_ERO_11.tif','tif');

BoundingBox=[270 140;280 140];
 image=zeros(size(Correlation_image));
 IndexInsideBB=[];
 for roi_nb=1:size(ROIs,2)
    ROI=reshape(full(ROIs(:,roi_nb)),imageSizeY,imageSizeX);%reshape cause of resize
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
 imwrite(uint16(image),'ROISFull_09131_ERO_11.tif','tif');
 
 BoundingBox=[270 140;280 140];
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
 imwrite(uint16(image),'GoodROISFull_09131_ERO_11.tif','tif');
 
%correction for 09111_ERO et 09112_ELO
ZS_temp=zscore(DenoisedTraces(idx_components(IndexInsideBB)+1,:)+Noise(idx_components(IndexInsideBB)+1,:),1,2);
ZS_denoised=zscore(DenoisedTraces(idx_components(IndexInsideBB)+1,:),1,2);
%ZS_temp=ZS_temp(:,21:end);

Stimuli=zeros(6,size(ZS_temp,2));
GCaMP6=[0,1.69644104899772,5.13796058542217,8.27886020152244,10.3756715204800,11.8173714529814,12.2425184714093,10.8571417354877,8.80831829681196,6.91339112244670,5.46959264663869,4.30868766622567,3.42533619066766,2.75378443486879,2.18017250852183,1.72816235135824,1.32732537295463,1.00684435500268,0.730210038304555,0.530242444093118,0.362253250339685,0.227668255288566,0.0869242416152502,0.000718266708050853,-0.0828334873368325]';
idxStart=40;
counter=0;
for i=1:18
    if mod(i,3)==1
        counter=counter+1;
    end
    Stimuli(counter,(idxStart+(i-1)*40):(idxStart+(i-1)*40)+size(GCaMP6,1)-1)=GCaMP6;
end


ModelMultipower2=[];
parfor i=1:size(ZS_temp,1)    
    mdl=fitlm(Stimuli',ZS_temp(i,:));%,'interactions');
    ModelMultipower2(i).coef=mdl.Coefficients;    
    ModelMultipower2(i).rsquared=mdl.Rsquared.Adjusted;
end
idx_rsq=find([ModelMultipower2.rsquared]>0.2);
figure;plot(ZS_temp(idx_rsq,:)');

figure;
for i=1:length(idx_rsq)
    subplot(length(idx_rsq),1,i);plot(ZS_temp(idx_rsq(i),:));
end

Fighandle=figure;
set(Fighandle, 'Position', [100, 100, 1600, 800]);
imagesc(ZS_temp(:,1:700),[0 3]);h=colorbar;set(h,'Ticklabels',[]);
set(gca,'XTickLabel',[]);set(gca,'YTickLabel',[]);

idx_select=[50 58];

figure;plot(ZS_temp(idx_select,:)');xlim([0 700]);

colors=[0.14 1 0.14; 0.7 0.4 1];
counter=1;
Fighandle=figure;
set(Fighandle, 'Position', [100, 100, 1000, 400]);
for i=idx_select
    subplot(length(idx_select),1,counter);
    plot(ZS_denoised(i,:),'color',colors(counter,:),'LineWidth',3);xlim([0 700]);ylim([-2 8]);set(gca,'XTickLabel',[]);set(gca,'YTickLabel',[]);
    counter=counter+1;    
end

colors=[0.14 1 0.14; 0.7 0.4 1];
counter=1;
Fighandle=figure;
set(Fighandle, 'Position', [100, 100, 1000, 400]);
for i=idx_select
    subplot(length(idx_select),1,counter);
    plot(ZS_denoised(i,:),'color',colors(counter,:),'LineWidth',3);hold on;plot(zscore(RAW_Fluorescence(21:end,counter)),'color','k','LineWidth',1);
    xlim([0 700]);ylim([-2 8]);set(gca,'XTickLabel',[]);set(gca,'YTickLabel',[]);
    counter=counter+1;    
end


colors=[0.14 1 0.14; 0.7 0.4 1];
counter=1;
Fighandle=figure;
set(Fighandle, 'Position', [100, 100, 800, 400]);
for i=idx_select
    subplot(length(idx_select),1,counter);
    plot(RAW_Fluorescence(:,counter),'color',colors(counter,:),'LineWidth',3);xlim([0 700]);set(gca,'XTickLabel',[]);set(gca,'YTickLabel',[]);
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
 imwrite(uint16(image),'SelectROIs_09131_ERO_11.tif','tif');
 
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
 
load('D:\Pictures\processed\Itia\New\_aMultipower.mat'); 
NumbForSup=[MatFiles(73).GoodNumber MatFiles(74).GoodNumber-1];
ZS_temp=ZS2(NumbForSup(1):NumbForSup(2),:);
figure;plot(ZS_temp(IndexInsideBB(50)+1,:));hold on;plot(ZS_temp(IndexInsideBB(58)+1,:));
figure;plot(ZS2(NumbForSup(1)+IndexInsideBB(50),:));
[ModelMultipower2(NumbForSup(1)+IndexInsideBB(50)).rsquared ModelMultipower2(NumbForSup(1)+IndexInsideBB(58)).rsquared]


figure;
coef=ModelMultipower2(NumbForSup(1)+IndexInsideBB(50)).coef;
temp=coef.Properties.RowNames;temp=regexp(temp,'x(\d+)','tokens');
fitted=zeros(1,760);
for coef_idx=1:height(coef)
    if coef.pValue(coef_idx)<0.05
        coefs(coef_idx)=coef.Estimate(coef_idx);
        if coef_idx==1
            fitted=fitted+coefs(coef_idx);
        else
            fitted=fitted+coefs(coef_idx)*Stimuli(coef_idx-1,:);
        end
    else
        coefs(coef_idx)=0;
    end
end
plot(ZS2(NumbForSup(1)+IndexInsideBB(50),:),'color',colors(1,:),'LineWidth',3);hold on;plot(fitted,'color','k','LineWidth',1);


coefficients={};
for idx=1:length(ModelMultipower2)
    coef=[ModelMultipower2(idx).coef];
    temp=coef.Properties.RowNames;temp=regexp(temp,'x(\d+)','tokens');
    if ~isempty(temp)
        %temp=[temp{:}];temp=[temp{:}];temp=[temp{:}];%temp=str2num(temp);
        for coef_idx=2:height(coef)
            if coef.pValue(coef_idx)<0.05
                coefficients{idx,str2num(temp{coef_idx}{1}{1})}=coef.Estimate(coef_idx);
                endexit
        end
    end
end
idxempty=cellfun('isempty',coefficients);
coefficients(idxempty)={0};
clearvars idxempty idx coef_idx coef temp
coefficients=cell2mat(coefficients);
rsq_list=[ModelMultipower2.rsquared];


colors=[0.14 1 0.14; 0.7 0.4 1];
counter=1;
Fighandle=figure;
set(Fighandle, 'Position', [100, 100, 1000, 400]);
for i=idx_select
    coef=ModelMultipower2(NumbForSup(1)+IndexInsideBB(i)).coef;
    temp=coef.Properties.RowNames;temp=regexp(temp,'x(\d+)','tokens');
    fitted=zeros(1,760);
    for coef_idx=1:height(coef)
        if coef.pValue(coef_idx)<0.05
            coefs(coef_idx)=coef.Estimate(coef_idx);
            if coef_idx==1
                fitted=fitted+coefs(coef_idx);
            else
                fitted=fitted+coefs(coef_idx)*Stimuli(coef_idx-1,:);
            end
        else
            coefs(coef_idx)=0;
        end
    end
    subplot(length(idx_select),1,counter);
    plot(ZS2(NumbForSup(1)+IndexInsideBB(i),:),'color',colors(counter,:),'LineWidth',3);hold on;plot(fitted,'color','k','LineWidth',1);
    xlim([0 700]);ylim([-2 8]);set(gca,'XTickLabel',[]);set(gca,'YTickLabel',[]);
    counter=counter+1;
end

idxKmeans_final_goodmemberInBrain_merge(NumbForSup(1)+IndexInsideBB(idx_select(1)))

Fighandle=figure;
set(Fighandle, 'Position', [100, 100, 1600, 800]);
imagesc(ZS2(find(idxKmeans_final_goodmemberInBrain_merge==idxKmeans_final_goodmemberInBrain_merge(NumbForSup(1)+IndexInsideBB(idx_select(1)))),:),[0 3]);
h=colorbar;set(h,'Ticklabels',[]);set(gca,'XTickLabel',[]);set(gca,'YTickLabel',[]);

colors=[0.14 1 0.14; 0.7 0.4 1];
counter=1;
Fighandle=figure;
set(Fighandle, 'Position', [100, 100, 1000, 200]);
plot(mean(ZS2(find(idxKmeans_final_goodmemberInBrain_merge==idxKmeans_final_goodmemberInBrain_merge(NumbForSup(1)+IndexInsideBB(idx_select(1)))),:),1),'color',colors(counter,:),'LineWidth',3);xlim([0 760]);
set(gca,'XTickLabel',[]);set(gca,'YTickLabel',[]);


Fighandle=figure;
set(Fighandle, 'Position', [100, 100, 1000, 400]);
imagesc(Stimuli,[0 10]);
set(gca,'XTickLabel',[]);set(gca,'YTickLabel',[]);%h=colorbar;set(h,'Ticklabels',[]);

Fighandle=figure;
set(Fighandle, 'Position', [100, 100, 1000, 600]);
x = linspace(0.25,760/4,760);
ha = tight_subplot(3,1,[.01 .01],[.01 .01],[.01 .01]);
colors=[0.14 1 0.14; 0.7 0.4 1;0 0.6 0.6];
for Beta=1:length(GoodBetas_merge)    
    axes(ha(Beta));
    idx_temp2=find(idxKmeans_final_goodmemberInBrain_merge==GoodBetas_merge(Beta));
    plot(detrend(mean(ZS2(idx_temp2,:),1)),'color',colors(Beta,:),'LineWidth',3);xlim([0 760]);ylim([-3 6]); 
    set(gca,'XTickLabel',[]);set(gca,'YTickLabel',[]);
    %print(Fighandle,strcat('__WB_Cluster-Multi',num2str(Beta)),'-dsvg','-r0');
    %close all;
end
Model_ZS2([GoodBetas_merge]).rsquared;

idx_temp=find([Model_ZS2.rsquared]<0.3);

idx_temp=[5 8 15];
Fighandle=figure;
set(Fighandle, 'Position', [100, 100, 1000, 600]);
x = linspace(0.25,760/4,760);
ha = tight_subplot(3,1,[.01 .01],[.01 .01],[.01 .01]);
colors=[0.14 1 0.14; 0.7 0.4 1;0 0.6 0.6];
for Beta=1:length(idx_temp)    
    axes(ha(Beta));
    idx_temp2=find(idxKmeans_final==idx_temp(Beta));
    plot(mean(ZS2(idx_temp2,:),1),'color',colors(Beta,:),'LineWidth',3);xlim([0 760]);ylim([-3 6]); 
    %set(gca,'XTickLabel',[]);set(gca,'YTickLabel',[]);
    %print(Fighandle,strcat('__WB_Cluster-Multi',num2str(Beta)),'-dsvg','-r0');
    %close all;
end
[Model_ZS2(idx_temp).rsquared]

figure;
plot(x,mean(ZS2(idx_temp2,:),1),'color',colors(Beta,:),'LineWidth',3);