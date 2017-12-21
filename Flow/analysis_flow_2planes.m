Threshold=0.3;
[Model_ZS,GoodBetas]=Test_Regress(Cmap_ZS,flow,idxKmeans_ZS,Threshold);

Corr_BasicCLust=[];
parfor i=1:size(ZS,1)
    corr_temp=[];    
    for j=1:size(Basic_Clusters,1)
        temp=corrcoef(ZS(i,:), Basic_Clusters(j,:));
        corr_temp(j)=temp(1,2);        
    end
    Corr_BasicCLust(:,i)=corr_temp;
end

GoodBetas_select=[1 7 6 10 8 9 5 4];
idxKmeans=zeros(1,size(ZS,1));
idxKmeans(idx_rsq)=idxKmeans_ZS;
Cmap=Cmap_ZS;

Fighandle=figure;
set(Fighandle, 'Position', [100, 100, 1400, 900]);
counter=1;xplot=floor(sqrt(length(GoodBetas_select)));yplot=ceil(length(GoodBetas_select)/xplot);
for i=GoodBetas_select
idx=find(idxKmeans==i);
NumberOfCells=length(idx);
subplot(length(GoodBetas_select),1,counter);plot(Cmap_ZS(i,:));%,'color',colors(counter,:)/256);title(num2str(NumberOfCells))
counter=counter+1;
end

Fighandle=figure;
set(Fighandle, 'Position', [100, 100, 1400, 900]);
counter=1;xplot=floor(sqrt(length(GoodBetas_select)));yplot=ceil(length(GoodBetas_select)/xplot);
x = linspace(0.2,size(Cmap,2)/5,size(Cmap,2));
for i=GoodBetas_select
idx=find(idxKmeans==i);
NumberOfCells=length(idx);
subplot(length(GoodBetas_select),3,counter);plot(x,mean(ZS(idx,:),1));title(num2str(NumberOfCells));axis([0 131 -1 3]);rectangle('FaceColor','r','Position',[10 -1 10 0.2]);rectangle('FaceColor','r','Position',[50 -1 10 0.2]);rectangle('FaceColor','r','Position',[90 -1 10 0.2]);rectangle('FaceColor','b','Position',[30 -1 10 0.2]);rectangle('FaceColor','b','Position',[70 -1 10 0.2]);rectangle('FaceColor','b','Position',[110 -1 10 0.2]);
subplot(length(GoodBetas_select),3,counter+1);imagesc(zscore(GoodCalcium(idx,:),1,2),[-1 3]);colormap hot
subplot(length(GoodBetas_select),3,counter+2);histogram(idx_Plane(idx),[0:20:300]);
counter=counter+3;
end

Fighandle=figure;
set(Fighandle, 'Position', [100, 100, 1400, 900]);
counter=1;xplot=floor(sqrt(length(GoodBetas_select)));yplot=ceil(length(GoodBetas_select)/xplot);
BasicSpikes=zeros(size(Basic_Clusters));
for i=GoodBetas_select
idx=find(idxKmeans==i);
NumberOfCells=length(idx);
subplot(length(GoodBetas_select),1,counter);plot(mean(zscore(GoodSpikes(idx,:),1,2),1));%,'color',colors(counter,:)/256);title(num2str(NumberOfCells))
BasicSpikes(i,:)=mean(zscore(GoodSpikes(idx,:),1,2),1);
counter=counter+1;
end

All_ROIs=[];
ROIs_idx=[];
for i = 1:length(MatFiles)
    name=strcat(MatFiles(i).name);
    Rs=load(name, 'ROIs');
    Rs=Rs.ROIs;
    F=load(name, 'idx_components');
    F=F.idx_components+1;
    Rs=Rs(:,F);
    All_ROIs{i}=Rs;
    if i==1
        ROIs_idx(i)=length(F);
    else
        ROIs_idx(i)=ROIs_idx(i-1)+length(F);
    end
end
clearvars GC C S F N name i;

Numbers=[0 [ROIs_idx]];
temp=[];
counter=1;
for i=GoodBetas_select
    temp{counter}=find(idxKmeans==i);
    %tempidx=find(idxKmeans==idx);
    %temp{counter}=GoodClusters_goodmembers(counter).idx;
    counter=counter+1;    
end

Start=min(cellfun(@min, temp));Start=find(Numbers<Start,1,'last');
filename=MatFiles(Start).name;

%colors = [1,0,0;0,1,0;0,0,1;1,0.600000000000000,0;1,0,0.7];
%colors = distinguishable_colors(9,[1 1 1; 0 0 0]);
colors = [0         0    1.0000
         0    0.5000    1.0000
    1.0000         0         0
    1.0000    0.1034    0.7241
    1.0000    0.5000    0.3000
         0    0.7000    0.2000
    0.5000    0.5000         0
         0    0.5000    0.5000];
colors = colors*256;
for idx=Start:length(MatFiles)
    filename=MatFiles(idx).name;
    ROIsNb=[];ClusterNb=[];
    %for k = 1 : length(temp)
    for k = 1 : length(temp)
        tempROIsNb=find([temp{k}]<=Numbers(idx+1));
        if tempROIsNb            
            ROIsNb=[ROIsNb temp{k}(tempROIsNb)];
            temp{k}(tempROIsNb)=[];
            ClusterNb=[ClusterNb ; repmat(k,length(tempROIsNb),1)];
        end
    end
    if ROIsNb
        imagename=regexp(filename,'_output_analysis','split');
        %imagename=regexp(imagename,'_output_analysis_matlab2.mat','split');
        imagename=strcat('AVG_',imagename{1},'.tif');
        image=double(imread(imagename));image=image/max(max(image));image=image*90;
        image=uint8(image);
        image2=zeros(size(image(:,:,1)));
        image3=repmat(image,1,1,3);
        ROIs=All_ROIs{idx};       
        ROIsNb=ROIsNb-Numbers(idx);
        ROIs=ROIs(:,ROIsNb);
        for k = 1 : size(ROIs,2)
            image2=zeros(size(image(:,:,1)));
            ROI=full(reshape(ROIs(:,k),size(image,1),size(image,2)));            
            image2=image2+ROI;image2=(image2/max(max(image2)));image2=uint8(image2);
            for j=1:3
                image3(:,:,j)=image3(:,:,j)+image2*colors(ClusterNb(k),j)*0.7;
            end
        end
        %image3(:,:,3)=image;
            name=strcat('_Kmeans_',imagename(4:end));
    imwrite(image3,name,'tif');
    end
    %image3=uint8(image3);

end
clearvars idx i temp tempFileNb fileNb AVG_files filename image counter Numbers image2 image3 k ROI ROIs ROIsNb Start tempROIsNb name imagename tempidx Raster


GoodClustersData=[];
for i=1:length(GoodBetas_select)
    GoodClustersData(i).ZS=ZS(idxKmeans==GoodBetas_select(i),:);
    GoodClustersData(i).Mean=mean(GoodClustersData(i).ZS,1);
    GoodClustersData(i).STD=std(GoodClustersData(i).ZS,1,1);
end
