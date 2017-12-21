MatFiles=dir('*analysis_matlab.mat');i=1;
name=strcat(MatFiles(1).name);
Calcium=load(name, 'DenoisedTraces');
Calcium=Calcium.DenoisedTraces;
MatFiles(1).number=size(Calcium,1);
Spikes=load(name, 'Spikes');
Spikes=Spikes.Spikes;
Noise=load(name, 'Noise');
Noise=Noise.Noise;
Fitness=load(name, 'idx_components');
Fitness=Fitness.idx_components+1;
GoodCalcium=Calcium(Fitness,:);%should be (Fitness,:)
GoodSpikes=Spikes(Fitness,:);
MatFiles(1).GoodNumber=length(Fitness);
MatFiles(1).GC=GoodCalcium;
Cn=load(strrep(name,'analysis_matlab','correlation'), 'Correlation_image');
MatFiles(1).Cn=Cn.Correlation_image;
for i = 2:length(MatFiles)
name=strcat(MatFiles(i).name);
C=load(name, 'DenoisedTraces');
C=C.DenoisedTraces;
%     if i==3
%         C=[C(:,1) C(:,1) C(:,1:58)];
%     end
S=load(name, 'Spikes');
S=S.Spikes;
N=load(name, 'Noise');
N=N.Noise;
F=load(name, 'idx_components');
F=F.idx_components+1;
GC=C(F,:);
GS=S(F,:);
Cn=load(strrep(name,'analysis_matlab','correlation'), 'Correlation_image');
Noise=vertcat(Noise,N);
Calcium=vertcat(Calcium,C);
Spikes=vertcat(Spikes,S);
Fitness=horzcat(Fitness,F);
GoodCalcium=vertcat(GoodCalcium,GC);
GoodSpikes=vertcat(GoodSpikes,GS);
MatFiles(i).number=size(Calcium,1);
MatFiles(i).GoodNumber=MatFiles(i-1).GoodNumber+length(F);
MatFiles(i).GC=GC;
MatFiles(i).Cn=Cn.Correlation_image;
end
clearvars GC C S F N name i GS DF Cn;

ZS=zscore(GoodCalcium,1,2);

Numbers=[1 [MatFiles.GoodNumber]];
counter=1;
idx_Plane=nan(length(GoodCalcium),1);
idx_Position=nan(length(GoodCalcium),1);
name=strcat(MatFiles(1).name);
for i=1:length(MatFiles)	
    name=strcat(MatFiles(i).name);
    [Plane,~]=regexp(name,'Slice(\d+)_','tokens','match');Plane=str2num(Plane{1}{1});
    if strfind(name,'spontaneous')
        Fish=1;
    elseif strfind(name,'loom')
        Fish=2;
    end
    idx_Plane(Numbers(i):Numbers(i+1))=Plane;
    idx_Position(Numbers(i):Numbers(i+1))=Fish;
end
clearvars i Fish Plane name counter

figure;
for i=1:size(Cmap_ZS2,1)
    plot(Cmap_ZS2(i,:));title(num2str(i));pause;
end

idxKmeans_ZS=idxKmeans_ZS2;
GoodBetas=[3 4 6 7 9];
rows=length(GoodBetas);
counter=1;
Fighandle=figure;
set(Fighandle, 'Position', [100, 100, 1300, 900]);
for i=GoodBetas
    idx_temp=find(idxKmeans_ZS==i);
    subplot(rows,4,counter);plot(mean(ZS(idx_temp,:),1));
    subplot(rows,4,counter+1);imagesc(ZS(idx_temp,:),[0 3]);
    subplot(rows,4,counter+2);histogram(idx_Plane(idx_temp));
    subplot(rows,4,counter+3);histogram(idx_Position(idx_temp));
    counter=counter+4;
end

figure;imagesc(ZS(find(idx_Position==2),:),[0 4]);colormap hot

GoodClustersData=[];
for i=1:length(GoodBetas)
    GoodClustersData(i).ZS=ZS(idxKmeans_ZS==GoodBetas(i),:)*100;
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
idxKmeans_ZS_goodmembers=zeros(1,size(GoodCalcium,1));
for i=1:length(GoodBetas)
%GoodClusters_goodmembers(i).Spikes=GoodClustersData(i).Spikes(find(GoodClustersData(i).CorrCoef>=0.5),:);
%GoodClusters_goodmembers(i).ZS=zscore(GoodClustersData(i).DF(find(GoodClustersData(i).CorrCoef>=0.5),:),1,2);
GoodClusters_goodmembers(i).ZS=GoodClustersData(i).ZS(find(GoodClustersData(i).CorrCoef>=Threshold),:);
temp=find(idxKmeans_ZS==GoodBetas(i));
GoodClusters_goodmembers(i).idx=temp(find(GoodClustersData(i).CorrCoef>=0.5));
GoodClusters_goodmembers(i).mean=mean(GoodClusters_goodmembers(i).ZS,1);
GoodClusters_goodmembers(i).STD=std(GoodClusters_goodmembers(i).ZS,1,1);
idx=find(idxKmeans_ZS==GoodBetas(i));
idx=idx(find(GoodClustersData(i).CorrCoef>=Threshold));
idxKmeans_ZS_goodmembers(idx)=GoodBetas(i);
%GoodClusters_goodmembers(i).Fish=idx_Fish(idx);
end

Fighandle=figure;
set(Fighandle, 'Position', [100, 100, 1400, 900]);
counter=1;counter2=1;xplot=floor(sqrt(length(GoodBetas)));yplot=ceil(length(GoodBetas)/xplot);
for i=1:length(GoodBetas)  
    subplot(xplot,yplot,counter);plot(GoodClusters_goodmembers(i).mean,'color',colors(counter2,:)/256);
    counter=counter+1;
    counter2=counter2+1;
end


All_ROIs=[];
ROIs_idx=[];
i=1;
name=strcat(MatFiles(i).name);
Rs=load(name, 'ROIs');
Rs=Rs.ROIs;
F=load(name, 'idx_components');
F=F.idx_components+1;
Rs=Rs(:,F);
All_ROIs{1}=Rs;
ROIs_idx(1)=size(Rs,2);
for i = 2:length(MatFiles)
    name=strcat(MatFiles(i).name);
    Rs=load(name, 'ROIs');
    Rs=Rs.ROIs;
    F=load(name, 'idx_components');
    F=F.idx_components+1;
    Rs=Rs(:,F);
    All_ROIs{i}=Rs;
    ROIs_idx(i)=ROIs_idx(i-1)+length(F);    
end
clearvars GC C S F N name i;

Numbers=[0 [ROIs_idx]];
temp=[];
counter=1;
for i=GoodBetas
    temp{counter}=find(idxKmeans_ZS_goodmembers==i);    
    counter=counter+1;    
end

Start=min(cellfun(@min, temp));Start=find(Numbers<Start,1,'last');
filename=MatFiles(Start).name;

%colors = [1,0,0;0,1,0;0,0,1;1,0.600000000000000,0;1,0,0.7];
colors = distinguishable_colors(length(GoodBetas),[1 1 1; 0 0 0]);
% colors = [0         0    1.0000
%          0    0.5000    1.0000
%     1.0000         0         0
%     1.0000    0.1034    0.7241
%     1.0000    0.5000    0.3000
%          0    0.7000    0.2000
%     0.5000    0.5000         0
%          0    0.5000    0.5000];
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
        imagename=regexp(filename,'_output','split');
        %imagename=regexp(imagename,'_output_analysis_matlab2.mat','split');
        imagename=strcat(imagename{1},'_mean.tif');
        image=double(imread(imagename));image=image/max(max(image));image=image*128;
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
                image3(:,:,j)=image3(:,:,j)+image2*colors(ClusterNb(k),j);
            end
        end
        %image3(:,:,3)=image;
        name=strcat('_Kmeans_',imagename(4:end));
    imwrite(image3,name,'tif');
    end
    %image3=uint8(image3);

end
clearvars idx i temp tempFileNb fileNb AVG_files filename image counter Numbers image2 image3 k ROI ROIs ROIsNb Start tempROIsNb name imagename tempidx Raster

Fighandle=figure;
set(Fighandle, 'Position', [100, 100, 1400, 900]);x = linspace(0.2,size(ZS,2)/5,size(ZS,2));
counter=1;counter2=1;xplot=floor(sqrt(length(GoodBetas)));yplot=ceil(length(GoodBetas)/xplot);
for i=GoodBetas
%     if counter==3
%         counter=counter+1;
%     end
    idx_temp=find(idxKmeans_ZS==i);    
    subplot(xplot,yplot,counter);plot(x,mean(ZS(idx_temp,:),1),'color',colors(counter2,:)/256);%axis([0 131 -1 4]);rectangle('FaceColor','r','Position',[11 -1 10 0.25]);rectangle('FaceColor','r','Position',[51 -1 10 0.25]);rectangle('FaceColor','r','Position',[91 -1 10 0.25]);rectangle('FaceColor','b','Position',[31 -1 10 0.25]);rectangle('FaceColor','b','Position',[71 -1 10 0.25]);rectangle('FaceColor','b','Position',[111 -1 10 0.25]);    
    counter=counter+1;
    counter2=counter2+1;
end

Correla=zeros(1,size(ZS,1));
for i=1:size(ZS,1)
    corr_temp=corrcoef(ZS(i,:), Temp');    
    Correla(i)=corr_temp(1,2);
end
figure;imagesc(ZS(find(Correla>0.3),:),[0 4]);

GCaMP6=[0,5.13796058542217,10.3756715204800,12.2425184714093,8.80831829681196,5.46959264663869,3.42533619066766,1.32732537295463]';
Sound_start=286;interval=23;
Temp_sound=zeros(1,length(OMR));
for i=1:5
    Temp_sound(Sound_start+(i-1)*interval:(Sound_start+(i-1)*interval)+length(GCaMP6)-1)=GCaMP6';
end

Correla_sound=zeros(1,size(ZS,1));
for i=1:size(ZS,1)
    corr_temp=corrcoef(ZS(i,:), Temp_sound');    
    Correla_sound(i)=corr_temp(1,2);
end
figure;imagesc(ZS(find(Correla_sound>0.3),:),[0 4]);

Temp_video=zeros(1,length(OMR));
GCaMP6=[0,1.69644104899772,5.13796058542217,8.27886020152244,10.3756715204800,11.8173714529814,12.2425184714093,10.8571417354877,8.80831829681196,6.91339112244670,5.46959264663869,4.30868766622567,3.42533619066766,2.75378443486879,2.18017250852183,1.72816235135824,1.32732537295463,1.00684435500268,0.730210038304555,0.530242444093118,0.362253250339685,0.227668255288566,0.0869242416152502]';
Video_start=71;interval=35;
for i=1:5
    Temp_video(Video_start+(i-1)*interval:(Video_start+(i-1)*interval)+length(GCaMP6)-1)=GCaMP6';
end

