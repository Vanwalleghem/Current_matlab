MatFiles=dir('*matlab*.mat');
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
GoodCalcium=Calcium(Fitness,:);
GoodSpikes=Spikes(Fitness,:);
GoodNoise=Noise(Fitness,:);
MatFiles(1).GoodNumber=length(Fitness);
%MatFiles(1).GC=GoodCalcium;
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
    GN=N(F,:);
    Noise=vertcat(Noise,N);
    Calcium=vertcat(Calcium,C);
    Spikes=vertcat(Spikes,S);
    Fitness=horzcat(Fitness,F);
    GoodCalcium=vertcat(GoodCalcium,GC);
    GoodSpikes=vertcat(GoodSpikes,GS);
    GoodNoise=vertcat(GoodNoise,GN);
    MatFiles(i).number=size(Calcium,1);
    MatFiles(i).GoodNumber=MatFiles(i-1).GoodNumber+length(F);
    %MatFiles(i).GC=GC;
end
clearvars GC C S F N name i GS;

MatFiles=dir('*matlab*.mat');
name=strcat(MatFiles(1).name);
Noise=load(name, 'Noise');
Noise=Noise.Noise;
Fitness=load(name, 'idx_components');
Fitness=Fitness.idx_components+1;
GoodNoise=Noise(Fitness,:);
for i = 2:length(MatFiles)
    name=strcat(MatFiles(i).name);
    N=load(name, 'Noise');
    N=N.Noise;
    F=load(name, 'idx_components');
    F=F.idx_components+1;
    GN=N(F,:);
    GoodNoise=vertcat(GoodNoise,GN);
end
clearvars GC C S F N name i GS;

ZS=zscore(GoodCalcium,1,2);
ZS2=zscore(GoodCalcium+GoodNoise,1,2);
x = linspace(0.2,size(ZS,2)/5,size(ZS,2));y = linspace(1,size(ZS,1),size(ZS,1));
Fighandle=figure;
set(Fighandle, 'Position', [100, 100, 1500, 1500]);
imagesc(x,y,ZS_rsq(randperm(size(ZS_rsq,1)),:),[-0.5 5]);colormap hot;set(gca,'YTickLabel',[]);
print(Fighandle,strcat('D:\Pictures\processed\Flow\Basic\Figure\FullRaster_rsq.svg'),'-dsvg','-r0');


Fighandle=figure;
set(Fighandle, 'Position', [100, 100, 1500, 1500]);
imagesc(x,y,ZS(randperm(size(ZS,1)),:),[-0.5 5]);colormap hot;set(gca,'YTickLabel',[]);
print(Fighandle,strcat('D:\Pictures\processed\Flow\Basic\Figure\FullRaster_all.svg'),'-dsvg','-r0');

Numbers=[0 [MatFiles.GoodNumber]];
counter=1;
idx_Plane=nan(length(GoodCalcium),1);
idx_Fish=nan(length(GoodCalcium),1);
name=strcat(MatFiles(1).name);
[Plane,~]=regexp(name,'\d+_(\d+)_','tokens','match');Plane=str2num(Plane{1}{1});
% [Plane,~]=regexp(name,'\d\D(\d+)um','tokens','match');Plane=str2num(Plane{1}{1});
[Fish,~]=regexp(name,'(\d+)_\d+_','tokens','match');Fish=str2num(Fish{1}{1});
% [Fish,~]=regexp(name,'(\d)\D\d+um','tokens','match');Fish=str2num(Fish{1}{1});
% idx_Plane(1:Numbers(2))=Plane;
% idx_Fish(1:Numbers(2))=Fish;
for i=1:length(MatFiles)
	%[Fish,~]=regexp(files{i},'(\d+)_','tokens','match');Fish=str2num(Fish{1}{1});
    name=strcat(MatFiles(i).name);
    if findstr(MatFiles(i).name,'2planes')
        [Plane,~]=regexp(name,'f\d-(\d+)um_','tokens','match');Plane=str2num(Plane{1}{1});
        [Fish,~]=regexp(name,'f(\d)-\d+um_','tokens','match');Fish=str2num(Fish{1}{1});
    else
        [Plane,~]=regexp(name,'\d+_(\d+)_','tokens','match');Plane=str2num(Plane{1}{1});
        [Fish,~]=regexp(name,'(\d+)_\d+_','tokens','match');Fish=str2num(Fish{1}{1});
    end
    %[Plane,~]=regexp(name,'\d+_(\d+)_','tokens','match');Plane=str2num(Plane{1}{1});
    %[Fish,~]=regexp(name,'(\d+)_\d+_','tokens','match');Fish=str2num(Fish{1}{1});
   
    idx_Plane(Numbers(i)+1:Numbers(i+1))=Plane;
    idx_Fish(Numbers(i)+1:Numbers(i+1))=Fish;
end
clearvars i Fish Plane name counter

Fighandle=figure;
set(Fighandle, 'Position', [100, 100, 1400, 900]);
counter=1;xplot=floor(sqrt(length(GoodBetas_select)));yplot=ceil(length(GoodBetas_select)/xplot);
x = linspace(0.2,size(Cmap,2)/5,size(Cmap,2));
for i=GoodBetas_select
idx=find(idxKmeans==i);
NumberOfCells=length(idx);
subplot(length(GoodBetas_select),3,counter);plot(x,mean(ZS(idx,:),1));title(num2str(NumberOfCells));axis([0 131 -1 3]);rectangle('FaceColor','r','Position',[10 -1 10 0.2]);rectangle('FaceColor','r','Position',[50 -1 10 0.2]);rectangle('FaceColor','r','Position',[90 -1 10 0.2]);rectangle('FaceColor','b','Position',[30 -1 10 0.2]);rectangle('FaceColor','b','Position',[70 -1 10 0.2]);rectangle('FaceColor','b','Position',[110 -1 10 0.2]);
subplot(length(GoodBetas_select),3,counter+1);imagesc(ZS(idx(randperm(length(idx))),:),[-0.5 5]);colormap hot
subplot(length(GoodBetas_select),3,counter+2);histogram(idx_Plane(idx),[0:20:300]);
counter=counter+3;
end

GoodClustersData=[];
for i=1:length(GoodBetas_select)
    GoodClustersData(i).ZS=ZS(idxKmeans==GoodBetas_select(i),:);
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


GoodClusters_goodmembers=[];Threshold=0.4;
idxKmeans_ZS_goodmembers=zeros(1,size(ZS,1));
for i=1:length(GoodBetas_select)
    %GoodClusters_goodmembers(i).Spikes=GoodClustersData(i).Spikes(find(GoodClustersData(i).CorrCoef>=0.5),:);
    %GoodClusters_goodmembers(i).ZS=zscore(GoodClustersData(i).DF(find(GoodClustersData(i).CorrCoef>=0.5),:),1,2);
    GoodClusters_goodmembers(i).ZS=GoodClustersData(i).ZS(find(GoodClustersData(i).CorrCoef>=Threshold),:);
    temp=find(idxKmeans==GoodBetas_select(i));
    GoodClusters_goodmembers(i).idx=temp(find(GoodClustersData(i).CorrCoef>=Threshold));
    GoodClusters_goodmembers(i).mean=mean(GoodClusters_goodmembers(i).ZS,1);
    GoodClusters_goodmembers(i).STD=std(GoodClusters_goodmembers(i).ZS,1,1);
    idx=find(idxKmeans==GoodBetas_select(i));
    idx=idx(find(GoodClustersData(i).CorrCoef>=Threshold));
    idxKmeans_ZS_goodmembers(idx)=GoodBetas_select(i);
    %GoodClusters_goodmembers(i).Fish=idx_Fish(idx);
end


Fighandle=figure;
set(Fighandle, 'Position', [100, 100, 1400, 900]);x = linspace(0.2,size(Cmap,2)/5,size(Cmap,2));
counter=1;counter2=1;xplot=floor(sqrt(length(GoodBetas_select)));yplot=ceil(length(GoodBetas_select)/xplot);
for i=1:length(GoodBetas_select)  
    if counter==3
        counter=counter+1;
    end
    %subplot(3,3,counter);plot(x,GoodClusters_goodmembers(i).mean,'color',colors(counter2,:)/256);hold on;plot(x,Basic_Clusters(i,:));axis([0 131 -1 4]);rectangle('FaceColor','r','Position',[11 -1 10 0.25]);rectangle('FaceColor','r','Position',[51 -1 10 0.25]);rectangle('FaceColor','r','Position',[91 -1 10 0.25]);rectangle('FaceColor','b','Position',[31 -1 10 0.25]);rectangle('FaceColor','b','Position',[71 -1 10 0.25]);rectangle('FaceColor','b','Position',[111 -1 10 0.25]);
    subplot(3,3,counter);plot(x,GoodClusters_goodmembers(i).mean,'color',colors(counter2,:)/256);axis([0 131 -1 4]);rectangle('FaceColor','r','Position',[11 -1 10 0.25]);rectangle('FaceColor','r','Position',[51 -1 10 0.25]);rectangle('FaceColor','r','Position',[91 -1 10 0.25]);rectangle('FaceColor','b','Position',[31 -1 10 0.25]);rectangle('FaceColor','b','Position',[71 -1 10 0.25]);rectangle('FaceColor','b','Position',[111 -1 10 0.25]);
    counter=counter+1;
    counter2=counter2+1;
end

Fighandle=figure;
set(Fighandle, 'Position', [100, 100, 1400, 900]);
counter=1;xplot=floor(sqrt(length(GoodBetas_select)));yplot=ceil(length(GoodBetas_select)/xplot);
x = linspace(0.2,size(Cmap,2)/5,size(Cmap,2));counter2=1;
for i=GoodBetas_select
idx=GoodClusters_goodmembers(counter2).idx;
NumberOfCells=length(idx);
subplot(length(GoodBetas_select),3,counter);plot(x,GoodClusters_goodmembers(counter2).mean);title(num2str(NumberOfCells));axis([0 131 -1 3]);rectangle('FaceColor','r','Position',[10 -1 10 0.2]);rectangle('FaceColor','r','Position',[50 -1 10 0.2]);rectangle('FaceColor','r','Position',[90 -1 10 0.2]);rectangle('FaceColor','b','Position',[30 -1 10 0.2]);rectangle('FaceColor','b','Position',[70 -1 10 0.2]);rectangle('FaceColor','b','Position',[110 -1 10 0.2]);
subplot(length(GoodBetas_select),3,counter+1);imagesc(ZS(idx(randperm(length(idx))),:),[-0.5 5]);colormap hot
subplot(length(GoodBetas_select),3,counter+2);histogram(idx_Plane(idx),[0:20:300]);
counter=counter+3;
counter2=counter2+1;
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
    temp{counter}=find(idxKmeans_ZS_goodmembers==i);
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
            ROIsNb=[ROIsNb  temp{k}(tempROIsNb)];
            temp{k}(tempROIsNb)=[];
            ClusterNb=[ClusterNb ; repmat(k,length(tempROIsNb),1)];
        end
    end
    if ROIsNb
        imagename=regexp(filename,'_output_analysis','split');
        %imagename=regexp(imagename,'_output_analysis_matlab2.mat','split');
        imagename=strcat('AVG_',imagename{1},'.tif');
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

image_ROI=zeros(size(image));
for i=1:size(ROI,3)
    image_ROI(squeeze(ROI(:,:,i))>0)=i;
end



for idx=Start:length(MatFiles)
    filename=MatFiles(idx).name;
    ROIsNb=[];ClusterNb=[];
    %for k = 1 : length(temp)
    for k = 1 : length(temp)
        tempROIsNb=find([temp{k}]<=Numbers(idx+1));
        if tempROIsNb            
            ROIsNb=[ROIsNb  temp{k}(tempROIsNb)];
            temp{k}(tempROIsNb)=[];
            ClusterNb=[ClusterNb ; repmat(k,length(tempROIsNb),1)];
        end
    end
    if ROIsNb
        imagename=regexp(filename,'_output_analysis','split');
        %imagename=regexp(imagename,'_output_analysis_matlab2.mat','split');
        imagename=strcat('AVG_',imagename{1},'.tif');
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


Fighandle=figure;
set(Fighandle, 'Position', [100, 100, 1400, 900]);
counter=1;xplot=floor(sqrt(length(GoodBetas_select)));yplot=ceil(length(GoodBetas_select)/xplot);
x = linspace(0.2,size(Cmap,2)/5,size(Cmap,2));counter2=1;
h = histogram(idx_Plane,[0:20:300]);
baseline=h.Values;
for i=GoodBetas_select
    idx=GoodClusters_goodmembers(counter2).idx;
    NumberOfCells=length(idx);
    subplot(length(GoodBetas_select),3,counter);plot(x,GoodClusters_goodmembers(counter2).mean);title(num2str(NumberOfCells));axis([0 131 -1 3]);rectangle('FaceColor','r','Position',[10 -1 10 0.2]);rectangle('FaceColor','r','Position',[50 -1 10 0.2]);rectangle('FaceColor','r','Position',[90 -1 10 0.2]);rectangle('FaceColor','b','Position',[30 -1 10 0.2]);rectangle('FaceColor','b','Position',[70 -1 10 0.2]);rectangle('FaceColor','b','Position',[110 -1 10 0.2]);
    subplot(length(GoodBetas_select),3,counter+1);imagesc(ZS(idx(randperm(length(idx))),:),[-0.5 5]);colormap hot
    subplot(length(GoodBetas_select),3,counter+2);h = histogram(idx_Plane(idx),[0:20:300]);bar([0:20:280],(h.Values./baseline)*100);
    
    counter=counter+3;
    counter2=counter2+1;
end

Fighandle=figure;
set(Fighandle, 'Position', [0, 0, 1280, 1024]);
set(findall(Fighandle,'type','text'),'fontSize',12,'fontWeight','bold','FontName','Arial')
counter=1;counter2=1;xplot=length(GoodBetas_select);yplot=4;
back=[55 255 455];
fwd=[155 355 555];
StimLength=100;
x = linspace(0.2,StimLength/5,StimLength);
for i=GoodBetas_select
    subplot(xplot,yplot,counter2);
    imagesc(GoodClusters_goodmembers(counter).ZS,[-0.5 4]);colormap hot
    tempPlot=GoodClusters_goodmembers(counter).mean;
    BackPlot=zeros(3,StimLength);
    FwdPlot=zeros(3,StimLength);
    for j=1:3
        BackPlot(j,:)=tempPlot(back(j):back(j)+99);
        FwdPlot(j,:)=tempPlot(fwd(j):fwd(j)+99);
    end
    temp=mean(BackPlot,1);std_temp=std(BackPlot,1,1);
    subplot(xplot,yplot,counter2+1);
    H=shadedErrorBar(x, temp, std_temp);axis([0 20 -1 3]);
    H.mainLine.Color=colors(counter,:)/256;
    H.patch.FaceColor=colors(counter,:)/256;
    H.edge(1).Color=colors(counter,:)/256;
    H.edge(2).Color=colors(counter,:)/256;
    temp=mean(FwdPlot,1);std_temp=std(FwdPlot,1,1);
    subplot(xplot,yplot,counter2+2);
    H=shadedErrorBar(x, temp, std_temp);axis([0 20 -1 3]);
    H.mainLine.Color=colors(counter,:)/256;
    H.patch.FaceColor=colors(counter,:)/256;
    H.edge(1).Color=colors(counter,:)/256;
    H.edge(2).Color=colors(counter,:)/256;    
    subplot(xplot,yplot,counter2+3);
    RespBWD=max(BackPlot,[],2);RespBWD=RespBWD+abs(min(RespBWD));
    RespFWD=max(FwdPlot,[],2);RespFWD=RespFWD+abs(min(RespFWD));
    %bar(mean((RespFWD-RespBWD)./(RespFWD+RespBWD)),'FaceColor',colors(counter,:)/256);hold on; ylim([-1 1]);
    scatter(zeros(1,13)+1,DSI(:,counter,1),[],colors(counter,:)/256);hold on;ylim([-1.1 1.1]);set(gca,'xtick',[]);set(gca,'xcolor','none')
    errorbar(mean((RespFWD-RespBWD)./(RespFWD+RespBWD)),nanstd(DSI(:,counter,1)),'.','LineWidth',2,'MarkerEdgeColor',colors(counter,:)/256,'MarkerFaceColor',colors(counter,:)/256,'Color',colors(counter,:)/256);view([90 -90]);hold off    
    counter=counter+1;
    counter2=counter2+yplot;
end

back=[55 255 455];
fwd=[155 355 555];
StimLength=100;
DS_raw=NaN(length(unique(idx_Fish)),3,2,length(GoodBetas_select));counter2=1;
for j=GoodBetas_select
    idx_temp=find(idxKmeans_ZS_goodmembers==j);
    idx_fish_temp=idx_Fish(idx_temp);
    counter=1;
    for fish=(unique(idx_Fish)')
        if find(idx_fish_temp==fish)
            tempPlot=mean(ZS(idx_temp(find(idx_fish_temp==fish)),:),1);
            BackPlot=zeros(3,StimLength);
            FwdPlot=zeros(3,StimLength);
            for i=1:3
                BackPlot(i,:)=tempPlot(back(i):back(i)+99);
                FwdPlot(i,:)=tempPlot(fwd(i):fwd(i)+99);
            end
            RespBWD=max(BackPlot,[],2);RespBWD=RespBWD+abs(min(RespBWD));
            RespFWD=max(FwdPlot,[],2);RespFWD=RespFWD+abs(min(RespFWD));
            DS_raw(counter,:,1,counter2)=RespBWD;
            DS_raw(counter,:,2,counter2)=RespFWD;
        end
        counter=counter+1;
    end
    counter2=counter2+1;
end

DSI=zeros(length(unique(idx_Fish)),length(GoodBetas_select),2);
for i=1:size(DSI,1)
    for j=1:size(DSI,2)
        RespFWD=DS_raw(i,:,2,j);
        RespBWD=DS_raw(i,:,1,j);
        DSI(i,j,1)=nanmean((RespFWD-RespBWD)./(RespFWD+RespBWD));
        DSI(i,j,2)=nanstd((RespFWD-RespBWD)./(RespFWD+RespBWD));
    end
end

back=[55 255 455];back=back-5;
fwd=[155 355 555];fwd=fwd-5;
StimLength=100;
DS_raw_inhib=NaN(length(unique(idx_Fish)),3,2,length(3));counter2=1;
for j=[1 3 4]
    idx_temp=GoodInhib_goodmembers(j).idx;
    idx_fish_temp=idx_Fish(idx_temp);
    counter=1;
    for fish=(unique(idx_Fish)')
        if find(idx_fish_temp==fish)
            tempPlot=mean(ZS(idx_temp(find(idx_fish_temp==fish)),:),1);
            BackPlot=zeros(3,StimLength);
            FwdPlot=zeros(3,StimLength);
            for i=1:3
                BackPlot(i,:)=tempPlot(back(i):back(i)+99);
                FwdPlot(i,:)=tempPlot(fwd(i):fwd(i)+99);
            end
            RespBWD=min(BackPlot,[],2);RespBWD=RespBWD-abs(max(RespBWD));
            RespFWD=min(FwdPlot,[],2);RespFWD=RespFWD-abs(max(RespFWD));
            DS_raw_inhib(counter,:,1,counter2)=RespBWD;
            DS_raw_inhib(counter,:,2,counter2)=RespFWD;
        end
        counter=counter+1;
    end
    counter2=counter2+1;
end

DSI_inhib=zeros(length(unique(idx_Fish)),3,2);
for i=1:size(DSI_inhib,1)
    for j=1:size(DSI_inhib,2)
        RespFWD=DS_raw_inhib(i,:,2,j);
        RespBWD=DS_raw_inhib(i,:,1,j);
        DSI_inhib(i,j,1)=nanmean((RespFWD-RespBWD)./(RespFWD+RespBWD));
        DSI_inhib(i,j,2)=nanstd((RespFWD-RespBWD)./(RespFWD+RespBWD));
    end
end

back=[55 255 455];
fwd=[155 355 555];
StimLength=100;
TimeToPeak_raw_inhib=NaN(length(unique(idx_Fish)),3,2,length(GoodBetas_select));counter2=1;
for j=[1 3 4]
    idx_temp=GoodInhib_goodmembers(j).idx;
    idx_fish_temp=idx_Fish(idx_temp);
    counter=1;
    for fish=(unique(idx_Fish)')
        if find(idx_fish_temp==fish)
            tempPlot=mean(ZS(idx_temp(find(idx_fish_temp==fish)),:),1);
            BackPlot=zeros(3,StimLength);
            FwdPlot=zeros(3,StimLength);
            for i=1:3
                BackPlot(i,:)=tempPlot(back(i):back(i)+99);
                FwdPlot(i,:)=tempPlot(fwd(i):fwd(i)+99);
            end
            [~,RespBWD]=min(BackPlot,[],2);
            [~,RespFWD]=min(FwdPlot,[],2);
            TimeToPeak_raw_inhib(counter,:,1,counter2)=RespBWD;
            TimeToPeak_raw_inhib(counter,:,2,counter2)=RespFWD;
        end
        counter=counter+1;
    end
    counter2=counter2+1;
end

TimeToPeak_AVG_inhib=zeros(length(unique(idx_Fish)),3,2);
for i=1:size(TimeToPeak_AVG_inhib,1)
    for j=1:size(TimeToPeak_AVG_inhib,2)
        RespFWD=TimeToPeak_raw_inhib(i,:,2,j);
        RespBWD=TimeToPeak_raw_inhib(i,:,1,j);
        TimeToPeak_AVG_inhib(i,j,1)=nanmean(RespFWD);
        TimeToPeak_AVG_inhib(i,j,2)=nanmean(RespBWD);
    end
end

TimeToPeak_AVG_inhib=TimeToPeak_AVG_inhib/5;

%Spikes are not worth it
% Fighandle=figure;
% set(Fighandle, 'Position', [100, 100, 1400, 900]);
% counter=1;counter2=1;xplot=floor(sqrt(length(GoodBetas_select)));yplot=ceil(length(GoodBetas_select)/xplot);
% back=[55 255 455];
% fwd=[155 355 555];
% StimLength=100;
% x = linspace(0.2,StimLength/5,StimLength);
% for i=GoodBetas_select
%     idx_temp=find(idxKmeans_ZS_goodmembers==i);
%     tempPlot=mean(GoodSpikes(idx_temp,:),1);
%     BackPlot=zeros(3,StimLength);
%     FwdPlot=zeros(3,StimLength);
%     for i=1:3
%         BackPlot(i,:)=tempPlot(back(i):back(i)+99);
%         FwdPlot(i,:)=tempPlot(fwd(i):fwd(i)+99);
%     end
%     temp=mean(BackPlot,1);std_temp=std(BackPlot,1,1);
%     subplot(length(GoodBetas_select),2,counter2);
%     H=shadedErrorBar(x, temp, std_temp);%axis([0 20 -1 3]);
%     H.mainLine.Color=colors(counter,:)/256;
%     H.patch.FaceColor=colors(counter,:)/256;
%     H.edge(1).Color=colors(counter,:)/256;
%     H.edge(2).Color=colors(counter,:)/256;
%     temp=mean(FwdPlot,1);std_temp=std(FwdPlot,1,1);
%     subplot(length(GoodBetas_select),2,counter2+1);
%     H=shadedErrorBar(x, temp, std_temp);%axis([0 20 -1 3]);
%     H.mainLine.Color=colors(counter,:)/256;
%     H.patch.FaceColor=colors(counter,:)/256;
%     H.edge(1).Color=colors(counter,:)/256;
%     H.edge(2).Color=colors(counter,:)/256;    
%     counter=counter+1;
%     counter2=counter2+2;
% end
% 
% for i=1:8;
% figure;subplot(1,2,1);plot(mean(ZS(find(idxKmeans_ZS_goodmembers==GoodBetas_select(i)),:),1))
% subplot(1,2,2);plot(mean(GoodSpikes(find(idxKmeans_ZS_goodmembers==GoodBetas_select(i)),:),1))
% end

for i=1:length(GoodBetas_select)
temp=find(idxKmeans_ZS_goodmembers==GoodBetas_select(i));
GoodClusters_goodmembers(i).plane=idx_Plane(temp);
GoodClusters_goodmembers(i).Fish=idx_Fish(temp);
end

% Fighandle=figure;
% set(Fighandle, 'Position', [0, 0, 1000, 1000]);
% set(findall(Fighandle,'type','text'),'fontSize',12,'fontWeight','bold','FontName','Arial')
% xplot=floor(sqrt(length(GoodBetas_select)));yplot=ceil(length(GoodBetas_select)/xplot);
% for i=1:length(GoodBetas_select)
%     subplot(xplot,yplot,i);histogram(GoodClusters_goodmembers(i).Fish);
% end

Planes_perFish=cell(length(unique(idx_Fish)),length(GoodBetas_select));
FishList=unique(idx_Fish);
for j=1:size(Planes_perFish,2)
    planes=GoodClusters_goodmembers(j).plane;
    fish_idx=GoodClusters_goodmembers(j).Fish;
    for i=1:size(Planes_perFish,1)        
        Planes_perFish{i,j}=planes(find(fish_idx==FishList(i)));
    end
end
clearvars i j planes fish_idx FishList

Fighandle=figure;
set(Fighandle, 'Position', [0, 0, 1000, 1000]);
set(findall(Fighandle,'type','text'),'fontSize',12,'fontWeight','bold','FontName','Arial')
xplot=floor(sqrt(length(GoodBetas_select)));yplot=ceil(length(GoodBetas_select)/xplot);
yplot=1;xplot=8;
for j=1:size(Planes_perFish,2)
    subplot(xplot,yplot,j);
    for i=1:size(Planes_perFish,1)   
        histogram(Planes_perFish{i, j},[0:20:300]);hold on;
    end
    hold off;
end

Planes_KStest=zeros(3,length(GoodBetas_select),length(GoodBetas_select));
for i=1:length(GoodBetas_select)
    for j=1:length(GoodBetas_select)
        [Planes_KStest(1,i,j),Planes_KStest(2,i,j),Planes_KStest(3,i,j)] = kstest2(GoodClusters_goodmembers(i).plane,GoodClusters_goodmembers(j).plane);        
    end
end
clearvars i j
Planes_KS_pvalues=triu(squeeze(Planes_KStest(2,:,:)));Planes_KS_pvalues(Planes_KS_pvalues==0)=nan;
[pVals_KS KS_idx]=sort(Planes_KS_pvalues(:));
Multiple_comparison=((length(GoodBetas_select)-1)^2+(length(GoodBetas_select)-1))/2;
pval=0.05;
for i=1:length(pVals_KS)
    if pVals_KS(i) >= pval/(Multiple_comparison-(i-1))
        break
    end
end
pVals_KS(i:end)=nan;
KS_idx(isfinite(pVals_KS))=[];
Planes_KS_pvalues(KS_idx)=1;

load('D:\Pictures\processed\Flow\Basic\transformation_crop_to_full.mat')
Twoplanes_transfo=load('D:\Pictures\processed\Flow\Basic\transformation_crop_to_full_2planes.mat')
Errored_ROI={};
progressbar(0,0,0);
for fish_nb=1:length(Fish_list)    
    IndexC=find(idx_Fish==Fish_list(fish_nb));
    IndexC=ismember(ROIs_idx,IndexC);
    MatFiles_fish = find(IndexC>0);
    if MatFiles_fish(1)==1
        numbersForROIs=[1 [MatFiles(MatFiles_fish).GoodNumber]];
    else
        numbersForROIs=[MatFiles(MatFiles_fish(1)-1).GoodNumber [MatFiles(MatFiles_fish).GoodNumber]];
    end
    Centroids=[];
    counter=1;
    for plane=1:length(MatFiles_fish)
        filename=MatFiles(MatFiles_fish(plane)).name;
        if findstr(filename,'2planes')
            [slice,~]=regexp(filename,'f\d-(\d+)um_','tokens','match');slice=str2num(slice{1}{1});
            idx_name=strcat(num2str(Fish_list(fish_nb)),'-',num2str(slice));
            Transfo_idx = ~cellfun('isempty',strfind(cellstr(Twoplanes_transfo.FileName),idx_name));Transfo_idx=find(Transfo_idx>0);
        else
            [slice,~]=regexp(filename,'\d+_(\d+)_','tokens','match');slice=str2num(slice{1}{1});
            idx_name=strcat(num2str(Fish_list(fish_nb)),'_',num2str(slice));
            Transfo_idx = ~cellfun('isempty',strfind(cellstr(FileName),strcat(idx_name,'.tif')));Transfo_idx=find(Transfo_idx>0);
        end
        ROI=All_ROIs{MatFiles_fish(plane)};
        imagename=regexp(filename,'_output_analysis','split');
        imagename=strcat('AVG_',imagename{1},'.tif');
        image=double(imread(imagename));image=image/max(max(image));image=image*128;
        ROI=reshape(full(ROI),size(image,1),size(image,2),size(ROI,2));        
        if isempty(Transfo_idx)
            FileName
            break
        end
        if findstr(filename,'2planes')
            Transfo_tmp=Twoplanes_transfo.TransfoMatrix(Transfo_idx,:);
             if Transfo_tmp(1)<-5
                Transfo_tmp(1)=830+Transfo_tmp(1);
            end
            if Transfo_tmp(2)<-5
                Transfo_tmp(2)=1120+Transfo_tmp(2);
            end
        else
            Transfo_tmp=TransfoMatrix(Transfo_idx,:);
            if Transfo_tmp(1)<-5
                Transfo_tmp(1)=1080+Transfo_tmp(1);
            end
            if Transfo_tmp(2)<-5
                Transfo_tmp(2)=1280+Transfo_tmp(2);
            end
        end
        for roi_nb=1:size(ROI,3)
            test=ROI(:,:,roi_nb);
            [M I]=max(test(:));
            [I_row I_col]=ind2sub(size(test),I);            
            Centroids(counter,5)=counter;
            temp=[I_col I_row];
            Centroids(counter,1:2)=temp+Transfo_tmp([2 1]);            
            Centroids(counter,3)=(slice/20)+1;
            progressbar([],[],roi_nb/size(ROI,3));
            counter=counter+1;
        end
        if findstr(filename,'2planes')
            test=Centroids(:,1:2);
            test(:,1)=Centroids(:,2);
            test(:,2)=1120-Centroids(:,1);
            Centroids(:,1:2)=test;
        end
        progressbar([],fish_nb/length(MatFiles_fish));
    end
    if iscell(Fish_list)
        image_name=strcat('_ROIsFish_',Fish_list{fish_nb},'.csv');
    else
        image_name=strcat('_ROIsFish_',num2str(Fish_list(fish_nb)),'.csv');
    end
    csvwrite(image_name,Centroids);
    progressbar(fish_nb/length(Fish_list));
end
clearvars i j k ROI Centroids M I I_row I_col test temp image_name filename fish fish_nb


idx_name=strcat(num2str(Fish_list(fish_nb)),'_',num2str(slice));
idx = all(ismember(FileName,idx_name),2);find(idx>0)

strfind(FileName,idx_name)


%Make individual figures for FWD and BWD in each cluster
Fighandle=figure;
set(Fighandle, 'Position', [100, 100, 1200, 400]);
set(findall(Fighandle,'type','text'),'fontSize',12,'fontWeight','bold','FontName','Arial')
counter=1;counter2=1;xplot=1;yplot=2;
back=[55 255 455];
fwd=[155 355 555];
StimLength=100;
x = linspace(0.2,StimLength/5,StimLength);
ha = tight_subplot(xplot,yplot,[.01 .01],[.01 .01],[.01 .01]);
for i=GoodBetas_select    
    tempPlot=GoodClusters_goodmembers(counter).mean;
    BackPlot=zeros(3,StimLength);
    FwdPlot=zeros(3,StimLength);
    for j=1:3
        BackPlot(j,:)=tempPlot(back(j):back(j)+99);
        FwdPlot(j,:)=tempPlot(fwd(j):fwd(j)+99);
    end
    temp=mean(BackPlot,1);std_temp=std(BackPlot,1,1);
    %subplot(xplot,yplot,1);
    axes(ha(1));
    H=shadedErrorBar(x, temp, std_temp);axis([0 20 -1 3]);
    H.mainLine.Color=colors(counter,:)/256;
    H.patch.FaceColor=colors(counter,:)/256;
    H.edge(1).Color=colors(counter,:)/256;
    H.edge(2).Color=colors(counter,:)/256;
    H.mainLine.LineWidth=3;
    set(gca,'XTickLabel',[]);set(gca,'YTickLabel',[]);
    temp=mean(FwdPlot,1);std_temp=std(FwdPlot,1,1);
    %subplot(xplot,yplot,2);
    axes(ha(2));
    H=shadedErrorBar(x, temp, std_temp);axis([0 20 -1 3]);
    H.mainLine.Color=colors(counter,:)/256;
    H.patch.FaceColor=colors(counter,:)/256;
    H.edge(1).Color=colors(counter,:)/256;
    H.edge(2).Color=colors(counter,:)/256;    
    H.mainLine.LineWidth=3;
    set(gca,'XTickLabel',[]);set(gca,'YTickLabel',[]);
%     subplot(xplot,yplot,3);
%     RespBWD=max(BackPlot,[],2);RespBWD=RespBWD+abs(min(RespBWD));
%     RespFWD=max(FwdPlot,[],2);RespFWD=RespFWD+abs(min(RespFWD));
%     %bar(mean((RespFWD-RespBWD)./(RespFWD+RespBWD)),'FaceColor',colors(counter,:)/256);hold on; ylim([-1 1]);
%     scatter(zeros(1,13)+1,DSI(:,counter,1),[],colors(counter,:)/256);hold on;ylim([-1.1 1.1]);set(gca,'xtick',[]);set(gca,'xcolor','none')
%     errorbar(mean((RespFWD-RespBWD)./(RespFWD+RespBWD)),nanstd(DSI(:,counter,1)),'.','LineWidth',2,'MarkerEdgeColor',colors(counter,:)/256,'MarkerFaceColor',colors(counter,:)/256,'Color',colors(counter,:)/256);view([90 -90]);hold off    
%     set(gca,'XTickLabel',[]);set(gca,'YTickLabel',[]);
    print(Fighandle,strcat('D:\Pictures\processed\Flow\Basic\Figure\Cluster_basic_',num2str(counter),'.svg'),'-dsvg','-r0');    
    counter=counter+1;
end

%Make figure of full length + STD + rasterhisto
Fighandle=figure;
set(Fighandle, 'Position', [100, 100, 800, 3000]);
set(findall(Fighandle,'type','text'),'fontSize',12,'fontWeight','bold','FontName','Arial')
counter=1;counter2=1;xplot=8;yplot=2;
StimLength=655;
x = linspace(0.2,StimLength/5,StimLength);
ha = tight_subplot(xplot,yplot,[.01 .01],[.01 .01],[.01 .01]);
for i=GoodBetas_select    
    temp=GoodClusters_goodmembers(counter).mean;
    std_temp=GoodClusters_goodmembers(counter).STD;        
    %subplot(xplot,yplot,(2*counter)-1);
    axes(ha((2*counter)-1));
    H=shadedErrorBar(x, temp, std_temp);axis([0 130 -1 4]);
    H.mainLine.Color=colors(counter,:)/256;
    H.patch.FaceColor=colors(counter,:)/256;
    H.edge(1).Color=colors(counter,:)/256;
    H.edge(2).Color=colors(counter,:)/256;
    H.mainLine.LineWidth=3;    
    set(gca,'XTickLabel',[]);set(gca,'YTickLabel',[]);  
    %subplot(xplot,yplot,(2*counter));
    axes(ha((2*counter)));
    imagesc(GoodClusters_goodmembers(counter).ZS(randperm(size(GoodClusters_goodmembers(counter).ZS,1)),:),[-1 3]);colormap hot
    set(gca,'XTickLabel',[]);set(gca,'YTickLabel',[]);    
    counter=counter+1;
end
print(Fighandle,strcat('D:\Pictures\processed\Flow\Basic\Figure\FullLength_basic_clusters.svg'),'-dsvg','-r0');    

%Make individual figures of full length + STD + rasterhisto
Fighandle=figure;
set(Fighandle, 'Position', [100, 100, 1600, 400]);
set(findall(Fighandle,'type','text'),'fontSize',12,'fontWeight','bold','FontName','Arial')
counter=1;counter2=1;xplot=1;yplot=2;
StimLength=655;
x = linspace(0.2,StimLength/5,StimLength);
for i=GoodBetas_select    
    ha = tight_subplot(xplot,yplot,[.01 .01],[.01 .01],[.01 .01]);
    temp=GoodClusters_goodmembers(counter).mean;
    std_temp=GoodClusters_goodmembers(counter).STD;        
    %subplot(xplot,yplot,(2*counter)-1);
    axes(ha(1));
    H=shadedErrorBar(x, temp, std_temp);axis([0 130 -2 5]);
    H.mainLine.Color=colors(counter,:)/256;
    H.patch.FaceColor=colors(counter,:)/256;
    H.edge(1).Color=colors(counter,:)/256;
    H.edge(2).Color=colors(counter,:)/256;
    H.mainLine.LineWidth=3;    
    set(gca,'XTickLabel',[]);set(gca,'YTickLabel',[]);  
    %subplot(xplot,yplot,(2*counter));
    axes(ha(2));
    imagesc(GoodClusters_goodmembers(counter).ZS(randperm(size(GoodClusters_goodmembers(counter).ZS,1)),:),[-1 3]);colormap hot
    set(gca,'XTickLabel',[]);set(gca,'YTickLabel',[]);    
    
    print(Fighandle,strcat('D:\Pictures\processed\Flow\Basic\Figure\FullLength_basic_clusters',num2str(counter),'.svg'),'-dsvg','-r0');    
    counter=counter+1;
end

%Make individual figures of full length + STD + rectangles
Fighandle=figure;
set(Fighandle, 'Position', [100, 100, 1600, 800]);
set(findall(Fighandle,'type','text'),'fontSize',12,'fontWeight','bold','FontName','Arial')
counter=1;counter2=1;xplot=1;yplot=1;
StimLength=655;
x = linspace(0.2,StimLength/5,StimLength);
for i=GoodBetas_select    
    ha = tight_subplot(xplot,yplot,[.01 .01],[.01 .01],[.01 .01]);
    temp=GoodClusters_goodmembers(counter).mean;
    std_temp=GoodClusters_goodmembers(counter).STD;        
    %subplot(xplot,yplot,(2*counter)-1);
    
    H=shadedErrorBar(x, temp, std_temp);axis([0 130 -2 5]);
    H.mainLine.Color=colors(counter,:)/256;
    H.patch.FaceColor=colors(counter,:)/256;
    H.edge(1).Color=colors(counter,:)/256;
    H.edge(2).Color=colors(counter,:)/256;
    H.mainLine.LineWidth=3;    
    set(gca,'XTickLabel',[]);set(gca,'YTickLabel',[]);
   %subplot(xplot,yplot,(2*counter));
    rectangle('FaceColor','m','Position',[11 -2 10 0.25]);rectangle('FaceColor','m','Position',[51 -2 10 0.25]);rectangle('FaceColor','m','Position',[91 -2 10 0.25]);rectangle('FaceColor','g','Position',[31 -2 10 0.25]);rectangle('FaceColor','g','Position',[71 -2 10 0.25]);rectangle('FaceColor','g','Position',[111 -2 10 0.25]);
    
    print(Fighandle,strcat('D:\Pictures\processed\Flow\Basic\Figure\FullPlot_basic_clusters',num2str(counter),'.svg'),'-dsvg','-r0');    
    
    
    imagesc(GoodClusters_goodmembers(counter).ZS(randperm(size(GoodClusters_goodmembers(counter).ZS,1)),:),[-1 4]);colormap hot
    set(gca,'XTickLabel',[]);set(gca,'YTickLabel',[]);    
    
    print(Fighandle,strcat('D:\Pictures\processed\Flow\Basic\Figure\FullRaster_basic_clusters',num2str(counter),'.svg'),'-dsvg','-r0');    
    counter=counter+1;
end
colorbar;


parfor i=1:size(ZS,1)
    mdl=fitlm(flow',ZS(i,:));
    model_DF_Thr5(i).coef=mdl.Coefficients;
    model_DF_Thr5(i).rsquared=mdl.Rsquared.Adjusted;
end

coefficients={};
for idx=1:length(model_DF_Thr5)
    coef=[model_DF_Thr5(idx).coef];
    temp=coef.Properties.RowNames;temp=regexp(temp,'x(\d+)','tokens');
    if ~isempty(temp)
        %temp=[temp{:}];temp=[temp{:}];temp=[temp{:}];%temp=str2num(temp);
        for coef_idx=2:height(coef)
            if coef.pValue(coef_idx)<0.05
                coefficients{idx,str2num(temp{coef_idx}{1}{1})}=coef.Estimate(coef_idx);
            end
        end
    end
end
idxempty=cellfun('isempty',coefficients);
coefficients(idxempty)={0};
clearvars idxempty idx coef_idx coef temp
coefficients=cell2mat(coefficients);

CSV_Files=dir('_2WarpedLong*.csv');
ROIs=struct();
for i=1:length(CSV_Files);
    temp=csvread(CSV_Files(i).name,0);   
    Fishname=regexp(CSV_Files(i).name,'_2WarpedLong(\d+).csv','tokens');Fishname=Fishname{1}{1};    
    ROIs(i).name=Fishname;    
    ROIs(i).coord=temp(:,1:3);
    ROIs(i).idx=temp(:,5);
    size(temp,1)==sum(idx_Fish==str2num(Fishname))
end
clearvars i temp CSV_Files Fishname

i=1;ROI_pool=ROIs(i).coord;
for i=2:length(ROIs);
    ROI_pool=[ROI_pool; ROIs(i).coord];
end

Sort_ROIs=[];temp_nb=0;
MatFiles_names={MatFiles.name};
for fish_nb=1:length(Fish_list)
    if fish_nb <3
        fish_name=strcat('f',num2str(Fish_list(fish_nb)),'-');
        IndexC=strfind(MatFiles_names,fish_name);
    else
        fish_name=strcat(num2str(Fish_list(fish_nb)),'_');
        IndexC=regexp(MatFiles_names,strcat('^',fish_name));
    end            
    MatFiles_fish = find(not(cellfun('isempty', IndexC)));
    for file_nb=1:length(MatFiles_fish)
        if MatFiles_fish(file_nb)==1
            numbersForROIs=[1 MatFiles(1).GoodNumber];
        else
            numbersForROIs=[MatFiles(MatFiles_fish(file_nb)-1).GoodNumber+1 MatFiles(MatFiles_fish(file_nb)).GoodNumber];
        end
        if ismember(numbersForROIs,Sort_ROIs)
            fish_name
            break
        end
        Sort_ROIs=[Sort_ROIs numbersForROIs(1):1:numbersForROIs(2)];        
    end    
    if ~length(Sort_ROIs)-temp_nb==sum(idx_Fish==(Fish_list(fish_nb)))
        fish_name
        break
    end
    temp_nb=length(Sort_ROIs);
end
clearvars slice fish_nb roi_nb ROI Centroids IndexC file_nb

ROI_fish(Sort_ROIs,:)=ROI_pool;
clearvars ROI_pool
% 
% for i=1:length(GoodBetas_select)
%     idx_temp=GoodClusters_goodmembers(i).idx;    
%     CSV_temp=ROI_fish(idx_temp,:);
%     CSV_temp(:,3)=CSV_temp(:,3);
%     CSV_temp(:,4)=1;
%     filename=strcat('__Coords_clust_Basic',num2str(i),'.csv');
%     csvwrite(filename,CSV_temp);
% end
% 
for i=1:length(ROIs)
    CSV_temp=ROIs(i).coord;    
    CSV_temp(:,4)=1;
    filename=strcat('__WB_F',ROIs(i).name,'.csv');    
    csvwrite(filename,CSV_temp);
end

% % 
% % CSV_temp=ROIs(2).coord;
% CSV_temp(:,3)=CSV_temp(:,3);
% CSV_temp(:,4)=1;
% filename=strcat('__WB_F13','.csv');
% csvwrite(filename,CSV_temp);
% 
% CSV_temp=ROIs(3).coord;
% CSV_temp(:,3)=CSV_temp(:,3);
% CSV_temp(:,4)=1;
% filename=strcat('__WB_F14','.csv');
% csvwrite(filename,CSV_temp);

idx_inhib=find(min(coefficients,[],2)<-1);
figure;imagesc(ZS(idx_inhib,:),[-2 2]);colormap hot;options = statset('UseParallel',1);
[idxKmeans_inhib Cmap_inhib]=kmeans(ZS(idx_inhib,:),4,'Options',options,'Distance','cityblock','Replicates',3,'MaxIter',1000,'Display','final');
[~,GoodBetas_inhib]=Test_Regress(Cmap_inhib,flow,idxKmeans_inhib,0.4);

GoodInhib_goodmembers=[];Threshold=0.4;
for i=1:length(GoodBetas_inhib)    
    ZS_temp=ZS(idx_inhib(find(idxKmeans_inhib==GoodBetas_inhib(i))),:);
    corr_temp=zeros(1,length(ZS_temp));
    for jj=1:size(ZS_temp,1)
        temp_corr=corrcoef(Cmap_inhib(GoodBetas_inhib(i),:), ZS_temp(jj,:));
        corr_temp(jj)=temp_corr(1,2);
    end    
    GoodInhib_goodmembers(i).ZS=ZS_temp(find(corr_temp>Threshold),:);
    idx_temp=idx_inhib(find(idxKmeans_inhib==GoodBetas_inhib(i)));
    GoodInhib_goodmembers(i).idx=idx_temp(find(corr_temp>Threshold));
    GoodInhib_goodmembers(i).mean=mean(GoodInhib_goodmembers(i).ZS,1);
    GoodInhib_goodmembers(i).STD=std(GoodInhib_goodmembers(i).ZS,1,1);       
end


for i=1:length(RegionList_select)
Fighandle=figure;set(Fighandle, 'Position', [0, 0, 1000, 500]);
regionName=RegionList_select{i};
idx_temp=PerBrainRegions.(regionName).idx;
idx_temp=idxKmeans_ZS_goodinhib(idx_temp);
idx_temp(idx_temp==0)=nan;
h=histogram(idx_temp,'Normalization','probability');
h.FaceColor = [0.5 0.5 0.5];
h.EdgeColor = 'k';
set(gca,'XTick',[1:1:length(RegionList_select)]);
set(gca,'YTick',[0.1:0.1:0.5]);set(gca,'YTickLabel',[10:10:50]);
box off
ax1 = gca;
ax2 = axes('Position', get(ax1, 'Position'), 'FontSize', 10,...
'Color','None','XColor','k','YColor','k', 'LineWidth', 1,...
'XAxisLocation','top', 'XTick', [],...
'YAxisLocation','right', 'YTick', []);
linkaxes([ax1, ax2]);
print(Fighandle,strcat('D:\Pictures\processed\Flow\Basic\Figure\InhibClustersDistrib_',regionName,'.png'),'-dpng','-r0');
close;
end


colors_inhib=colors([7 1 2 4],:);
%Make individual figures of full length + STD + rectangles
Fighandle=figure;
set(Fighandle, 'Position', [100, 100, 1600, 800]);
set(findall(Fighandle,'type','text'),'fontSize',12,'fontWeight','bold','FontName','Arial')
counter=1;counter2=1;xplot=1;yplot=1;
StimLength=655;
x = linspace(0.2,StimLength/5,StimLength);
for i=GoodBetas_inhib
ha = tight_subplot(xplot,yplot,[.01 .01],[.01 .01],[.01 .01]);
temp=GoodInhib_goodmembers(counter).mean;
std_temp=GoodInhib_goodmembers(counter).STD;
%subplot(xplot,yplot,(2*counter)-1);
H=shadedErrorBar(x, temp, std_temp);axis([0 130 -2 5]);
H.mainLine.Color=colors_inhib(counter,:)/256;
H.patch.FaceColor=colors_inhib(counter,:)/256;
H.edge(1).Color=colors_inhib(counter,:)/256;
H.edge(2).Color=colors_inhib(counter,:)/256;
H.mainLine.LineWidth=3;
set(gca,'XTickLabel',[]);set(gca,'YTickLabel',[]);
%subplot(xplot,yplot,(2*counter));
rectangle('FaceColor','m','Position',[11 -2 10 0.25]);rectangle('FaceColor','m','Position',[51 -2 10 0.25]);rectangle('FaceColor','m','Position',[91 -2 10 0.25]);rectangle('FaceColor','g','Position',[31 -2 10 0.25]);rectangle('FaceColor','g','Position',[71 -2 10 0.25]);rectangle('FaceColor','g','Position',[111 -2 10 0.25]);
print(Fighandle,strcat('D:\Pictures\processed\Flow\Basic\Figure\FullPlot_inhib_clusters',num2str(counter),'.svg'),'-dsvg','-r0');
imagesc(GoodInhib_goodmembers(counter).ZS(randperm(size(GoodInhib_goodmembers(counter).ZS,1)),:),[-1 4]);colormap hot
set(gca,'XTickLabel',[]);set(gca,'YTickLabel',[]);
print(Fighandle,strcat('D:\Pictures\processed\Flow\Basic\Figure\FullRaster_inhib_clusters',num2str(counter),'.svg'),'-dsvg','-r0');
counter=counter+1;
end

%All brain regions
load('V:\Gilles\Zbrain_Masks.mat', 'Zbrain_Masks')
Zbrain_AllMask=vertcat(Zbrain_Masks{[1:1:77 79:1:294],3});
Zbrain_AllMask=unique(Zbrain_AllMask,'rows');


%Need to rotate ROI_temp 180
ROI_temp=round(ROI_fish);
dims_Zbrain=[1406 621];
ROI_rotated=ROI_temp;
ROI_rotated(:,2)=abs(ROI_temp(:,2)-dims_Zbrain(2));
ROI_rotated(:,1)=abs(ROI_temp(:,1)-dims_Zbrain(1));
ROI_rotated(:,3)=ROI_temp(:,3)/2;
clearvars ROI_temp


for i=1:length(GoodBetas_select)
    idx_temp=GoodClusters_goodmembers(i).idx;    
    CSV_temp=ROI_rotated(idx_temp,:);
    CSV_temp(:,3)=CSV_temp(:,3);
    CSV_temp(:,4)=1;
    filename=strcat('__Coords_clust_Basic',num2str(i),'b.csv');
    csvwrite(filename,CSV_temp);
end

for i=1:length(GoodBetas_inhib)
    idx_temp=GoodInhib_goodmembers(i).idx;   
    CSV_temp=ROI_rotated(idx_temp,:);
    CSV_temp(:,3)=CSV_temp(:,3);
    CSV_temp(:,4)=1;
    filename=strcat('__Coords_clust_BasicInhib',num2str(i),'.csv');
    csvwrite(filename,CSV_temp);
end

PerBrainRegions=struct();
RegionList={'Thalamus','Cerebellum','NucMLF','Semicircularis','Telencephalon','Tectum','Longitudinalis','Tegmentum','Habenula','Pretectum','MON','Hindbrain','pLLG'};
progressbar;
for i=1:length(RegionList)
    progressbar(i/length(RegionList));
    regionName=RegionList{i};
    if strcmp(regionName,'Telencephalon')
        Mask=Zbrain_Masks{294,3};
    elseif strcmp(regionName,'Hindbrain')
        Hindbrain_Mask=Zbrain_Masks{259,3};
        Mask=Zbrain_Masks{131,3};
        IsInEyes_temp=ismember(Hindbrain_Mask,Mask,'rows');IsInEyes_temp=find(IsInEyes_temp==1);%remove cerebellum
        Hindbrain_Mask(IsInEyes_temp,:)=[];
        Mask=Zbrain_Masks{295,3};
        IsInEyes_temp=ismember(Hindbrain_Mask,Mask,'rows');IsInEyes_temp=find(IsInEyes_temp==1);%remove MON
        Hindbrain_Mask(IsInEyes_temp,:)=[];
        Mask=Hindbrain_Mask;
        clearvars Hindbrain_Mask IsInEyes_temp;
	elseif strcmp(regionName,'pLLG')
        Mask=Zbrain_Masks{90,3};
    else
        Mask=[];
        IndexC=strfind({Zbrain_Masks{:,2}}, regionName);
        IndexC=find(not(cellfun('isempty', IndexC)));
        for j=IndexC
            if isempty(Mask)
                Mask=Zbrain_Masks{j,3};
            else
                Mask=vertcat(Mask,Zbrain_Masks{j,3});
            end
        end
    end
    Mask=unique(Mask,'rows');
    IsInBrainRegion=ismember(round(ROI_rotated),Mask,'rows');
    PerBrainRegions.(regionName).idx=find(IsInBrainRegion==1);
end
clearvars Mask IsInEyes_temp regionName i j k IndexC IsInBrainRegion

for j=1:length(RegionList)
    regionName=RegionList{j};
    temp=ZS(PerBrainRegions.(regionName).idx,:);
    coefficients=struct();
    rsquared=zeros(1,length(temp));
    parfor i=1:size(temp,1)
        mdl=fitlm(flow',temp(i,:));%,'interactions');
        coefficients(i).coef=mdl.Coefficients;
        rsquared(i)=mdl.Rsquared.Adjusted;
    end    
    LinReg.(regionName).coef=coefficients;
    LinReg.(regionName).rsquared=rsquared;    
end

%Kmeans of all
options = statset('UseParallel',0);
for j=1:length(RegionList)
    regionName=RegionList{j};
    temp=ZS(PerBrainRegions.(regionName).idx,:);
    idx_temp=find(LinReg.(regionName).rsquared>0.1);
    LinReg.(regionName).ZS_rsq=temp(idx_temp,:);
    if length(idx_temp)>1000
        [idxKmeans Cmap]=kmeans(LinReg.(regionName).ZS_rsq,20,'Replicates',5,'MaxIter',1000,'Display','final');
        LinReg.(regionName).KmeansCenter=Cmap;
        LinReg.(regionName).KmeansIdx=idxKmeans;
    elseif length(idx_temp)>100
        [idxKmeans Cmap]=kmeans(LinReg.(regionName).ZS_rsq,10,'Replicates',5,'MaxIter',1000,'Display','final');
        LinReg.(regionName).KmeansCenter=Cmap;
        LinReg.(regionName).KmeansIdx=idxKmeans;
    elseif length(idx_temp)>10
        [idxKmeans Cmap]=kmeans(LinReg.(regionName).ZS_rsq,2,'Replicates',5,'MaxIter',1000,'Display','final');
        LinReg.(regionName).KmeansCenter=Cmap;
        LinReg.(regionName).KmeansIdx=idxKmeans;
    else
        LinReg.(regionName).KmeansCenter=mean(temp,1);
        LinReg.(regionName).KmeansIdx=idx_temp;
    end
end

%Get the good ones
for j=1:length(RegionList)
    regionName=RegionList{j};
    [~,LinReg.(regionName).GoodBetas]=Test_Regress( LinReg.(regionName).KmeansCenter,flow,LinReg.(regionName).KmeansIdx,0.2);
end


%---------------------------------------------
%Remove ROIs from outside the brain
Zbrain_AllMask=vertcat(Zbrain_Masks{[1:1:77 79:1:294],3});
Zbrain_AllMask=unique(Zbrain_AllMask,'rows');
%Removing the eyes
idx_brain=ismember(round(ROI_rotated),Zbrain_AllMask,'rows');
idx_brain=find(idx_brain>0);

for i=1:length(GoodBetas_select)
    idx_temp=GoodClusters_goodmembers(i).idx;    
    CSV_temp=ROI_rotated(idx_temp,:);
    CSV_temp(:,3)=CSV_temp(:,3);
    CSV_temp(:,4)=1;
    filename=strcat('__Coords_clust_Basic',num2str(i),'b.csv');
    csvwrite(filename,CSV_temp);
end

CSV_temp=ROI_rotated;
CSV_temp(:,4)=1;
csvwrite('__All_fish.csv',CSV_temp);


for i=1:length(GoodBetas_inhib)
    idx_temp=GoodInhib_goodmembers(i).idx;   
    CSV_temp=ROI_rotated(idx_temp,:);
    CSV_temp(:,3)=CSV_temp(:,3);
    CSV_temp(:,4)=1;
    filename=strcat('__Coords_clust_BasicInhib',num2str(i),'.csv');
    csvwrite(filename,CSV_temp);
end


%---------------------------------------------
idx_PerBrainRegions=nan(size(idx_Fish));
RegionList_select={'Thalamus','Cerebellum','Semicircularis','Telencephalon','Tectum','Tegmentum','Habenula','Pretectum','MON','Hindbrain','pLLG'};
for i=1:length(RegionList_select)
    regionName=RegionList_select{i};
    idx_PerBrainRegions(PerBrainRegions.(regionName).idx)=i;
end
clearvars Mask IsInEyes_temp regionName i j k IndexC IsInBrainRegion

RegionList_select2={'Thalamus','Cereb','TS','Telencephalon','OT','Tegmentum','Habenula','Pretectum','MON','Hindbrain','pLLG'};
for i=1:length(GoodBetas_select)
    idx_temp=GoodClusters_goodmembers(i).idx;
    Fighandle=figure;set(Fighandle, 'Position', [0, 0, 1000, 500]);
    h=histogram(idx_PerBrainRegions(idx_temp),'Normalization','probability');
    h.FaceColor = colors(i,:)/256;
    h.EdgeColor = 'k';
    set(gca,'XTickLabel',RegionList_select2);set(gca,'XTick',[1:1:length(RegionList_select)]);set(gca,'XTickLabelRotation',-90)
    set(gca,'YTick',[0.1:0.1:0.5]);set(gca,'YTickLabel',[10:10:50]);ylim([0 45]);
    box off
    ax1 = gca;
    ax2 = axes('Position', get(ax1, 'Position'), 'FontSize', 10,...
        'Color','None','XColor','k','YColor','k', 'LineWidth', 1,...
        'XAxisLocation','top', 'XTick', [],...
        'YAxisLocation','right', 'YTick', []);
    linkaxes([ax1, ax2]);
    print(Fighandle,strcat('D:\Pictures\processed\Flow\Basic\Figure\BrainDistrib_clusters',num2str(i),'.png'),'-dpng','-r0');
    print(Fighandle,strcat('D:\Pictures\processed\Flow\Basic\Figure\BrainDistrib_clusters',num2str(i),'.svg'),'-dsvg','-r0');
    close;
end

idxKmeans_ZS_goodmembers_1To8=idxKmeans_ZS_goodmembers;
for i=1:length(GoodBetas_select)
    idxKmeans_ZS_goodmembers_1To8(idxKmeans_ZS_goodmembers==GoodBetas_select(i))=i;
end


for i=1:length(RegionList_select)
    Fighandle=figure;set(Fighandle, 'Position', [0, 0, 1000, 500]);
    regionName=RegionList_select{i};
    idx_temp=PerBrainRegions.(regionName).idx;
    idx_temp=idxKmeans_ZS_goodmembers_1To8(idx_temp);
    idx_temp(idx_temp==0)=nan;
    h=histogram(idx_temp,'Normalization','probability');
    h.FaceColor = [0.5 0.5 0.5];
    h.EdgeColor = 'k';
    set(gca,'XTick',[1:1:length(RegionList_select)]);
    set(gca,'YTick',[0.1:0.1:0.5]);set(gca,'YTickLabel',[10:10:50]);
    box off
    ax1 = gca;
    ax2 = axes('Position', get(ax1, 'Position'), 'FontSize', 10,...
        'Color','None','XColor','k','YColor','k', 'LineWidth', 1,...
        'XAxisLocation','top', 'XTick', [],...
        'YAxisLocation','right', 'YTick', []);
    linkaxes([ax1, ax2]);
    print(Fighandle,strcat('D:\Pictures\processed\Flow\Basic\Figure\ClustersDistrib_',regionName,'.png'),'-dpng','-r0');
    close;
end


CSV_temp=ROI_rotated(idx_PerBrainRegions==0,:);
CSV_temp(:,4)=1;
csvwrite('__All_fish.csv',CSV_temp);

CSV_temp=ROI_rotated(idx_PerBrainRegions>0,:);
CSV_temp(:,4)=1;
csvwrite('__All_fishb.csv',CSV_temp);

edges=[0:10:200];
figure;histogram(ROI_rotated(idx_PerBrainRegions>0,3),edges,'Normalization','probability');hold on;histogram(ROI_rotated(idx_PerBrainRegions==0,3),edges,'Normalization','probability');
figure;histogram(ROI_rotated(:,3),'Normalization','probability');
figure;histogram(ROI_rotated(idx_PerBrainRegions>0,3),edges,'Normalization','probability');hold on;histogram(Zbrain_AllMask(:,3),edges,'Normalization','probability');

DSI=zeros(length(unique(idx_Fish)),length(GoodBetas_select),2);
for i=1:size(DSI,1)
    for j=1:size(DSI,2)
        RespFWD=DS_raw(i,:,2,j);
        RespBWD=DS_raw(i,:,1,j);
        DSI(i,j,1)=nanmean((RespFWD-RespBWD)./(RespFWD+RespBWD));
        DSI(i,j,2)=nanstd((RespFWD-RespBWD)./(RespFWD+RespBWD));
    end
end

back=[55 255 455];
fwd=[155 355 555];
StimLength=100;
TimeToPeak_raw=NaN(length(unique(idx_Fish)),3,2,length(GoodBetas_select));counter2=1;
for j=GoodBetas_select
    idx_temp=find(idxKmeans_ZS_goodmembers==j);
    idx_fish_temp=idx_Fish(idx_temp);
    counter=1;
    for fish=(unique(idx_Fish)')
        if find(idx_fish_temp==fish)
            tempPlot=mean(ZS(idx_temp(find(idx_fish_temp==fish)),:),1);
            BackPlot=zeros(3,StimLength);
            FwdPlot=zeros(3,StimLength);
            for i=1:3
                BackPlot(i,:)=tempPlot(back(i):back(i)+99);
                FwdPlot(i,:)=tempPlot(fwd(i):fwd(i)+99);
            end
            [~,RespBWD]=max(BackPlot,[],2);
            [~,RespFWD]=max(FwdPlot,[],2);
            TimeToPeak_raw(counter,:,1,counter2)=RespBWD;
            TimeToPeak_raw(counter,:,2,counter2)=RespFWD;
        end
        counter=counter+1;
    end
    counter2=counter2+1;
end

TimeToPeak_AVG=zeros(length(unique(idx_Fish)),length(GoodBetas_select),2);
for i=1:size(TimeToPeak_AVG,1)
    for j=1:size(TimeToPeak_AVG,2)
        RespFWD=TimeToPeak_raw(i,:,2,j);
        RespBWD=TimeToPeak_raw(i,:,1,j);
        TimeToPeak_AVG(i,j,1)=nanmean(RespFWD);
        TimeToPeak_AVG(i,j,2)=nanmean(RespBWD);
    end
end


for i=1:length(GoodBetas_select)
    idx_temp=GoodClusters_goodmembers(i).idx;
    SpatialLoc_basic(i).ROIs=ROI_rotated(idx_temp,:);
    [SpatialLoc_basic(i).coeff,SpatialLoc_basic(i).score,~,~,SpatialLoc_basic(i).explained,~] = pca(ROI_rotated(idx_temp,:));
    SpatialLoc_basic(i).explained
end

for i=1:length(GoodBetas_select)
SpatialLoc_basic(i).explained(1)
end

ksTest_PCA_xyz_basic=struct();
for i=1:length(GoodBetas_select)
    for j=1:length(GoodBetas_select)
        [ksTest_PCA_xyz_basic.h(i,j),ksTest_PCA_xyz_basic.p(i,j),ksTest_PCA_xyz_basic.ks2stat(i,j)]=kstest2(SpatialLoc_basic(i).score(:,1),SpatialLoc_basic(j).score(:,1));
    end
end

p_val_list=[];counter=1;
for i=1:length(GoodBetas_select)
    for j=1:length(GoodBetas_select)
        if j>i
            p_val_list(counter)=ksTest_PCA_xyz_basic.p(i,j);
            counter=counter+1;
        end
    end
end
[corr_p_val_list,h]=bonf_holm(p_val_list);



counter=1;
for i=1:length(GoodBetas_select)
    for j=1:length(GoodBetas_select)
        if j>i
             counter=counter+1;
        end
    end
end

regionName='pLLG';
idx_temp=PerBrainRegions.(regionName).idx;
idx_K=idxKmeans_ZS_goodmembers_1To8(idx_temp);
GB_pLLG=unique(idx_temp);
GB_pLLG(GB_pLLG==0)=[];
figure;
pLLG_clusters=[];
for i=1:length(GB_pLLG)
    idx_Plot=find(idx_K==GB_pLLG(i));
    pLLG_clusters(i,:)=mean(ZS(idx_temp(idx_Plot),:),1);
    
end

figure;imagesc(ZS(idx_temp,:));