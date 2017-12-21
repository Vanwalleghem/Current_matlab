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
MatFiles(1).GoodNumber=length(Fitness);
MatFiles(1).GC=GoodCalcium;
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
    Noise=vertcat(Noise,N);
    Calcium=vertcat(Calcium,C);
    Spikes=vertcat(Spikes,S);
    Fitness=horzcat(Fitness,F);
    GoodCalcium=vertcat(GoodCalcium,GC);
    GoodSpikes=vertcat(GoodSpikes,GS);
    MatFiles(i).number=size(Calcium,1);
    MatFiles(i).GoodNumber=MatFiles(i-1).GoodNumber+length(F);
    MatFiles(i).GC=GC;
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
set(Fighandle, 'Position', [100, 100, 1200, 900]);
imagesc(x,y,ZS(randperm(size(ZS,1)),:),[-0.5 5]);colormap hot;set(gca,'YTickLabel',[]);

Numbers=[1 [MatFiles.GoodNumber]];
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
   
    idx_Plane(Numbers(i):Numbers(i+1))=Plane;
    idx_Fish(Numbers(i):Numbers(i+1))=Fish;
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


GoodClusters_goodmembers=[];Threshold=0.5;
idxKmeans_ZS_goodmembers=zeros(1,size(GoodCalcium,1));
for i=1:length(GoodBetas_select)
%GoodClusters_goodmembers(i).Spikes=GoodClustersData(i).Spikes(find(GoodClustersData(i).CorrCoef>=0.5),:);
%GoodClusters_goodmembers(i).ZS=zscore(GoodClustersData(i).DF(find(GoodClustersData(i).CorrCoef>=0.5),:),1,2);
GoodClusters_goodmembers(i).ZS=GoodClustersData(i).ZS(find(GoodClustersData(i).CorrCoef>=Threshold),:);
temp=find(idxKmeans==GoodBetas_select(i));
GoodClusters_goodmembers(i).idx=temp(find(GoodClustersData(i).CorrCoef>=0.5));
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

