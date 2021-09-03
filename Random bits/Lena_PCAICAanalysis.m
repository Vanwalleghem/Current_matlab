filelist=dir('*_2.tif.tif');
for File=1:length(filelist)
    [mixedsig, mixedfilters, CovEvals, covtrace, movm, movtm] = CellsortPCA(filelist(File).name,[],250,1,'C:\Temp\CROP');
    [ica_sig, ica_filters, ica_A, numiter] = CellsortICA(mixedsig, mixedfilters, CovEvals, [1:size(mixedsig,1)],  0.8,size(mixedsig,1),randn(size(mixedsig,1), size(mixedsig,1)),1e-6,50000);
    %PCA_ICA_results(File).ica_sig=ica_sig;
    %PCA_ICA_results(File).ica_filters=ica_filters;
    %[PCA_ICA_results(File).ica_segments, PCA_ICA_results(File).segmentlabel, PCA_ICA_results(File).segcentroid] = CellsortSegmentation(ica_filters, 2, 2, 20, 0)
    [PCA_ICA_results(File).ica_segments, PCA_ICA_results(File).segmentlabel, PCA_ICA_results(File).segcentroid] = CellsortSegmentation(ica_filters, 2, 2, 20, 0);    
    PCA_ICA_results(File).cell_sig = CellsortApplyFilter(filelist(File).name, PCA_ICA_results(File).ica_segments);
end
clearvars mixedsig ica_sig ica_filters ica_A numiter mixedsig mixedfilters CovEvals covtrace movm movtm
save('Results_CA8_2','-v7.3');clear all;


AllSegTraces=[];
idx=1;AllSegTraces=PCA_ICA_results(idx).cell_sig;
for idx=2:length(PCA_ICA_results)
    if size(PCA_ICA_results(idx).cell_sig,2)==750
        AllSegTraces=vertcat(AllSegTraces,PCA_ICA_results(idx).cell_sig(:,51:705));
    else
        AllSegTraces=vertcat(AllSegTraces,PCA_ICA_results(idx).cell_sig);
    end
end
%DF=DeltaF2(AllSegTraces,5,11);
ZS=zscore(AllSegTraces,1,2);
ZS=detrend(ZS')';

%idx_DF_5to200=find(max(DF,[],2)<5 & max(DF,[],2)>0.01);
%DF_select=DF(idx_DF_5to200,:);

%options = statset('UseParallel',1); [idxKmeans Cmap]=kmeans(DF_select,10,'Options',options,'Distance','cityblock','Replicates',5,'MaxIter',1000,'Display','final');
options = statset('UseParallel',1); [idxKmeans_ZS Cmap_ZS]=kmeans(ZS,20,'Options',options,'Distance','cityblock','Replicates',5,'MaxIter',1000,'Display','final');

figure;
for i=1:size(Cmap_ZS,1)
    plot(Cmap_ZS(i,:));pause;    
end

spike=[0.231927459, 1.649442638, 4.980908923, 1.490147774, 0.548712094, 0.234322333, 0.165871831, 0.039443932];
Lena_regressor=zeros(4,300);
for i = [68, 82, 115,150,157,189,200,218,232]
    Lena_regressor(1,i-1:i-1+length(spike)-1)=spike';
end
for i = [75, 89, 108,164,171,182,207,225,239]
    Lena_regressor(2,i-1:i-1+length(spike)-1)=spike';
end
for i = [64, 100, 122,157,174,178,189,196,210,214,228,232]
    Lena_regressor(3,i-1:i-1+length(spike)-1)=spike';
end
for i = [60, 96, 104,150,167,178,185,196,200,214,218,242]
    Lena_regressor(4,i-1:i-1+length(spike)-1)=spike';
end

[Model_ZS,GoodBetas_ZS]=Test_Regress(Cmap_ZS,Lena_regressor,idxKmeans_ZS,0.05);

for File=1:length(filelist)
    PCA_ICA_results(File).Numbers=size(PCA_ICA_results(File).cell_sig,1);
end

Numbers=cumsum([PCA_ICA_results.Numbers]);
Numbers=[1 Numbers];
counter=1;
idx_Plane=nan(length(DF),1);
name=strcat(filelist(1).name);
[Plane,~]=regexp(name,'_(\d+)','tokens','match');Plane=str2num(Plane{1}{1});
idx_Plane(1:Numbers(2))=Plane;
for i=2:length(filelist)    
    name=strcat(filelist(i).name);    
    [Plane,~]=regexp(name,'Slice(\d+)_','tokens','match');
    Plane=str2num(Plane{1}{1});
    idx_Plane(Numbers(i):Numbers(i+1))=Plane;    
end



GoodBetas=[1 3 8 10];
GoodBetas_select=GoodBetas;
colors = distinguishable_colors(length(GoodBetas),[1 1 1; 0 0 0]);

GoodClustersData=[];
for i=1:length(GoodBetas)
    GoodClustersData(i).ZS=ZS(idxKmeans_DF==GoodBetas(i),:);
    GoodClustersData(i).DF=DF(idxKmeans_DF==GoodBetas(i),:)*100;
    GoodClustersData(i).Mean=mean(GoodClustersData(i).DF,1);
    GoodClustersData(i).STD=std(GoodClustersData(i).DF,1,1);
end

Fighandle=figure;
set(Fighandle, 'Position', [100, 100, 1400, 900]);x = linspace(0.2,size(ZS,2)/5,size(ZS,2));
counter=1;counter2=1;xplot=floor(sqrt(length(GoodBetas_select)));yplot=ceil(length(GoodBetas_select)/xplot);
for i=1:length(GoodBetas_select)    
    subplot(length(GoodBetas_select),3,counter);plot(x,mean(zscore(GoodClustersData(i).DF,1,2),1),'color',colors(i,:),'LineWidth',2);
    subplot(length(GoodBetas_select),3,counter+1);imagesc(zscore(GoodClustersData(i).DF,1,2),[0 10]);colormap hot
    subplot(length(GoodBetas_select),3,counter+2);histogram(idx_Plane(idxKmeans_DF==GoodBetas(i)));
    counter=counter+3;
    counter2=counter2+1;
end

%idxKmeans_final=zeros(1,size(DF,1));
%idxKmeans_final(idx_rsq)=idxKmeans2;
idxKmeans_final=idxKmeans_DF;

File=1
load(strcat(filelist(File).name(1:end-3),'mat'), 'PCA_ICA_results');
Traces=PCA_ICA_results.Cell_sig;AllSegTraces=Traces;
Numbers=size(PCA_ICA_results.Cell_sig,1);
for File=2:length(filelist)
    load(strcat(filelist(File).name(1:end-3),'mat'), 'PCA_ICA_results');
    Traces=PCA_ICA_results.Cell_sig;
    Number=size(PCA_ICA_results.Cell_sig,1);
    Numbers=[Numbers Number];
    if size(Traces,2)==750
        AllSegTraces=vertcat(AllSegTraces,Traces(:,51:705));
    else
        AllSegTraces=vertcat(AllSegTraces,Traces);
    end
end

Number=Numbers;Number(1)=1;
colors = distinguishable_colors(length(GoodBetas),[1 1 1; 0 0 0]);
%colors = colors*256;
%Number=[1 Number];
%idxKmeans_final(idx_rsq_ZS)=idxKmeans_ZS_goodmembers;
for File=1:length(filelist)
    %load(strcat(filelist(File).name(1:end-3),'mat'), 'PCA_ICA_results');
    ROIs=PCA_ICA_results(File).ica_segments;
    idx_temp=idxKmeans_final(Number(File)+1:Number(File+1));
    imagename=strcat('AVG_',filelist(File).name);
    image=double(imread(imagename));image(image==max(max(image)))=0;image=image/max(prctile(image,95));image=image*64;
    image=uint8(image);    
    image3=repmat(image,1,1,3);    
    for i=GoodBetas
        idx_ROI=find(idx_temp==i);
        image2=squeeze(sum(ROIs(idx_ROI,:,:),1));
        image2=(image2/max(max(image2)))*200;image2=uint8(image2);
        for j=1:3
            image3(:,:,j)=image3(:,:,j)+image2*colors(find(GoodBetas==i),j);
        end
    end
    name=strcat('_Kmeans_good',imagename(4:end));
    imwrite(image3,name,'tif');
end

GoodBetas=[1 3 5 2 4];
x = linspace(0.2,size(ZS,2)/5,size(ZS,2));
Fighandle=figure;
set(Fighandle, 'Position', [100, 100, 1300, 900]);
counter=1;xplot=floor(sqrt(length(GoodBetas)));yplot=ceil(length(GoodBetas)/xplot);
for i=GoodBetas    
    NumberOfCells=length(find(idxKmeans_final==i));
    %NumberOfCells=length(find(idxKmeans_ZS_goodmembers==i));
    %subplot(5,1,counter);plot(x,Cmap_s(i,:),'color',colors(counter,:));title(num2str(NumberOfCells))
    subplot(5,1,counter);plot(x,GoodClusters_goodmembers(i).mean,'color',colors(counter,:),'LineWidth',2);axis([0 131 -2 4]);rectangle('FaceColor','r','Position',[11 -2 10 0.25]);rectangle('FaceColor','r','Position',[51 -2 10 0.25]);rectangle('FaceColor','r','Position',[91 -2 10 0.25]);rectangle('FaceColor','b','Position',[31 -2 10 0.25]);rectangle('FaceColor','b','Position',[71 -2 10 0.25]);rectangle('FaceColor','b','Position',[111 -2 10 0.25]);;title(num2str(NumberOfCells))
%    xlim([0 size(Cmap_ZS,2)])
    counter=counter+1;
end
clearvars idx i temp tempFileNb fileNb AVG_files filename image counter Numbers image2 image3 k ROI ROIs ROIsNb Start tempROIsNb name imagename tempidx Raster

