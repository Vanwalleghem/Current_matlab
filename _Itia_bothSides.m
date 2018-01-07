MatFiles=dir('*analysis_matlab.mat');
name=strcat(MatFiles(1).name);
Calcium=load(name, 'DenoisedTraces');
Calcium=Calcium.DenoisedTraces;
MatFiles(1).number=size(Calcium,1);
%Spikes=load(name, 'Spikes');
%Spikes=Spikes.Spikes;
Noise=load(name, 'Noise');
Noise=Noise.Noise;
Fitness=load(name, 'idx_components');
Fitness=Fitness.idx_components+1;
GoodCalcium=Calcium(Fitness,:);
%GoodSpikes=Spikes(Fitness,:);
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
    %S=load(name, 'Spikes');
    %S=S.Spikes;
    N=load(name, 'Noise');
    N=N.Noise;
    F=load(name, 'idx_components');
    F=F.idx_components+1;
    GC=C(F,:);
    GN=N(F,:);
    %GS=S(F,:);
    Noise=vertcat(Noise,N);
    Calcium=vertcat(Calcium,C);
    %Spikes=vertcat(Spikes,S);
    Fitness=horzcat(Fitness,F);
    GoodCalcium=vertcat(GoodCalcium,GC);
    GoodNoise=vertcat(GoodNoise,GN);
    %GoodSpikes=vertcat(GoodSpikes,GS);
    MatFiles(i).number=size(Calcium,1);
    MatFiles(i).GoodNumber=MatFiles(i-1).GoodNumber+length(F);
    %MatFiles(i).GC=GC;
end
clearvars GC C S F N name i GS Calcium Noise Fitness

ZS=zscore(GoodCalcium,1,2);
x = linspace(0.2,size(ZS,2)/5,size(ZS,2));y = linspace(1,size(ZS,1),size(ZS,1));
Fighandle=figure;
set(Fighandle, 'Position', [100, 100, 1200, 900]);
imagesc(x,y,ZS(randperm(size(ZS,1)),:),[-0.5 5]);colormap hot;set(gca,'YTickLabel',[]);

Numbers=[1 [MatFiles.GoodNumber]];
counter=1;
idx_Plane=nan(length(ZS),1);
idx_Position=nan(length(ZS),1);
idx_Fish=nan(length(ZS),1);
name=strcat(MatFiles(1).name);
for i=1:length(MatFiles)	
    name=strcat(MatFiles(i).name);
    [Plane,~]=regexp(name,'_(\d+)_','tokens','match');Plane=str2num(Plane{1}{1});
    if strfind(name,'Center')
        Fish=1;
    elseif strfind(name,'ELO')
        Fish=2;
    elseif strfind(name,'ERO')
        Fish=3;    
    elseif strfind(name,'Out')
        Fish=4;
    end
    idx_Plane(Numbers(i):Numbers(i+1))=Plane;
    idx_Position(Numbers(i):Numbers(i+1))=Fish;
    [Fish,~]=regexp(name,'Fish2017(\d+)_','tokens','match');Fish=str2num(Fish{1}{1});
    idx_Fish(Numbers(i):Numbers(i+1))=Fish;
end
clearvars i Fish Plane name counter

options = statset('UseParallel',1); [idxKmeans_ZS Cmap_ZS]=kmeans(ZS,10,'Options',options,'Replicates',5,'MaxIter',1000,'Display','final');
figure;imagesc(Cmap_ZS,[-0.5 3]);colormap hot
%ZS_max=max(ZS,[],2);idx_ZS_max=find(ZS_max>1);
%options = statset('UseParallel',1); [idxKmeans_ZS2 Cmap_ZS2]=kmeans(ZS(idx_ZS_max,:),20,'Options',options,'Distance','Correlation','Replicates',3,'MaxIter',1000,'Display','final');

ZS_8=ZS(find(idxKmeans_ZS==8),:);
options = statset('UseParallel',1); [idxKmeans_ZS_8 Cmap_ZS_8]=kmeans(ZS_8,10,'Options',options,'Distance','Correlation','Replicates',3,'MaxIter',1000,'Display','final');
% tSNE_Bothsides=tsne(ZS_9,[],3,20,100);
% options = statset('UseParallel',1); [idxKmeans_ZS_9_t Cmap_ZS_9_t]=kmeans(tSNE_Bothsides,5,'Options',options,'Replicates',3,'MaxIter',1000,'Display','final');
% Cmap_ZS_9_t=[];
% for i=1:max(idxKmeans_ZS_9_t)
%     Cmap_ZS_9_t(i,:)=mean(ZS_9(find(idxKmeans_ZS_9_t==i),:),1);
% end

% ZS_diff_rsq = zeros(size(ZS_9));
% parfor i=1:size(ZS_9,1)
%     temp = TVRegDiff(ZS_9(i,:), 50, 1e-1,[],'large',1e-8,[],0,0);
%     ZS_diff_rsq(i,:)=temp;
% end
% options = statset('UseParallel',1); [idxKmeans_ZS_9_t Cmap_ZS_9_t]=kmeans(ZS_diff_rsq,10,'Distance','Correlation','Options',options,'Replicates',3,'MaxIter',1000,'Display','final');

figure;
for i=1:size(Cmap_ZS_8,1)
    plot(Cmap_ZS_8(i,:));pause
end

Stimuli=zeros(3,size(ZS,2));
start=4;
spike=[0,1.69644104899772,5.13796058542217,8.27886020152244,10.3756715204800,12.2425184714093,10.8571417354877,8.80831829681196,6.91339112244670,5.46959264663869,4.30868766622567,3.42533619066766,2.75378443486879,2.18017250852183,1.72816235135824,1.32732537295463,1.00684435500268,0.730210038304555,0.530242444093118,0.362253250339685,0.0869242416152502,0.000718266708050853]';
for i = 1:3
    for j = 0:2
        Stimuli(i,start+j*60:start+j*60+length(spike)-1)=spike';
    end
    start=start+20;
end

ModelBothSides=[];
parfor i=1:length(ZS)
    %mdl=stepwiselm(Regressor',DataSet(i,:),'linear','Criterion','adjrsquared','Intercept',false,'Upper','interactions','Verbose',0);
    mdl=fitlm(Stimuli',ZS(i,:),'interactions');
    ModelBothSides(i).coef=mdl.Coefficients;
    %ModelResultsSeg_ZS(i).MSE=mdl.MSE;
    %ModelResultsSeg_ZS(i).Fitted=mdl.Fitted;
    ModelBothSides(i).rsquared=mdl.Rsquared.Adjusted;
end
idx_rsq=find([ModelBothSides.rsquared]>0.15);
ZS_rsq=ZS(idx_rsq,:);

Fighandle=figure;
set(Fighandle, 'Position', [100, 100, 1200, 900]);
imagesc(x,y,ZS_rsq(randperm(size(ZS_rsq,1)),:),[-0.5 5]);colormap hot;set(gca,'YTickLabel',[]);

coefficients={};
for idx=1:length(ModelBothSides)
    coef=[ModelBothSides(idx).coef];
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

options = statset('UseParallel',1); [idxKmeans_ZS_rsq Cmap_ZS_rsq]=kmeans(ZS_rsq,30,'Options',options,'Distance','cityblock','Replicates',3,'MaxIter',1000,'Display','final');

GoodBetas=[5 7 17 19 23 24 26 28];

idxKmeans_final=zeros(size(ZS,1),1);
idxKmeans_final(idx_rsq)=idxKmeans_ZS_rsq;

counter=1;
x = linspace(0.5,size(ZS,2)/2,size(ZS,2));
rows=length(GoodBetas);
counter=1;
Fighandle=figure;
set(Fighandle, 'Position', [100, 100, 1300, 900]);
for i=GoodBetas
    idx_temp=find(idxKmeans_final==i);
    subplot(rows,5,counter);plot(mean(ZS(idx_temp,:),1));
    subplot(rows,5,counter+1);imagesc(ZS(idx_temp,:),[-0.5 4]);
    subplot(rows,5,counter+2);histogram(idx_Plane(idx_temp));
    subplot(rows,5,counter+3);histogram(idx_Position(idx_temp),[0.5:1:4.5]);h = gca;h.XTick=[1 2 3 4];h.XTickLabel={'Center','ELO-IRO','ERO-ILO','Out'};
    subplot(rows,5,counter+4);histogram(idx_Fish(idx_temp));h = gca;h.XTickLabel={'1','2'};
    counter=counter+5;
end

GoodClustersData=[];
for i=1:length(GoodBetas)
    GoodClustersData(i).ZS=ZS(idxKmeans_final==GoodBetas(i),:);
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
temp=find(idxKmeans_final==GoodBetas(i));
GoodClusters_goodmembers(i).idx=temp(find(GoodClustersData(i).CorrCoef>=0.5));
GoodClusters_goodmembers(i).mean=mean(GoodClusters_goodmembers(i).ZS,1);
GoodClusters_goodmembers(i).STD=std(GoodClusters_goodmembers(i).ZS,1,1);
idx=find(idxKmeans_final==GoodBetas(i));
idx=idx(find(GoodClustersData(i).CorrCoef>=Threshold));
idxKmeans_ZS_goodmembers(idx)=GoodBetas(i);
%GoodClusters_goodmembers(i).Fish=idx_Fish(idx);
end

counter=1;
x = linspace(0.5,size(ZS,2)/2,size(ZS,2));
rows=length(GoodBetas);
counter=1;
Fighandle=figure;
set(Fighandle, 'Position', [100, 100, 1300, 900]);
for i=GoodBetas
    idx_temp=find(idxKmeans_ZS_goodmembers==i);
    subplot(rows,5,counter);plot(mean(ZS(idx_temp,:),1));
    subplot(rows,5,counter+1);imagesc(ZS(idx_temp,:),[-0.5 4]);
    subplot(rows,5,counter+2);histogram(idx_Plane(idx_temp));
    subplot(rows,5,counter+3);histogram(idx_Position(idx_temp),[0.5:1:4.5]);h = gca;h.XTick=[1 2 3 4];h.XTickLabel={'Center','ELO-IRO','ERO-ILO','Out'};
    subplot(rows,5,counter+4);histogram(idx_Fish(idx_temp));h = gca;h.XTickLabel={'1','2'};
    counter=counter+5;
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
for i=GoodBetas
    temp{counter}=find(idxKmeans_ZS_goodmembers==i);
    %tempidx=find(idxKmeans==idx);
    %temp{counter}=GoodClusters_goodmembers(counter).idx;
    counter=counter+1;    
end

colors = distinguishable_colors(length(GoodBetas),[1 1 1; 0 0 0]);
colors = colors*256;
for idx=1:length(MatFiles)
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
set(Fighandle, 'Position', [100, 100, 500, 1400]);x = linspace(0.5,size(ZS,2)/2,size(ZS,2));
counter=1;counter2=1;xplot=1;yplot=length(GoodBetas);%yplot=ceil(length(GoodBetas)/xplot);
coloring={'r','g','y'};
for i=GoodBetas
    idx_temp=find(idxKmeans_ZS_goodmembers==i);
    subplot(xplot,yplot,counter);plot(x,mean(ZS(idx_temp,:),1),'color',colors(counter2,:)/256);axis([0 131 -1 4]);
    start=30;
    for k=1:3        
        for j=1:3            
            rectangle('FaceColor',coloring{k},'Position',[start+(j-1)*30 -1 2 0.25]);
        end
        start=10+start;
    end    
    counter=counter+1;
    counter2=counter2+1;
end
