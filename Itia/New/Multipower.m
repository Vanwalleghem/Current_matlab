MatFiles=dir('*Audio*analysis_matlab.mat');
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

% MatFiles=dir('*analysis_matlab.mat');
% name=strcat(MatFiles(1).name);
% Baseline=load(name, 'Baseline');
% Baseline=Baseline.Baseline;
% Baseline=cell2mat(Baseline);
% Fitness=load(name, 'idx_components');
% Fitness=Fitness.idx_components+1;
% Baseline=Baseline(Fitness);
% for i = 2:length(MatFiles)
%     name=strcat(MatFiles(i).name);
%     B=load(name, 'Baseline');
%     B=B.Baseline;
%     B=cell2mat(B);
%     Fitness=load(name, 'idx_components');
%     Fitness=Fitness.idx_components+1;
%     B=B(Fitness);
%     Baseline=horzcat(Baseline,B);    
%     %MatFiles(i).Baseline=Baseline;    
% end
% clearvars GC C S F N name i GS Calcium Noise Fitness

ZS=zscore(GoodCalcium,1,2);
x = linspace(0.2,size(ZS,2)/5,size(ZS,2));y = linspace(1,size(ZS,1),size(ZS,1));
Fighandle=figure;
set(Fighandle, 'Position', [100, 100, 1200, 900]);
imagesc(x,y,ZS(randperm(size(ZS,1)),:),[-0.5 5]);colormap hot;set(gca,'YTickLabel',[]);

ZS2=zscore(GoodCalcium+GoodNoise,1,2);

%correction for 09111_ERO et 09112_ELO

IndexC=strfind({MatFiles.name}, Fish_list{1});
MatFiles_fish = find(not(cellfun('isempty', IndexC)));
if MatFiles_fish(1)==1
        numbersForROIs=[1 [MatFiles(MatFiles_fish).GoodNumber]];
    else
        numbersForROIs=[MatFiles(MatFiles_fish(1)-1).GoodNumber [MatFiles(MatFiles_fish).GoodNumber]];
end    
ZS_temp=ZS2(numbersForROIs(1):numbersForROIs(end)-1,:);
ZS_temp=ZS_temp(:,21:end);
ZS_temp(:,741:760)=ZS_temp(:,end-19:end);
ZS2(numbersForROIs(1):numbersForROIs(end)-1,:)=ZS_temp;

IndexC=strfind({MatFiles.name}, Fish_list{2});
MatFiles_fish = find(not(cellfun('isempty', IndexC)));
if MatFiles_fish(1)==1
        numbersForROIs=[1 [MatFiles(MatFiles_fish).GoodNumber]];
    else
        numbersForROIs=[MatFiles(MatFiles_fish(1)-1).GoodNumber [MatFiles(MatFiles_fish).GoodNumber]];
end    
ZS_temp=ZS2(numbersForROIs(1):numbersForROIs(end)-1,:);
ZS_temp=ZS_temp(:,21:end);
ZS_temp(:,741:760)=ZS_temp(:,end-19:end);
ZS2(numbersForROIs(1):numbersForROIs(end)-1,:)=ZS_temp;


Numbers=[1 [MatFiles.GoodNumber]];
counter=1;
idx_Plane=nan(length(ZS2),1);
idx_Position=nan(length(ZS2),1);
idx_Fish=nan(length(ZS2),1);
name=strcat(MatFiles(1).name);
for i=1:length(MatFiles)	
    name=strcat(MatFiles(i).name);
    [Plane,~]=regexp(name,'_(\d+)_','tokens','match');Plane=str2num(Plane{1}{1});
    if strfind(name,'ELO')
        Fish=1;
    elseif strfind(name,'ERO')
        Fish=2;  
    end
    idx_Plane(Numbers(i):Numbers(i+1))=Plane;
    idx_Position(Numbers(i):Numbers(i+1))=Fish;
    [Fish,~]=regexp(name,'Fish2017(\d+)_','tokens','match');Fish=str2num(Fish{1}{1});
    idx_Fish(Numbers(i):Numbers(i+1))=Fish;
end
clearvars i Fish Plane name counter

temp=unique(idx_Fish);
for i=1:length(temp);
    idx_Fish(find(idx_Fish==temp(i)))=i;
end

% options = statset('UseParallel',1); [idxKmeans_ZS Cmap_ZS]=kmeans(ZS,10,'Options',options,'Replicates',5,'MaxIter',1000,'Display','final');
% figure;imagesc(Cmap_ZS,[-0.5 3]);colormap hot

Stimuli=zeros(6,size(ZS2,2));
GCaMP6=[0,1.69644104899772,5.13796058542217,8.27886020152244,10.3756715204800,11.8173714529814,12.2425184714093,10.8571417354877,8.80831829681196,6.91339112244670,5.46959264663869,4.30868766622567,3.42533619066766,2.75378443486879,2.18017250852183,1.72816235135824,1.32732537295463,1.00684435500268,0.730210038304555,0.530242444093118,0.362253250339685,0.227668255288566,0.0869242416152502,0.000718266708050853,-0.0828334873368325]';
idxStart=40;
counter=0;
for i=1:18
    if mod(i,3)==1
        counter=counter+1;
    end
    Stimuli(counter,(idxStart+(i-1)*40):(idxStart+(i-1)*40)+size(GCaMP6,1)-1)=GCaMP6;
end

% ModelMultipower=[];
% parfor i=1:length(ZS)
%     %mdl=stepwiselm(Regressor',DataSet(i,:),'linear','Criterion','adjrsquared','Intercept',false,'Upper','interactions','Verbose',0);
%     mdl=fitlm(Stimuli',ZS(i,:));%,'interactions');
%     ModelMultipower(i).coef=mdl.Coefficients;
%     %ModelResultsSeg_ZS(i).MSE=mdl.MSE;
%     %ModelResultsSeg_ZS(i).Fitted=mdl.Fitted;
%     ModelMultipower(i).rsquared=mdl.Rsquared.Adjusted;
% end

ModelMultipower2=[];
parfor i=1:length(ZS2)
    %mdl=stepwiselm(Regressor',DataSet(i,:),'linear','Criterion','adjrsquared','Intercept',false,'Upper','interactions','Verbose',0);
    mdl=fitlm(Stimuli',ZS2(i,:));%,'interactions');
    ModelMultipower2(i).coef=mdl.Coefficients;
    %ModelResultsSeg_ZS(i).MSE=mdl.MSE;
    %ModelResultsSeg_ZS(i).Fitted=mdl.Fitted;
    ModelMultipower2(i).rsquared=mdl.Rsquared.Adjusted;
end

% coefficients={};
% for idx=1:length(ModelMultipower)
%     coef=[ModelMultipower(idx).coef];
%     temp=coef.Properties.RowNames;temp=regexp(temp,'x(\d+)','tokens');
%     if ~isempty(temp)
%         %temp=[temp{:}];temp=[temp{:}];temp=[temp{:}];%temp=str2num(temp);
%         for coef_idx=2:height(coef)
%             if coef.pValue(coef_idx)<0.05
%                 coefficients{idx,str2num(temp{coef_idx}{1}{1})}=coef.Estimate(coef_idx);
%             end
%         end
%     end
% end
% idxempty=cellfun('isempty',coefficients);
% coefficients(idxempty)={0};
% clearvars idxempty idx coef_idx coef temp
% coefficients=cell2mat(coefficients);

% figure;histogram([ModelMultipower.rsquared]);
% idx_rsq=find([ModelMultipower.rsquared]>0.1);
% options = statset('UseParallel',1); [idxKmeans_ZS_rsq Cmap_ZS_rsq]=kmeans(ZS(idx_rsq,:),30,'Options',options,'Distance','cityblock','Replicates',5,'MaxIter',1000,'Display','final');
% [Model_ZS,GoodBetas_ZS]=Test_Regress(Cmap_ZS_rsq,Stimuli,idxKmeans_ZS_rsq,0.4);

coefficients={};
for idx=1:length(ModelMultipower2)
    coef=[ModelMultipower2(idx).coef];
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
rsq_list=[ModelMultipower2.rsquared];


idx_rsq=find([ModelMultipower2.rsquared]>0.2);
ZS2_rsq=ZS2(idx_rsq,:);

options = statset('UseParallel',1); [idxKmeans_ZS_rsq Cmap_ZS_rsq]=kmeans(ZS2_rsq,10,'Options',options,'Distance','cityblock','Replicates',5,'MaxIter',1000,'Display','final');
[Model_ZS,GoodBetas_ZS]=Test_Regress(Cmap_ZS_rsq,Stimuli,idxKmeans_ZS_rsq,0.4);

options = statset('UseParallel',1); [idxKmeans_ZS_rsq2 Cmap_ZS_rsq2]=kmeans(ZS2_rsq,15,'Options',options,'Distance','cityblock','Replicates',5,'MaxIter',1000,'Display','final');
[Model_ZS2,GoodBetas_ZS2]=Test_Regress(Cmap_ZS_rsq2,Stimuli,idxKmeans_ZS_rsq2,0.4);

eva = evalclusters(ZS2(idx_rsq,:),'kmeans','CalinskiHarabasz','KList',[1:20]);
eva2 = evalclusters(ZS2(idx_rsq,:),'kmeans','DaviesBouldin','KList',[1:20]);


coef_stats=zeros(2,size(coefficients,2))
coef_stats(1,:)=mean(coefficients,1);
coef_stats(2,:)=std(coefficients,1,1);

% min(coefficients(:,1));
% idx_neg=find(coefficients(:,1)<-0.05);
% figure;imagesc(ZS2(idx_neg,:));
% options = statset('UseParallel',1); [idxKmeans_ZS_neg Cmap_ZS_neg]=kmeans(ZS2(idx_neg,:),2,'Options',options,'Replicates',5,'MaxIter',1000,'Display','final');
% figure;
% for i=1:5
%     plot(Cmap_ZS_neg(i,:));pause
% end


% ZS_AVG=zeros(size(ZS,1),246);
% ZS_STD=zeros(size(ZS,1),246);
% parfor idx_ZS=1:size(ZS,1)
%     start=30;
%     AVG=[];
%     for i=1:3
%         AVG(i,:)=ZS(idx_ZS,start:start+40);
%         start=start+40;
%     end
%     STD=std(AVG,1,1);
%     AVG=mean(AVG,1);
%     AVG=AVG-min(AVG);
%     j=1;
%     for j=2:6
%         for i=1:3
%             temp(i,:)=ZS(idx_ZS,start:start+40);
%             start=start+40;
%         end
%         temp_STD=std(temp,1,1);
%         temp=mean(temp,1);
%         temp=temp-min(temp);
%         STD=[STD temp_STD];
%         AVG=[AVG temp];
%     end
%     ZS_AVG(idx_ZS,:)=AVG;
%     ZS_STD(idx_ZS,:)=STD;
% end

ZS_AVG2=zeros(size(ZS2,1),246);
parfor idx_ZS=1:size(ZS2,1)
    start=30;
    AVG=[];
    for i=1:3
        AVG(i,:)=ZS2(idx_ZS,start:start+40);
        start=start+40;
    end    
    AVG=mean(AVG,1);
    AVG=AVG-min(AVG);
    j=1;
    for j=2:6
        for i=1:3
            temp(i,:)=ZS2(idx_ZS,start:start+40);
            start=start+40;
        end        
        temp=mean(temp,1);
        temp=temp-min(temp);       
        AVG=[AVG temp];
    end
    ZS_AVG2(idx_ZS,:)=AVG;   
end

Stimuli_AVG=zeros(6,size(ZS_AVG2,2));
GCaMP6=[0,1.69644104899772,5.13796058542217,8.27886020152244,10.3756715204800,11.8173714529814,12.2425184714093,10.8571417354877,8.80831829681196,6.91339112244670,5.46959264663869,4.30868766622567,3.42533619066766,2.75378443486879,2.18017250852183,1.72816235135824,1.32732537295463,1.00684435500268,0.730210038304555,0.530242444093118,0.362253250339685,0.227668255288566,0.0869242416152502,0.000718266708050853,-0.0828334873368325]';
idxStart=11;
for i=1:6    
    Stimuli_AVG(i,(idxStart+(i-1)*41):(idxStart+(i-1)*41)+size(GCaMP6,1)-1)=GCaMP6;
end
% 
% figure;plot(mean(ZS_AVG(idx_neg(find(idxKmeans_ZS_neg==2)),:),1));
% figure;imagesc(ZS_AVG(idx_neg(find(idxKmeans_ZS_neg==1)),:));
% temp=mean(ZS_AVG(idx_neg(find(idxKmeans_ZS_neg==2)),:),1);
% start=1;
% for i=1:6
%     temp(start:start+39)=temp(start:start+39)-temp(start);
%     start=start+40;
% end

% Fighandle=figure;
% set(Fighandle, 'Position', [100, 100, 1200, 300]);
% start=1;
% for i=1:6
%     subplot(1,6,i);plot(temp(start:start+39));ylim([-0.7 0.1]);
%     if i>1
%         set(gca,'YTickLabel',[]);
%     end
%     start=start+40;    
% end

% Fighandle=figure;
% set(Fighandle, 'Position', [100, 100, 1200, 300]);
% start=1;
% for i=1:6
%     subplot(1,6,i);plot(mean(ZS_AVG(idx_rsq(find(idxKmeans_ZS_rsq==GoodBetas_ZS(1))),start:start+39),1));ylim([0 4]);
%     start=start+40;    
% end

% figure;plot(mean(ZS_AVG(idx_rsq(find(idxKmeans_ZS_rsq==GoodBetas_ZS(1))),:),1));
% figure;plot(temp);

% options = statset('UseParallel',1); [idxKmeans_ZSAVG_rsq Cmap_ZSAVG_rsq]=kmeans(ZS_AVG(idx_rsq,:),30,'Options',options,'Distance','cityblock','Replicates',5,'MaxIter',1000,'Display','final');
% [Model_ZS,GoodBetas_ZS]=Test_Regress(Cmap_ZSAVG_rsq,Stimuli_AVG,idxKmeans_ZSAVG_rsq,0.4);

idxKmeans_final=zeros(1,size(ZS2,1));
%idxKmeans_final(idx_rsq)=idxKmeans_ZSAVG_rsq;
idxKmeans_final(idx_rsq)=idxKmeans_ZS_rsq2;

%GoodBetas=GoodBetas_ZS([1 3 4 6 8 9 10]);
%GoodBetas=GoodBetas_ZS3([2 4 5 6 7]);
GoodBetas=GoodBetas_ZS2([6 3 4 8]);

counter=1;
framerate=4;
x = linspace(1/framerate,size(ZS_AVG2,2)/framerate,size(ZS_AVG2,2));
rows=length(GoodBetas);xplot=4;
counter=1;
Fighandle=figure;
set(Fighandle, 'Position', [100, 100, 1300, 900]);
for i=GoodBetas
    idx_temp=find(idxKmeans_final==i);
    subplot(rows,xplot,counter);shadedErrorBar(x, mean(ZS_AVG2(idx_temp,:),1), std(ZS_AVG2(idx_temp,:),1,1));axis([0 60 -1 6]);
    subplot(rows,xplot,counter+1);imagesc(ZS2(idx_temp,:),[-1 4]); colormap hot
    subplot(rows,xplot,counter+2);histogram(idx_Plane(idx_temp));
    subplot(rows,xplot,counter+3);histogram(idx_Position(idx_temp),[0.5:1:3.5]);h = gca;h.XTick=[1 2];h.XTickLabel={'ELO','ERO'};
    %subplot(rows,xplot,counter+4);histogram(idx_Fish(idx_temp));%h = gca;h.XTickLabel={'1','2','3','4','5','6','7'};
    counter=counter+xplot;
end

% counter=1;
% for i=GoodBetas
%     [Clusters{counter,1} Clusters{counter,2}]=kmeans(ZS_AVG(find(idxKmeans_final==i),:),2,'Options',options,'Replicates',5,'MaxIter',1000,'Display','final');
%     counter=counter+1;
% end
% 
% figure;
% for i=1:length(GoodBetas)
%     plot(Clusters{i,2}');pause
% end

GoodClustersData=[];
for i=1:length(GoodBetas)
    GoodClustersData(i).ZS=ZS2(idxKmeans_final==GoodBetas(i),:);
    GoodClustersData(i).Mean=mean(GoodClustersData(i).ZS,1);
    GoodClustersData(i).STD=std(GoodClustersData(i).ZS,1,1);
end

% for i=1:numel(GoodClustersData)
%     corr_temp=zeros(size(GoodClustersData(i).ZS,1),1);
%     parfor j=1:size(GoodClustersData(i).ZS,1)
%         temp=corrcoef(GoodClustersData(i).Mean, GoodClustersData(i).ZS(j,:));
%         corr_temp(j)=temp(1,2);
%     end
%     GoodClustersData(i).CorrCoef=corr_temp;
% end

Correlation_Clusters={};Threshold=0.4;
for i=1:length(GoodBetas)    
    idx_temp=find(idxKmeans_final==GoodBetas(i));
    corr_temp=zeros(size(idx_temp));
    parfor j=1:length(idx_temp)
        %temp=corrcoef(Cmap_ZSAVG_rsq(GoodBetas(i),:), ZS_AVG(idx_temp(j),:));
        temp=corrcoef(GoodClustersData(i).Mean, ZS2(idx_temp(j),:));
        corr_temp(j)=temp(1,2);
    end    
    GoodClusters_goodmembers(i).ZS=ZS_AVG2(idx_temp(find(corr_temp>=Threshold)),:);        
    %GoodClusters_goodmembers(i).STD=ZS_STD(idx_temp(find(corr_temp>=Threshold)),:);    
    Correlation_Clusters{i}=corr_temp;
    idxKmeans_final_goodmember(idx_temp(find(corr_temp>=Threshold)))=GoodBetas(i);
end

figure;
for i=1:length(Correlation_Clusters)
    subplot(length(Correlation_Clusters),1,i);histogram(Correlation_Clusters{i});
end

GoodClusters_goodmembers=[];Threshold=[0.5 0.5 0.55 0.2];
idxKmeans_final_goodmember=zeros(size(idxKmeans_final));
for i=1:length(GoodBetas)      
    corr_temp=Correlation_Clusters{i};
    idx_temp=find(idxKmeans_final==GoodBetas(i));
    GoodClusters_goodmembers(i).ZS=ZS_AVG2(idx_temp(find(corr_temp>=Threshold(i))),:);            
    idxKmeans_final_goodmember(idx_temp(find(corr_temp>=Threshold(i))))=GoodBetas(i);
end

counter=1;
framerate=4;
x = linspace(1/framerate,size(ZS_AVG2,2)/framerate,size(ZS_AVG2,2));
rows=length(GoodBetas);
counter=1;counter2=1;
Fighandle=figure;xplot=4;
set(Fighandle, 'Position', [100, 100, 1300, 900]);
for i=GoodBetas
    idx_temp=find(idxKmeans_final_goodmember==i);
    subplot(rows,xplot,counter);shadedErrorBar(x, mean(GoodClusters_goodmembers(counter2).ZS,1), std(GoodClusters_goodmembers(counter2).ZS,1,1));axis([0 60 0 5]);
    subplot(rows,xplot,counter+1);imagesc(GoodClusters_goodmembers(counter2).ZS,[0 5]); colormap hot
    subplot(rows,xplot,counter+2);histogram(idx_Plane(idx_temp));
    subplot(rows,xplot,counter+3);histogram(idx_Position(idx_temp),[0.5:1:2.5]);h = gca;h.XTick=[1 2];h.XTickLabel={'ELO','ERO'};
    %subplot(rows,xplot,counter+4);histogram(idx_Fish(idx_temp));%h = gca;h.XTickLabel={'1','2','3','4','5','6','7'};
    counter=counter+xplot;
    counter2=counter2+1;
end

% colors = distinguishable_colors(length(GoodBetas),[1 1 1; 0 0 0]);
% %colors = [1 0 0; 0 1 0; 1 0 1];
% colors = colors*256;
% Fighandle=figure;
% set(Fighandle, 'Position', [100, 100, 1400, 500]);x = linspace(0.5,size(ZS_AVG2,2)/2,size(ZS_AVG2,2));
% counter=1;counter2=1;xplot=1;yplot=length(GoodBetas);%yplot=ceil(length(GoodBetas)/xplot);
% for i=GoodBetas
%     idx_temp=find(idxKmeans_final==i);
%     subplot(xplot,yplot,counter);plot(x,mean(ZS_AVG2(idx_temp,:),1),'color',colors(counter2,:)/256);axis([0 123 -0.5 5]);
%     start=30;    
%     counter=counter+1;
%     counter2=counter2+1;
% end

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
    temp{counter}=find(idxKmeans_final_goodmember==i);
    %tempidx=find(idxKmeans==idx);
    %temp{counter}=GoodClusters_goodmembers(counter).idx;
    counter=counter+1;    
end
% 
% for idx=1:length(MatFiles)
%     filename=MatFiles(idx).name;
%     ROIsNb=[];ClusterNb=[];
%     %for k = 1 : length(temp)
%     for k = 1 : length(temp)
%         tempROIsNb=find([temp{k}]<=Numbers(idx+1));
%         if tempROIsNb            
%             ROIsNb=[ROIsNb temp{k}(tempROIsNb)];
%             temp{k}(tempROIsNb)=[];
%             ClusterNb=[ClusterNb ; repmat(k,length(tempROIsNb),1)];
%         end
%     end
%     if ROIsNb
%         imagename=regexp(filename,'_output_analysis','split');
%         %imagename=regexp(imagename,'_output_analysis_matlab2.mat','split');
%         imagename=strcat(imagename{1},'_mean.tif');
%         image=double(imread(imagename));image=image/max(max(image));image=image*128;
%         image=uint8(image);
%         image2=zeros(size(image(:,:,1)));
%         image3=repmat(image,1,1,3);
%         ROIs=All_ROIs{idx};       
%         ROIsNb=ROIsNb-Numbers(idx);
%         ROIs=ROIs(:,ROIsNb);
%         for k = 1 : size(ROIs,2)
%             image2=zeros(size(image(:,:,1)));
%             ROI=full(reshape(ROIs(:,k),size(image,1),size(image,2)));            
%             image2=image2+ROI;image2=(image2/max(max(image2)));image2=uint8(image2);
%             for j=1:3
%                 image3(:,:,j)=image3(:,:,j)+image2*colors(ClusterNb(k),j);
%             end
%         end
%         %image3(:,:,3)=image;
%             name=strcat('_Kmeans_',imagename(4:end));
%     imwrite(image3,name,'tif');
%     end
%     %image3=uint8(image3);
% 
% end
% clearvars idx i temp tempFileNb fileNb AVG_files filename image counter Numbers image2 image3 k ROI ROIs ROIsNb Start tempROIsNb name imagename tempidx Raster


% Fighandle=figure;
% set(Fighandle, 'Position', [0, 0, 1280, 1024]);
% set(findall(Fighandle,'type','text'),'fontSize',12,'fontWeight','bold','FontName','Arial')
% counter=1;counter2=1;xplot=length(GoodBetas);yplot=6;
% StimLength=40;
% x = linspace(0.2,StimLength+1/5,StimLength+1);
% for i=GoodBetas
%     idx_temp=find(idxKmeans_final_goodmember==i);start=1; 
%     ZS_temp=ZS_AVG(idx_temp,:);
%     for j=1:yplot             
%         subplot(xplot,yplot,counter);
%         if (i==GoodBetas(3) | i==GoodBetas(7))
%             ZS_temp2=bsxfun(@minus,ZS_temp(:,start:start+StimLength),max(ZS_temp(:,start:start+StimLength),[],2));
%             H=shadedErrorBar(x, mean(ZS_temp2,1), std(ZS_temp2,1,1));axis([0 40 -3 1]);
%         else
%             H=shadedErrorBar(x, mean(ZS_AVG(idx_temp,start:start+StimLength),1), std(ZS_AVG(idx_temp,start:start+StimLength),1,1));axis([0 40 -1 5]);
%         end
%         counter2=floor((counter-1)/6)+1;
%         H.mainLine.Color=colors(counter2,:)/256;
%         H.patch.FaceColor=colors(counter2,:)/256;
%         H.edge(1).Color=colors(counter2,:)/256;
%         H.edge(2).Color=colors(counter2,:)/256;
%         start=start+StimLength;
%         counter=counter+1;    
%     end  
% end

% Multipower_clusters=nan(length(GoodBetas)*6,length(unique(idx_Fish)));
% counter=0;
% for i=GoodBetas
%     idx_temp=find(idxKmeans_final_goodmember==i);
%     ZS_temp=ZS_AVG(idx_temp,:);
%     if (i==GoodBetas(3) | i==GoodBetas(7))
%         start=1;
%         for powerLevel=1:6
%             ZS_temp(:,start:start+StimLength)=bsxfun(@minus,ZS_temp(:,start:start+StimLength),max(ZS_temp(:,start:start+StimLength),[],2));
%         end
%         start=start+StimLength;
%     end
%     for fish=unique(idx_Fish(idx_temp))'
%         start=1;
%         idx_temp2=find(idx_Fish(idx_temp)==fish);
%         for powerLevel=1:6
%             if (i==GoodBetas(3) | i==GoodBetas(7))
%                 Multipower_clusters(counter+powerLevel,fish)=mean(min(ZS_temp(idx_temp2,start:start+StimLength),[],2),1);
%             else
%                 Multipower_clusters(counter+powerLevel,fish)=mean(max(ZS_temp(idx_temp2,start:start+StimLength),[],2),1);
%             end
%             start=start+StimLength;
%         end
%     end
%     counter=counter+6;
% end
% 
% Multipower_clusters2=nan(length(GoodBetas)*6,length(unique(idx_Fish)));
% counter=0;
% for i=GoodBetas
%     idx_temp=find(idxKmeans_final_goodmember==i);
%     ZS_temp=ZS_AVG(idx_temp,:);
%     if (i==GoodBetas(3) | i==GoodBetas(7))
%         start=1;
%         for powerLevel=1:6
%             ZS_temp(:,start:start+StimLength)=bsxfun(@minus,ZS_temp(:,start:start+StimLength),max(ZS_temp(:,start:start+StimLength),[],2));
%         end
%         start=start+StimLength;
%     end
%     for fish=unique(idx_Fish(idx_temp))'
%         start=1;
%         idx_temp2=find(idx_Fish(idx_temp)==fish);
%         for powerLevel=1:6
%             if (i==GoodBetas(3) | i==GoodBetas(7))
%                 Multipower_clusters2(counter+powerLevel,fish)=mean(trapz(ZS_temp(idx_temp2,start:start+StimLength),2),1);
%             else
%                 Multipower_clusters2(counter+powerLevel,fish)=mean(trapz(ZS_temp(idx_temp2,start:start+StimLength),2),1);
%             end
%             start=start+StimLength;
%         end
%     end
%     counter=counter+6;
% end


colors = distinguishable_colors(length(GoodBetas),[1 1 1; 0 0 0]);
colors = colors*256;
Fighandle=figure;
set(Fighandle, 'Position', [100, 100, 1400, 500]);x = linspace(0.5,size(ZS_AVG2,2)/2,size(ZS_AVG2,2));
counter=1;counter2=1;xplot=1;yplot=length(GoodBetas);%yplot=ceil(length(GoodBetas)/xplot);
for i=GoodBetas
    idx_temp=find(idxKmeans_final_goodmember==i);
    subplot(xplot,yplot,counter);plot(x,mean(ZS_AVG2(idx_temp,:),1),'color',colors(counter2,:)/256);axis([0 123 -0.5 5]);
    start=30;    
    counter=counter+1;
    counter2=counter2+1;
end

% Numbers=[0 [ROIs_idx]];
% temp=[];
% counter=1;
% for i=GoodBetas
% temp{counter}=find(idxKmeans_final_goodmember==i);
% %tempidx=find(idxKmeans==idx);
% %temp{counter}=GoodClusters_goodmembers(counter).idx;
% counter=counter+1;
% end

for idx=1:length(MatFiles)
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
        name=strcat('_ThrKmeans_',imagename(4:end));
        imwrite(image3,name,'tif');
    end
    %image3=uint8(image3);
end
clearvars idx i temp tempFileNb fileNb AVG_files filename image counter Numbers image2 image3 k ROI ROIs ROIsNb Start tempROIsNb name imagename tempidx Raster

%colorblind colors
colors=[256 50 0; 240 159 0; 0 158 115; 0 114 250];

%Correct Inhib AVG
idx_temp=find(idxKmeans_final_goodmemberInBrain==GoodBetas(4));
temp_inhib=ZS2(idx_temp,:);
parfor idx_ZS=1:size(temp_inhib,1)
    start=30;
    AVG=[];
    for i=1:3
        AVG(i,:)=temp_inhib(idx_ZS,start:start+40);
        start=start+40;
    end
    AVG=mean(AVG,1);
    AVG=AVG-max(AVG);
    for j=2:6
        for i=1:3
            temp(i,:)=temp_inhib(idx_ZS,start:start+40)-prctile(temp_inhib(idx_ZS,start:start+40),90);
            start=start+40;
        end
        temp=mean(temp,1);
        temp=temp;
        AVG=[AVG temp];
    end
    temp_inhib_AVG(idx_ZS,:)=AVG;
end


%Make Figures
Fighandle=figure;
set(Fighandle, 'Position', [100, 100, 500, 1400]);x = linspace(0.5,size(ZS_AVG2,2)/2,size(ZS_AVG2,2));
counter=1;counter2=1;yplot=1;xplot=length(GoodBetas);%yplot=ceil(length(GoodBetas)/xplot);
for i=GoodBetas
    idx_temp=find(idxKmeans_final_goodmemberInBrain==i);
    subplot(xplot,yplot,counter);
    if counter ==4        
        H=shadedErrorBar(x,mean(temp_inhib_AVG,1), std(temp_inhib_AVG,1,1));axis([0 123 -3 6.5]);
        H.mainLine.Color=colors(counter2,:)/256;
        H.patch.FaceColor=colors(counter2,:)/512;
        H.edge(1).Color=colors(counter2,:)/512;
        H.edge(2).Color=colors(counter2,:)/512;
    else
        H=shadedErrorBar(x,mean(ZS_AVG2(idx_temp,:),1), std(ZS_AVG2(idx_temp,:),1,1));axis([0 123 -1 7]);
        H.mainLine.Color=colors(counter2,:)/256;
        H.patch.FaceColor=colors(counter2,:)/512;
        H.edge(1).Color=colors(counter2,:)/512;
        H.edge(2).Color=colors(counter2,:)/512;
    end
    counter=counter+1;
    counter2=counter2+1;
end
