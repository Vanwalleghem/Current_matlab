%---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
% Read the output of ANTs to extract the ROIs centroid
%

CSV_Files=dir('_2Warped*.csv');
ROIs=struct();
for i=1:length(CSV_Files);
    temp=csvread(CSV_Files(i).name,1);   
    Fishname=regexp(CSV_Files(i).name,'FishResized(\d+_\w+).csv','tokens');Fishname=Fishname{1}{1};    
    ROIs(i).name=Fishname;    
    ROIs(i).coord=temp(:,1:3);
    ROIs(i).idx=temp(:,5);    
end
clearvars i temp CSV_Files Fishname

load('D:\Pictures\processed\Itia\_BrainRegAndROIs.mat','Zbrain_Masks');

MatFiles=dir('*analysis_matlab.mat');
name=strcat(MatFiles(1).name);
Calcium=load(name, 'DenoisedTraces');
Calcium=Calcium.DenoisedTraces;
MatFiles(1).number=size(Calcium,1);
Noise=load(name, 'Noise');
Noise=Noise.Noise;
for i = 2:length(MatFiles)
    name=strcat(MatFiles(i).name);
    C=load(name, 'DenoisedTraces');
    C=C.DenoisedTraces;
    N=load(name, 'Noise');
    N=N.Noise;   
    Noise=vertcat(Noise,N);
    Calcium=vertcat(Calcium,C);
    MatFiles(i).number=size(Calcium,1);    
end
clearvars GC C S F N name i GS Fitness
ZS=zscore(Calcium+Noise,1,2);
clearvars GoodCalcium GoodNoise Calcium Noise

idx_Fish_name={};
for i=1:length(MatFiles)
    name=strcat(MatFiles(i).name);
    [Fish,~]=regexp(name,'Fish2017(\d+_E\DO)_','tokens');Fish=Fish{1};
    if iscell(Fish)
        Fish=Fish{1};
    end
    idx_Fish_name{i}=Fish;
end
clearvars i Fish Plane name counter
Fish_list=unique(idx_Fish_name);

PerBrainRegions=struct();
ItiaList={'Thalamus','Cerebellum','NucMLF','Semicircularis','Telencephalon','Tectum','Longitudinalis','Tegmentum','Habenula'};
progressbar;
for i=1:9%length(ItiaList)
    progressbar(i/length(ItiaList),[]);
    regionName=ItiaList{i};
    Mask=[];
    IndexC=strfind({Zbrain_Masks{:,2}}, regionName);
    IndexC=find(not(cellfun('isempty', IndexC)));
    if i==5
        IndexC=294;
    end
    for j=IndexC
        if isempty(Mask)
            Mask=Zbrain_Masks{j,3};
        else
            Mask=vertcat(Mask,Zbrain_Masks{j,3});
        end
    end
    for fish_nb=1:length(Fish_list)
        progressbar([],fish_nb/length(Fish_list));
        if iscell(Fish_list)
            Fish_name=Fish_list{fish_nb};
        else
            Fish_name=num2str(Fish_list(fish_nb));
        end
        IndexC=strfind({MatFiles.name}, Fish_name);
        MatFiles_fish = find(not(cellfun('isempty', IndexC)));
        ROI_name=strsplit(Fish_name,'Fish2017');
        while iscell(ROI_name)
            ROI_name=ROI_name{1};
        end
        IndexC=strfind({ROIs.name},ROI_name);
        ROI_fish=find(not(cellfun('isempty', IndexC)));        
        ROI_fish=ROIs(ROI_fish).coord;ROI_fish(:,1:2)=round(ROI_fish(:,1:2));
        ROI_fish(:,3)=round(((ROI_fish(:,3)-1)*2)+24);
        if MatFiles_fish(1)==1            
            numbersForROIs=[1 [MatFiles(MatFiles_fish).number]];
        else            
            numbersForROIs=[MatFiles(MatFiles_fish(1)-1).number+1 [MatFiles(MatFiles_fish).number]];
        end
        IsInBrainRegion=[];        
        IsInBrainRegion=ismember(ROI_fish,Mask,'rows');
        PerBrainRegions(fish_nb).(regionName).ROIsCent=ROI_fish(IsInBrainRegion,:);        
        temp=ZS(numbersForROIs(1):numbersForROIs(end),:);IsInBrainRegion=temp(IsInBrainRegion,:);        
        PerBrainRegions(fish_nb).(regionName).ZS=IsInBrainRegion;        
    end
end

ItiaList={'Thalamus','Cerebellum','NucMLF','Semicircularis','Telencephalon','Tectum','Longitudinalis','Tegmentum','Habenula','Hindbrain','MON'};
i=11;
regionName=ItiaList{i};
Mask=MON_Masks{1,3};
for fish_nb=1:length(Fish_list)
    progressbar([],fish_nb/length(Fish_list));
    if iscell(Fish_list)
        Fish_name=Fish_list{fish_nb};
    else
        Fish_name=num2str(Fish_list(fish_nb));
    end
    IndexC=strfind({MatFiles.name}, Fish_name);
    MatFiles_fish = find(not(cellfun('isempty', IndexC)));
    ROI_name=strsplit(Fish_name,'Fish2017');
    while iscell(ROI_name)
        ROI_name=ROI_name{1};
    end
    IndexC=strfind({ROIs.name},ROI_name);
    ROI_fish=find(not(cellfun('isempty', IndexC)));    
    ROI_fish=ROIs(ROI_fish).coord;ROI_fish(:,1:2)=round(ROI_fish(:,1:2));
    ROI_fish(:,3)=round(((ROI_fish(:,3)-1)*2)+24);
    if MatFiles_fish(1)==1        
        numbersForROIs=[1 [MatFiles(MatFiles_fish).number]];
    else        
        numbersForROIs=[MatFiles(MatFiles_fish(1)-1).number+1 [MatFiles(MatFiles_fish).number]];
    end
    IsInBrainRegion=[];
    IsInBrainRegion=ismember(ROI_fish,Mask,'rows');
    PerBrainRegions(fish_nb).(regionName).ROIsCent=ROI_fish(IsInBrainRegion,:);
    temp=ZS(numbersForROIs(1):numbersForROIs(end),:);IsInBrainRegion=temp(IsInBrainRegion,:);    
    PerBrainRegions(fish_nb).(regionName).ZS=IsInBrainRegion;
end

Hindbrain_Mask=Zbrain_Masks{259,3};
%Removes cerebellum
Mask=Zbrain_Masks{131,3};
IsInEyes_temp=ismember(Hindbrain_Mask,Mask,'rows');IsInEyes_temp=find(IsInEyes_temp==1);
Hindbrain_Mask(IsInEyes_temp,:)=[];
%Removes MON
Mask=MON_Masks{1,3};
IsInEyes_temp=ismember(Hindbrain_Mask,Mask,'rows');IsInEyes_temp=find(IsInEyes_temp==1);
Hindbrain_Mask(IsInEyes_temp,:)=[];
clearvars i j fish_nb Mask

regionName='Hindbrain';
 for fish_nb=1:length(Fish_list)
        progressbar(fish_nb/length(ItiaList));        
        if iscell(Fish_list)
            Fish_name=Fish_list{fish_nb};
        else
            Fish_name=num2str(Fish_list(fish_nb));
        end
        IndexC=strfind({MatFiles.name}, Fish_name);
        MatFiles_fish = find(not(cellfun('isempty', IndexC)));
        ROI_name=strsplit(Fish_name,'Fish2017');
        while iscell(ROI_name)
            ROI_name=ROI_name{1};
        end
        IndexC=strfind({ROIs.name},ROI_name);
        ROI_fish=find(not(cellfun('isempty', IndexC)));        
        ROI_fish=ROIs(ROI_fish).coord;ROI_fish(:,1:2)=round(ROI_fish(:,1:2));
        ROI_fish(:,3)=round(((ROI_fish(:,3)-1)*2)+24);
        if MatFiles_fish(1)==1           
            numbersForROIs=[1 [MatFiles(MatFiles_fish).number]];
        else            
            numbersForROIs=[MatFiles(MatFiles_fish(1)-1).number+1 [MatFiles(MatFiles_fish).number]];
        end
        IsInBrainRegion=[];
        IsInBrainRegion=ismember(ROI_fish,Hindbrain_Mask,'rows');
        PerBrainRegions(fish_nb).(regionName).ROIsCent=ROI_fish(IsInBrainRegion,:);
        temp=ZS(numbersForROIs(1):numbersForROIs(end),:);IsInBrainRegion=temp(IsInBrainRegion,:);
        PerBrainRegions(fish_nb).(regionName).ZS=IsInBrainRegion;
 end     
        
%Now to make it into one big pool of data
ZS_Brain=struct();
for i=1:length(ItiaList)
    regionName=ItiaList{i};        
    for j=1:length(PerBrainRegions)
        if j==1;
            ZS_Brain.(regionName)=PerBrainRegions(j).(regionName).ZS;            
        else
            ZS_Brain.(regionName)=vertcat(ZS_Brain.(regionName),PerBrainRegions(j).(regionName).ZS);            
        end
    end
end
clearvars i j Mask fish_nb IndexC IsInBrainRegion IsInEyes_temp regionName temp Fish_name

%----------------------------------------------------------------------
%_--------------------------------------------------------------------
%Analysis of ELO_ERO_Both

Stimuli=zeros(3,size(ZS,2));
start=42;
spike=[-0.104392135015146,1.69644104899772,5.13796058542217,8.27886020152244,10.3756715204800,11.8173714529814,12.2425184714093,10.8571417354877,8.80831829681196,6.91339112244670,5.46959264663869,4.30868766622567,3.42533619066766,2.75378443486879,2.18017250852183,1.72816235135824,1.32732537295463,1.00684435500268,0.730210038304555,0.530242444093118,0.362253250339685,0.227668255288566,0.0869242416152502,0.000718266708050853,-0.0828334873368325]';
for i = 1:3
    for j = 0:2
        Stimuli(i,start+j*120:start+j*120+length(spike)-1)=spike';
    end
    start=start+40;
end

progressbar;
for j=1:length(ItiaList)
    progressbar(j/length(ItiaList),[]);
    regionName=ItiaList{j};
    temp=ZS_Brain.(regionName);	
    coefficients=struct();
    rsquared=zeros(1,length(temp));
    for i=1:size(temp,1)
        progressbar([],i/size(temp,1));
        mdl=fitlm(Stimuli',temp(i,:));%,'interactions');
        coefficients(i).coef=mdl.Coefficients;
        rsquared(i)=mdl.Rsquared.Adjusted;
    end    
    LinReg.(regionName).coef=coefficients;
    LinReg.(regionName).rsquared=rsquared;    
end

% %Kmeans of all
% options = statset('UseParallel',1); 
% for j=1:11
%     regionName=ItiaList{j};
%     idx_temp=find(LinReg.(regionName).rsquared>0.1);
%     LinReg.(regionName).ZS_rsq=ZS_Brain.(regionName)(idx_temp,:);
%     if j==3 | j==4
%         [idxKmeans Cmap]=kmeans(LinReg.(regionName).ZS_rsq,5,'Options',options,'Distance','cityblock','Replicates',3,'MaxIter',1000,'Display','final');
%     else
%         [idxKmeans Cmap]=kmeans(LinReg.(regionName).ZS_rsq,10,'Options',options,'Distance','cityblock','Replicates',3,'MaxIter',1000,'Display','final');
%     end    
%     LinReg.(regionName).KmeansCenter=Cmap;
%     LinReg.(regionName).KmeansIdx=idxKmeans;
% end
% close all;
% %Get the good ones
% for j=1:length(ItiaList)
%     regionName=ItiaList{j};
%      if j==3 | j==4
%         [~,LinReg.(regionName).GoodBetas]=Test_Regress( LinReg.(regionName).KmeansCenter,Stimuli,LinReg.(regionName).KmeansIdx,0.2);
%     else
%         [~,LinReg.(regionName).GoodBetas]=Test_Regress( LinReg.(regionName).KmeansCenter,Stimuli,LinReg.(regionName).KmeansIdx,0.4);
%      end    
% end

%Pool data
region_nb=1;regionName=ItiaList{region_nb};
 idx_temp=find(LinReg.(regionName).rsquared>0.1);
    LinReg.(regionName).ZS_rsq=ZS_Brain.(regionName)(idx_temp,:);
ZS_pool=LinReg.(regionName).ZS_rsq;
for region_nb=2:length(ItiaList)
    regionName=ItiaList{region_nb};
    idx_temp=find(LinReg.(regionName).rsquared>0.1);
    LinReg.(regionName).ZS_rsq=ZS_Brain.(regionName)(idx_temp,:);
    ZS_pool=[ZS_pool;LinReg.(regionName).ZS_rsq];    
end
options = statset('UseParallel',1);
% [idxKmeans Cmap]=kmeans(ZS_pool,20,'Options',options,'Distance','cityblock','Replicates',5,'MaxIter',1000,'Display','final');
% [~,GoodBetas]=Test_Regress(Cmap,Stimuli,idxKmeans,0.3);

%Building the AVG across 3 presentations
Firststart=30;interstimulus=120;Nb_stimuli=3;
for j=1:length(ItiaList)
    regionName=ItiaList{j};
    ZS2=LinReg.(regionName).ZS_rsq;
    ZS_AVG2=zeros(size(ZS2,1),Nb_stimuli*41);
    parfor idx_ZS=1:size(ZS2,1)
        start=Firststart;
        AVG=[];
        for i=1:3
            AVG(i,:)=ZS2(idx_ZS,start:start+40);
            start=start+interstimulus;
        end
        AVG=mean(AVG,1);
        AVG=AVG-mean(AVG(1:5));%-min(AVG);
        temp=[];
        for j=2:3
            start=Firststart+40*(j-1);
            for i=1:3
                temp(i,:)=ZS2(idx_ZS,start:start+40);
                start=start+interstimulus;
            end
            temp=mean(temp,1);
            temp=temp-mean(temp(1:5));%-min(temp);
            AVG=[AVG temp];
        end
        ZS_AVG2(idx_ZS,:)=AVG;
    end
    LinReg.(regionName).ZS_AVG=ZS_AVG2;
end

region_nb=1;regionName=ItiaList{region_nb};
ZSAVG_pool=LinReg.(regionName).ZS_AVG;
for region_nb=2:length(ItiaList)
    regionName=ItiaList{region_nb};
    ZSAVG_pool=[ZSAVG_pool;LinReg.(regionName).ZS_AVG];    
end

[idxKmeans Cmap]=kmeans(ZS_pool,25,'Options',options,'Distance','cityblock','Replicates',3,'MaxIter',1000,'Display','final');
[~,GoodBetas]=Test_Regress(Cmap,Stimuli,idxKmeans,0.3);

%Correlate back to Kmeans to remove crap
Threshold=0.4;
KmeansIdx_select=idxKmeans;
for j=1:length(GoodBetas)
    idx_g=find(idxKmeans==GoodBetas(j));    
    corr_temp=zeros(size(idx_g));
    Mean_temp=mean(ZS_pool(idx_g,:),1);
    for k=1:length(idx_g)
        temp_corr=corrcoef(Mean_temp, ZS_pool(idx_g(k),:));
        corr_temp(k)=temp_corr(1,2);
    end
    KmeansIdx_select(idx_g(find(corr_temp<=Threshold)))=0;
end
clearvars i j k temp ans coeff coefficients explained i Hindbrain_Mask idx_temp ish_nb

figure;
counter=1;counter2=1;xplot=floor(sqrt(length(GoodBetas)));yplot=ceil(length(GoodBetas)/xplot);coloring={'r','g','y'};
for j=1:length(GoodBetas)  
    idx_temp=find(KmeansIdx_select==GoodBetas(j));
    subplot(xplot,yplot,j);plot(mean(ZS_pool(idx_temp,:),1));ylim([-3 3]);title(num2str(j));
end

ROIsPerBrain=struct();
for i=1:length(ItiaList)
    regionName=ItiaList{i};
    for j=1:length(PerBrainRegions)
        if j==1;
            ROIsPerBrain.(regionName).ROIs=PerBrainRegions(j).(regionName).ROIsCent;
            ROIsPerBrain.(regionName).Numbers=length(PerBrainRegions(j).(regionName).ROIsCent);
            temp=length(PerBrainRegions(j).(regionName).ROIsCent);
        else
            ROIsPerBrain.(regionName).ROIs=vertcat(ROIsPerBrain.(regionName).ROIs, PerBrainRegions(j).(regionName).ROIsCent); 
            temp=temp+length(PerBrainRegions(j).(regionName).ROIsCent);
            ROIsPerBrain.(regionName).Numbers=horzcat(ROIsPerBrain.(regionName).Numbers,temp);
        end
    end    
end
%Pool the ROIs
region_nb=1;regionName=ItiaList{region_nb};
ROI_temp=ROIsPerBrain.(regionName).ROIs;
ROI_temp=ROI_temp(find(LinReg.(regionName).rsquared>0.1),:);
ROI_pool=ROI_temp;
for region_nb=2:length(ItiaList)
    regionName=ItiaList{region_nb};
    ROI_temp=ROIsPerBrain.(regionName).ROIs;
    ROI_temp=ROI_temp(find(LinReg.(regionName).rsquared>0.1),:);
    ROI_pool=[ROI_pool;ROI_temp];
end

% For ERO_ELO Select and Merge
% GoodBetas_select=GoodBetas([2 3 5 7 8 10 11 12 13 17 19 20 21 24]);
% figure;
% counter=1;counter2=1;xplot=floor(sqrt(length(GoodBetas_select)));yplot=ceil(length(GoodBetas_select)/xplot);coloring={'r','g','y'};
% for j=1:length(GoodBetas_select)  
%     idx_temp=find(KmeansIdx_select==GoodBetas_select(j));
%     subplot(xplot,yplot,j);plot(mean(ZS_pool(idx_temp,:),1));ylim([-3 3]);xlim([1 420]); title(num2str(j));
% end
GoodBetas_select=GoodBetas;

KmeansIdx_merge=KmeansIdx_select;
idx_temp=ismember(KmeansIdx_select,GoodBetas_select([8 5 7 20 21]));
KmeansIdx_merge(idx_temp)=GoodBetas_select(8);
idx_temp=ismember(KmeansIdx_select,GoodBetas_select([17 24 6 10 12 19]));
KmeansIdx_merge(idx_temp)=GoodBetas_select(17);
idx_temp=ismember(KmeansIdx_select,GoodBetas_select([2 16]));
KmeansIdx_merge(idx_temp)=GoodBetas_select(2);
idx_temp=ismember(KmeansIdx_select,GoodBetas_select([11 3]));
KmeansIdx_merge(idx_temp)=GoodBetas_select(11);

GoodBetas_merge=GoodBetas([8 17 2 11 13]);

%Correlate back to Kmeans to remove crap
Threshold=0.4;
for j=1:length(GoodBetas_merge)
    idx_g=find(KmeansIdx_select==GoodBetas_merge(j));    
    Mean_temp=mean(ZS_pool(idx_g,:),1);    
    idx_g=find(KmeansIdx_merge==GoodBetas_merge(j));    
    corr_temp=zeros(size(idx_g));
    for k=1:length(idx_g)
        temp_corr=corrcoef(Mean_temp, ZS_pool(idx_g(k),:));
        corr_temp(k)=temp_corr(1,2);
    end
    KmeansIdx_merge(idx_g(find(corr_temp<=Threshold)))=0;
end
clearvars i j k temp ans coeff coefficients explained i Hindbrain_Mask idx_temp ish_nb

figure;
counter=1;counter2=1;xplot=floor(sqrt(length(GoodBetas_merge)));yplot=ceil(length(GoodBetas_merge)/xplot);coloring={'r','g','y'};
for j=1:length(GoodBetas_merge)  
    idx_temp=find(KmeansIdx_merge==GoodBetas_merge(j));
    subplot(xplot,yplot,j);plot(mean(ZS_pool(idx_temp,:),1));ylim([-3 3]);xlim([1 420]);
end

colors=[0.14 1 0.14; 0.7 0.4 1; 0 0.6 0.6];
x = linspace(0.25,123/4,123);
for i=1:3
    Fighandle=figure;
    set(Fighandle, 'Position', [100, 100, 250, 250]);
    idx_temp=find(KmeansIdx_merge==GoodBetas_merge(i));
    temp=ZSAVG_pool(idx_temp,:);    
    for k=0:2           
            hold on;
            meanToPlot=mean(temp(:,4+(k*41):40+(k*41)),1);
            plot(x(4+(k*40):40+(k*40)),meanToPlot-mean(meanToPlot(1:4)),'color',colors(1,:),'LineWidth',3);axis([0 30 -3 3]);                       
            rectangle('EdgeColor','none','FaceColor',[0.25, 0.25, 0.25, 0.4],'Position',[x(9+(k*40)) -3 1 7]);
    end
	hold off;    
    print(Fighandle,strcat('_Cluster_ELO_ERO-',num2str(i)),'-dsvg','-r0');
end
i=4;
Fighandle=figure;
set(Fighandle, 'Position', [100, 100, 250, 250]);
idx_temp=find(KmeansIdx_merge==GoodBetas_merge(i));
temp=ZSAVG_pool(idx_temp,:);
for k=0:2
    hold on;
    meanToPlot=mean(temp(:,4+(k*41):40+(k*41)),1);
    plot(x(4+(k*40):40+(k*40)),meanToPlot-mean(meanToPlot(1:4)),'color',colors(1,:),'LineWidth',3);axis([0 30 -3 3]);
    rectangle('EdgeColor','none','FaceColor',[0.25, 0.25, 0.25, 0.4],'Position',[x(9+(k*40)) -3 1 7]);
end
i=5;
idx_temp=find(KmeansIdx_merge==GoodBetas_merge(i));
temp=ZSAVG_pool(idx_temp,:);
for k=0:2
    hold on;
    meanToPlot=mean(temp(:,4+(k*41):40+(k*41)),1);
    plot(x(4+(k*40):40+(k*40)),meanToPlot-mean(meanToPlot(1:4)),'color',colors(2,:),'LineWidth',3);axis([0 30 -3 3]);
    %rectangle('EdgeColor','none','FaceColor',[0.25, 0.25, 0.25, 0.4],'Position',[x(9+(k*40)) -3 1 7]);
end
% i=6;
% idx_temp=find(KmeansIdx_merge==GoodBetas_merge(i));
% temp=ZSAVG_pool(idx_temp,:);
% for k=0:2
%     hold on;
%     meanToPlot=mean(temp(:,4+(k*41):40+(k*41)),1);
%     plot(x(4+(k*40):40+(k*40)),meanToPlot-mean(meanToPlot(1:4)),'color',colors(3,:),'LineWidth',3);axis([0 30 -3 3]);
%     %rectangle('EdgeColor','none','FaceColor',[0.25, 0.25, 0.25, 0.4],'Position',[x(9+(k*40)) -3 1 7]);
% end

print(Fighandle,strcat('_Cluster-','4_5_6_ELO-ERO'),'-dsvg','-r0');


for i=1:length(GoodBetas_merge)    
    idx_temp=find(KmeansIdx_merge==GoodBetas_merge(i));
    CSV_temp=ROI_pool(idx_temp,:);
    CSV_temp(:,3)=CSV_temp(:,3)*2;
    CSV_temp(:,4)=1;
    filename=strcat('_ROIs_coord_clust_ELOERO',num2str(i),'.csv');
    csvwrite(filename,CSV_temp);
end

%Coeff pooled

coefficients_pool=[];
for i=1:length(ItiaList)
    regionName=ItiaList{i};
    coefficients={};
    for idx=1:length(LinReg.(regionName).coef)
        coef=[LinReg.(regionName).coef(idx).coef];
        temp=coef.Properties.RowNames;temp=regexp(temp,'x(\d+)','tokens');
        if ~isempty(temp)
            %temp=[temp{:}];temp=[temp{:}];temp=[temp{:}];%temp=str2num(temp);
            for coef_idx=2:height(coef)
                %             if coef.pValue(coef_idx)<0.05
                coefficients{idx,str2num(temp{coef_idx}{1}{1})}=coef.Estimate(coef_idx);
                %             end
            end
        end
    end
    idxempty=cellfun('isempty',coefficients);
    coefficients(idxempty)={0};
    clearvars idxempty idx coef_idx coef temp
    coefficients=cell2mat(coefficients);
    coefficients=coefficients(find(LinReg.(regionName).rsquared>0.1),:);
    if isempty(coefficients_pool)
        coefficients_pool=coefficients;
    else
        coefficients_pool=[coefficients_pool;coefficients];
    end
end
clearvars i j k temp ans coeff coefficients explained i Hindbrain_Mask idx_temp ish_nb
mean_coef=mean(coefficients_pool,1);
std_coef=std(coefficients_pool,1,1);    
idx_weird=find(coefficients_pool(:,1)>mean_coef(1)+std_coef(1) & coefficients_pool(:,2)<mean_coef(2)-std_coef(2));
idx_weird2=find(coefficients_pool(:,1)<mean_coef(1)-std_coef(1) & coefficients_pool(:,2)>mean_coef(2)+std_coef(2));

x = linspace(0.25,123/4,123);
Fighandle=figure;
set(Fighandle, 'Position', [100, 100, 250, 250]);
temp=ZSAVG_pool(idx_weird,:);
for k=0:2
    hold on;
    meanToPlot=mean(temp(:,4+(k*41):40+(k*41)),1);
    plot(x(4+(k*40):40+(k*40)),meanToPlot-mean(meanToPlot(1:4)),'color',colors(2,:),'LineWidth',3);axis([0 30 -3 3]);
    rectangle('EdgeColor','none','FaceColor',[0.25, 0.25, 0.25, 0.4],'Position',[x(9+(k*40)) -3 1 7]);
end
temp=ZSAVG_pool(idx_weird2,:);
for k=0:2
    hold on;
    meanToPlot=mean(temp(:,4+(k*41):40+(k*41)),1);
    plot(x(4+(k*40):40+(k*40)),meanToPlot-mean(meanToPlot(1:4)),'color',colors(1,:),'LineWidth',3);axis([0 30 -3 3]);
    %rectangle('EdgeColor','none','FaceColor',[0.25, 0.25, 0.25, 0.4],'Position',[x(9+(k*40)) -3 1 7]);
end
print(Fighandle,strcat('_Cluster-ELOERO','weird'),'-dsvg','-r0');

CSV_temp=ROI_pool(idx_weird,:);
CSV_temp(:,4)=2;
CSV_temp2=ROI_pool(idx_weird2,:);
CSV_temp2(:,4)=1;
CSV_temp=[CSV_temp;CSV_temp2];
CSV_temp(:,3)=CSV_temp(:,3)*2;
filename=strcat('_ROIs_coord_clust_','weird','.csv');
csvwrite(filename,CSV_temp);
clearvars i j k counter counter temp idx_temp idx_temp2 idx_weird idx_weird2 CSV_temp list_merge corr_temp CSV_temp2
clearvars counter2 Fighandle filename Firststart GN GCaMP6 GoodBet_temp idx_g idxStart interstimulus mean_coef std_coef Mean_temp
clearvars meanToPlot ModelCmap Nb_stimuli options region_nb regionName ROI_temp score Stim temp_corr temp_g temp_ZS Threshold x y xplot yplot

%adding pretectum

%Correlate back to Kmeans to remove crap
Threshold=0.4;
lengthZSpool_withoutpretectum=100797;
KmeansIdx_merge_pretec=KmeansIdx_merge;
Clusters_mean=zeros(length(GoodBetas_merge),size(ZS_pool,2));
for i=1:length(GoodBetas_merge)
    Clusters_mean(i,:)=mean(ZS_pool(find(KmeansIdx_merge==GoodBetas_merge(i)),:),1);
end

Threshold=0.4;
for j=(lengthZSpool_withoutpretectum+1):length(ZS_pool)    
    corr_temp=zeros(1,length(GoodBetas_merge));    
    for i=1:length(GoodBetas_merge)
        temp_corr=corrcoef(Clusters_mean(i,:), ZS_pool(j,:));   
        corr_temp(i)=temp_corr(1,2);
    end
    [max_temp,idx_temp]=max(corr_temp);
    if max_temp>Threshold
        KmeansIdx_merge_pretec(j)=GoodBetas_merge(idx_temp);
    else        
        KmeansIdx_merge_pretec(j)=0;
    end
end
clearvars i j k temp ans coeff coefficients explained i Hindbrain_Mask idx_temp ish_nb


figure;
counter=1;counter2=1;xplot=floor(sqrt(length(GoodBetas_merge)));yplot=ceil(length(GoodBetas_merge)/xplot);coloring={'r','g','y'};
for j=1:length(GoodBetas_merge)  
    idx_temp=find(KmeansIdx_merge_pretec==GoodBetas_merge(j));
    subplot(xplot,yplot,j);plot(mean(ZS_pool(idx_temp,:),1));ylim([-3 3]);xlim([1 420]);
end

colors=[0.14 1 0.14; 0.7 0.4 1; 0 0.6 0.6];
x = linspace(0.25,123/4,123);
for i=1:3
    Fighandle=figure;
    set(Fighandle, 'Position', [100, 100, 250, 250]);
    idx_temp=find(KmeansIdx_merge_pretec==GoodBetas_merge(i));
    temp=ZSAVG_pool(idx_temp,:);    
    for k=0:2           
            hold on;
            meanToPlot=mean(temp(:,4+(k*41):40+(k*41)),1);
            plot(x(4+(k*40):40+(k*40)),meanToPlot-mean(meanToPlot(1:4)),'color',colors(1,:),'LineWidth',3);axis([0 30 -3 3]);                       
            rectangle('EdgeColor','none','FaceColor',[0.25, 0.25, 0.25, 0.4],'Position',[x(9+(k*40)) -3 1 7]);
    end
	hold off;    
    print(Fighandle,strcat('_Cluster_ELO_ERO-',num2str(i)),'-dsvg','-r0');
end
i=4;
Fighandle=figure;
set(Fighandle, 'Position', [100, 100, 250, 250]);
idx_temp=find(KmeansIdx_merge_pretec==GoodBetas_merge(i));
temp=ZSAVG_pool(idx_temp,:);
for k=0:2
    hold on;
    meanToPlot=mean(temp(:,4+(k*41):40+(k*41)),1);
    plot(x(4+(k*40):40+(k*40)),meanToPlot-mean(meanToPlot(1:4)),'color',colors(1,:),'LineWidth',3);axis([0 30 -3 3]);
    rectangle('EdgeColor','none','FaceColor',[0.25, 0.25, 0.25, 0.4],'Position',[x(9+(k*40)) -3 1 7]);
end
i=5;
idx_temp=find(KmeansIdx_merge_pretec==GoodBetas_merge(i));
temp=ZSAVG_pool(idx_temp,:);
for k=0:2
    hold on;
    meanToPlot=mean(temp(:,4+(k*41):40+(k*41)),1);
    plot(x(4+(k*40):40+(k*40)),meanToPlot-mean(meanToPlot(1:4)),'color',colors(2,:),'LineWidth',3);axis([0 30 -3 3]);
    %rectangle('EdgeColor','none','FaceColor',[0.25, 0.25, 0.25, 0.4],'Position',[x(9+(k*40)) -3 1 7]);
end

print(Fighandle,strcat('_Cluster-','4_5_6_ELO-ERO'),'-dsvg','-r0');


for i=1:length(GoodBetas_merge)    
    idx_temp=find(KmeansIdx_merge_pretec==GoodBetas_merge(i));
    CSV_temp=ROI_pool(idx_temp,:);
    CSV_temp(:,3)=CSV_temp(:,3)*2;
    CSV_temp(:,4)=1;
    filename=strcat('_ROIs_coord_clust_ELOERO',num2str(i),'.csv');
    csvwrite(filename,CSV_temp);
end

%---------------------------------------------
%Remove ROIs from outside the brain
Zbrain_AllMask=vertcat(Zbrain_Masks{[1:1:77 79:1:294],3});
Zbrain_AllMask=unique(Zbrain_AllMask,'rows');
%Removing the eyes
idx_brain=ismember(ROI_rot_correct,Zbrain_AllMask,'rows');

