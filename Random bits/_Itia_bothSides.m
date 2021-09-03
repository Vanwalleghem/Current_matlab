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
% %     end
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

Stimuli=zeros(4,size(ZS,2));
start=43;
spike=[0,1.69644104899772,5.13796058542217,8.27886020152244,10.3756715204800,12.2425184714093,10.8571417354877,8.80831829681196,6.91339112244670,5.46959264663869,4.30868766622567,3.42533619066766,2.75378443486879,2.18017250852183,1.72816235135824,1.32732537295463,1.00684435500268,0.730210038304555,0.530242444093118,0.362253250339685,0.0869242416152502,0.000718266708050853]';
for i = 1:3
    for j = 0:2
        Stimuli(i,start+j*120:start+j*120+length(spike)-1)=spike';
    end
    start=start+20;
end
start=67;i=4;
for j = 0:2
    Stimuli(i,start+j*60:start+j*60+length(spike)-1)=spike';
end
clearvars i j k spike start

 load('ELO_IRO.mat','Stimuli')
 Stimuli(4,:)=Stimuli(2,:);
 Stimuli(4,5:end)=Stimuli(2,1:end-4);

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

CSV_Files=dir('_2Warped*ELO_IRO*.csv');
%CSV_Files=dir('_Warped*.csv');
ROIs=struct();
for i=1:length(CSV_Files);
    temp=csvread(CSV_Files(i).name,1);
    %Fishname=regexp(CSV_Files(i).name,'Fish(\d+_\w+).csv','tokens');Fishname=Fishname{1}{1};
    Fishname=regexp(CSV_Files(i).name,'FishResized(\d+_\w+).csv','tokens');Fishname=Fishname{1}{1};
    %Fishname=regexp(CSV_Files(i).name,'Betas(\d+_\w+).csv','tokens');Fishname=Fishname{1}{1};
    ROIs(i).name=Fishname;    
    ROIs(i).coord=temp(:,1:3);
    ROIs(i).idx=temp(:,5);
    %ROIs(i).cluster=temp(:,6);
end
clearvars i temp CSV_Files Fishname

load('__MultiBrainRegions_final.mat','MON_Masks')
ItiaList={'Thalamus','Cerebellum','NucMLF','Semicircularis','Telencephalon','Tectum','Longitudinalis','Tegmentum','Habenula','Hindbrain','MON'};
i=11;
regionName=ItiaList{i};
Mask=MON_Masks{10,3};
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


i=11;
regionName=ItiaList{i};
for j=1:length(PerBrainRegions)
    if j==1;
        ZS_Brain.(regionName)=PerBrainRegions(j).(regionName).ZS;
    else
        ZS_Brain.(regionName)=vertcat(ZS_Brain.(regionName),PerBrainRegions(j).(regionName).ZS);
    end
end

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
clearvars i j k ans Cmap coef coefficients counter Fish_name fish_nb GCaMP6 GoodBetas_nuc idxStart idx_Fish idx_temp idxKmeans IndexC image_temp IsInBrainRegion
clearvars IsInBrainRegion_good Mask MatFiles_fish mdl numbersForROIs_good numbersForROIs regionName ROI_fish ROI_fish_good ROI_name rsquared  temp

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

%Kmeans of all
options = statset('UseParallel',1); 
for j=1:11
    regionName=ItiaList{j};
    idx_temp=find(LinReg.(regionName).rsquared>0.1);
    LinReg.(regionName).ZS_rsq=ZS_Brain.(regionName)(idx_temp,:);
    if j==3 | j==4
        [idxKmeans Cmap]=kmeans(LinReg.(regionName).ZS_rsq,5,'Options',options,'Distance','cityblock','Replicates',3,'MaxIter',1000,'Display','final');
    else
        [idxKmeans Cmap]=kmeans(LinReg.(regionName).ZS_rsq,10,'Options',options,'Distance','cityblock','Replicates',3,'MaxIter',1000,'Display','final');
    end    
    LinReg.(regionName).KmeansCenter=Cmap;
    LinReg.(regionName).KmeansIdx=idxKmeans;
end
close all;
%Get the good ones
for j=1:length(ItiaList)
    regionName=ItiaList{j};
     if j==3 | j==4
        [~,LinReg.(regionName).GoodBetas]=Test_Regress( LinReg.(regionName).KmeansCenter,Stimuli,LinReg.(regionName).KmeansIdx,0.2);
    else
        [~,LinReg.(regionName).GoodBetas]=Test_Regress( LinReg.(regionName).KmeansCenter,Stimuli,LinReg.(regionName).KmeansIdx,0.4);
     end    
end

% %Kmeans of nucMLF (Noisy)
% options = statset('UseParallel',1); 
% for j=3:4
%     regionName=ItiaList{j};
%     idx_temp=find(LinReg.(regionName).rsquared>0.1);
%     LinReg.(regionName).ZS_rsq=ZS_Brain.(regionName)(idx_temp,:);
%     [idxKmeans Cmap]=kmeans(LinReg.(regionName).ZS_rsq,4,'Options',options,'Distance','cityblock','Replicates',3,'MaxIter',1000,'Display','final');
%     LinReg.(regionName).KmeansCenter=Cmap;
%     LinReg.(regionName).KmeansIdx=idxKmeans;
% end
% 
% %Get the good ones
% for j=3:4
%     regionName=ItiaList{j};
%     [~,LinReg.(regionName).GoodBetas]=Test_Regress( LinReg.(regionName).KmeansCenter,Stimuli,LinReg.(regionName).KmeansIdx,0.2);
% end

%Correlate back to Kmeans to remove crap
Threshold=0.4;
for i=1:length(ItiaList)
    regionName=ItiaList{i};
    temp=LinReg.(regionName).ZS_rsq;
    idx_temp=LinReg.(regionName).KmeansIdx;
    GoodBet_temp=LinReg.(regionName).GoodBetas;
    for j=1:length(GoodBet_temp)
        idx_g=find(idx_temp==GoodBet_temp(j));
        temp_g=temp(idx_g,:);
        corr_temp=zeros(size(idx_g));
        Mean_temp=LinReg.(regionName).KmeansCenter(GoodBet_temp(j),:);
        for k=1:length(idx_g)        
            temp_corr=corrcoef(Mean_temp, temp_g(k,:));
            corr_temp(k)=temp_corr(1,2);
        end        
        idx_temp(idx_g(find(corr_temp<=Threshold)))=0;
    end
    LinReg.(regionName).KmeansIdx_select=idx_temp;    
end
clearvars i j k temp ans Cmap coeff coefficients explained i Hindbrain_Mask idx_temp idxKmeans ish_nb 

Threshold=0.2;
for i=3:4
    regionName=ItiaList{i};
    temp=LinReg.(regionName).ZS_rsq;
    idx_temp=LinReg.(regionName).KmeansIdx;
    GoodBet_temp=LinReg.(regionName).GoodBetas;
    for j=1:length(GoodBet_temp)
        idx_g=find(idx_temp==GoodBet_temp(j));
        temp_g=temp(idx_g,:);
        corr_temp=zeros(size(idx_g));
        Mean_temp=LinReg.(regionName).KmeansCenter(GoodBet_temp(j),:);
        for k=1:length(idx_g)        
            temp_corr=corrcoef(Mean_temp, temp_g(k,:));
            corr_temp(k)=temp_corr(1,2);
        end        
        idx_temp(idx_g(find(corr_temp<=Threshold)))=0;
    end
    LinReg.(regionName).KmeansIdx_select=idx_temp;    
end
clearvars i j k temp ans Cmap coeff coefficients explained i Hindbrain_Mask idx_temp idxKmeans ish_nb 

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

%Kmeans of AVG
options = statset('UseParallel',1); 
for j=1:length(ItiaList)
    regionName=ItiaList{j};    
    if j==3
        [idxKmeans Cmap]=kmeans(LinReg.(regionName).ZS_AVG,5,'Distance','cityblock','Options',options,'Replicates',3,'MaxIter',1000,'Display','final');    
    else
        [idxKmeans Cmap]=kmeans(LinReg.(regionName).ZS_AVG,10,'Distance','cityblock','Options',options,'Replicates',3,'MaxIter',1000,'Display','final');
    end
    LinReg.(regionName).KmeansCenter_AVG=Cmap;
    LinReg.(regionName).KmeansIdx_AVG=idxKmeans;
end

%Define the AVG stimuli
Stimuli_AVG=zeros(3,123);
GCaMP6=[0,1.69644104899772,5.13796058542217,8.27886020152244,10.3756715204800,11.8173714529814,12.2425184714093,10.8571417354877,8.80831829681196,6.91339112244670,5.46959264663869,4.30868766622567,3.42533619066766,2.75378443486879,2.18017250852183,1.72816235135824,1.32732537295463,1.00684435500268,0.730210038304555,0.530242444093118,0.362253250339685,0.227668255288566,0.0869242416152502,0.000718266708050853,-0.0828334873368325]';
idxStart=10;
for i=1:3    
    Stimuli_AVG(i,(idxStart+(i-1)*41):(idxStart+(i-1)*41)+size(GCaMP6,1)-1)=GCaMP6;
end

%Get the good ones
for j=1:length(ItiaList)
    regionName=ItiaList{j};
    [~,LinReg.(regionName).GoodBetas_AVG]=Test_Regress( LinReg.(regionName).KmeansCenter_AVG,Stimuli_AVG,LinReg.(regionName).KmeansIdx_AVG,0.3);
end

figure;
counter=1;xplot=length(ItiaList);yplot=10;counter2=1;
for i=1:length(ItiaList)
    regionName=ItiaList{i};
    idx_temp=LinReg.(regionName).KmeansIdx_select;
    GoodBet_temp=LinReg.(regionName).GoodBetas;
    counter=counter2;
    for j=1:length(GoodBet_temp)
        idx_temp2=find(idx_temp==GoodBet_temp(j));
        subplot(xplot,yplot,counter+(j-1));plot(mean(LinReg.(regionName).ZS_AVG(idx_temp2,:),1));title(num2str(length(idx_temp2)));axis([0 120 -3 3]);
    end
    counter2=counter2+yplot;
end

i=1;
regionName=ItiaList{i};
temp_ZS=LinReg.(regionName).ZS_rsq;
idx_temp=LinReg.(regionName).KmeansIdx_select;
GoodBet_temp=LinReg.(regionName).GoodBetas;
idx_g=ismember(idx_temp,GoodBet_temp([1 235]));
temp_g=temp_ZS(idx_g,:);
corr_temp=zeros(size(idx_g));
for k=1:length(idx_g)
    temp=corrcoef(LinReg.(regionName).KmeansCenter(GoodBet_temp(2),:), temp_g(k,:));
    corr_temp(k)=temp(1,2);
end
idx_temp(idx_g(find(corr_temp<=Threshold)))=0;
idx_temp(idx_g(find(corr_temp>Threshold)))=GoodBet_temp(2);
LinReg.(regionName).KmeansIdx_select_merge=idx_temp;

%----------------------------------------------------------------------------------------------------------------------------------------------------------------
%----------------------------------------------------------------------------------------------------------------------------------------------------------------

%Kmeans of all
options = statset('UseParallel',1); 
for j=1:11
    regionName=ItiaList{j};
    if j==3 | j==4
        [idxKmeans Cmap]=kmeans(LinReg.(regionName).ZS_rsq,5,'Options',options,'Distance','cityblock','Replicates',3,'MaxIter',1000,'Display','final');
    else
        [idxKmeans Cmap]=kmeans(LinReg.(regionName).ZS_rsq,5,'Options',options,'Distance','cityblock','Replicates',3,'MaxIter',1000,'Display','final');
    end    
    LinReg.(regionName).KmeansCenter2=Cmap;
    LinReg.(regionName).KmeansIdx2=idxKmeans;
end
close all;

%Get the good ones
for j=1:length(ItiaList)
    regionName=ItiaList{j};
     if j==3 | j==4
        [~,LinReg.(regionName).GoodBetas2]=Test_Regress( LinReg.(regionName).KmeansCenter2,Stimuli,LinReg.(regionName).KmeansIdx2,0.2);
    else
        [~,LinReg.(regionName).GoodBetas2]=Test_Regress( LinReg.(regionName).KmeansCenter2,Stimuli,LinReg.(regionName).KmeansIdx2,0.4);
     end    
end

%Correlate back to Kmeans to remove crap
Threshold=0.4;
for i=1:length(ItiaList)
    regionName=ItiaList{i};
    temp=LinReg.(regionName).ZS_rsq;
    idx_temp=LinReg.(regionName).KmeansIdx2;
    GoodBet_temp=LinReg.(regionName).GoodBetas2;
    for j=1:length(GoodBet_temp)
        idx_g=find(idx_temp==GoodBet_temp(j));
        temp_g=temp(idx_g,:);
        corr_temp=zeros(size(idx_g));
        Mean_temp=LinReg.(regionName).KmeansCenter2(GoodBet_temp(j),:);
        for k=1:length(idx_g)        
            temp_corr=corrcoef(Mean_temp, temp_g(k,:));
            corr_temp(k)=temp_corr(1,2);
        end        
        idx_temp(idx_g(find(corr_temp<=Threshold)))=0;
    end
    LinReg.(regionName).KmeansIdx_select2=idx_temp;    
end
clearvars i j k temp ans Cmap coeff coefficients explained i Hindbrain_Mask idx_temp idxKmeans ish_nb 

Threshold=0.2;
for i=3:4
    regionName=ItiaList{i};
    temp=LinReg.(regionName).ZS_rsq;
    idx_temp=LinReg.(regionName).KmeansIdx2;
    GoodBet_temp=LinReg.(regionName).GoodBetas2;
    for j=1:length(GoodBet_temp)
        idx_g=find(idx_temp==GoodBet_temp(j));
        temp_g=temp(idx_g,:);
        corr_temp=zeros(size(idx_g));
        Mean_temp=LinReg.(regionName).KmeansCenter2(GoodBet_temp(j),:);
        for k=1:length(idx_g)        
            temp_corr=corrcoef(Mean_temp, temp_g(k,:));
            corr_temp(k)=temp_corr(1,2);
        end        
        idx_temp(idx_g(find(corr_temp<=Threshold)))=0;
    end
    LinReg.(regionName).KmeansIdx_select2=idx_temp;    
end
clearvars i j k temp ans Cmap coeff coefficients explained i Hindbrain_Mask idx_temp idxKmeans ish_nb 

figure;
counter=1;xplot=length(ItiaList);yplot=5;counter2=1;
for i=1:length(ItiaList)
    regionName=ItiaList{i};
    idx_temp=LinReg.(regionName).KmeansIdx_select2;
    GoodBet_temp=LinReg.(regionName).GoodBetas2;
    counter=counter2;
    for j=1:length(GoodBet_temp)
        idx_temp2=find(idx_temp==GoodBet_temp(j));
        subplot(xplot,yplot,counter+(j-1));plot(mean(LinReg.(regionName).ZS_AVG(idx_temp2,:),1));title(num2str(length(idx_temp2)));axis([0 120 -3 3]);
    end
    counter2=counter2+yplot;
end


%---------------------------------------------------------------------------------------------
%Coeff
progressbar;
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
    mean_coef=mean(coefficients,1);
    std_coef=std(coefficients,1,1);
    
    idx_weird=find(coefficients(:,1)>mean_coef(1)+std_coef(1) & coefficients(:,2)<mean_coef(2)-std_coef(2));
    idx_weird2=find(coefficients(:,1)<mean_coef(1)-std_coef(1) & coefficients(:,2)>mean_coef(2)+std_coef(2));
    %figure;plot(mean(ZS_Brain.(regionName)(idx_weird,:),1));hold on;plot(mean(ZS_Brain.(regionName)(idx_weird2,:),1));
    %figure;imagesc(ZS_Brain.(regionName)(idx_weird,:),[-3 3]);colormap hot;figure;imagesc(ZS_Brain.(regionName)(idx_weird2,:),[-3 3]);colormap hot;
    LinReg.(regionName).idx_rightInhib=idx_weird;
    LinReg.(regionName).idx_leftInhib=idx_weird2;
    progressbar(i/length(ItiaList));
end

figure;
counter=1;xplot=floor(sqrt(length(ItiaList)));yplot=ceil(length(ItiaList)/xplot);
for i=1:length(ItiaList)
    regionName=ItiaList{i};
    idx_weird=LinReg.(regionName).idx_rightInhib;
    idx_weird2=LinReg.(regionName).idx_leftInhib;
    subplot(xplot,yplot,i);
    plot(mean(ZS_Brain.(regionName)(idx_weird,:),1));hold on;plot(mean(ZS_Brain.(regionName)(idx_weird2,:),1));title(strcat(regionName,' - left :',num2str(length(idx_weird)),' right :',num2str(length(idx_weird2))));
end

i=2;
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
mean_coef=mean(coefficients,1);
std_coef=std(coefficients,1,1);

idx_weird=find(coefficients(:,1)>mean_coef(1)+std_coef(1) & coefficients(:,2)<mean_coef(2)-std_coef(2));
idx_weird2=find(coefficients(:,1)<mean_coef(1)-std_coef(1) & coefficients(:,2)>mean_coef(2)+std_coef(2));
figure;plot(mean(ZS_Brain.(regionName)(idx_weird,:),1));hold on;plot(mean(ZS_Brain.(regionName)(idx_weird2,:),1));
figure;imagesc(ZS_Brain.(regionName)(idx_weird,:),[-3 3]);colormap hot;figure;imagesc(ZS_Brain.(regionName)(idx_weird2,:),[-3 3]);colormap hot;

figure;histogram(LinReg.(regionName).rsquared([idx_weird ; idx_weird2]));
figure;histogram(idx_weird);hold on;histogram(idx_weird2);


%Merge clusters together
Threshold=0.4;
i=1;
regionName=ItiaList{i};
temp_ZS=LinReg.(regionName).ZS_rsq;
idx_temp=LinReg.(regionName).KmeansIdx_select2;
GoodBet_temp=LinReg.(regionName).GoodBetas2;
list_merge=[3 4];
idx_g=ismember(idx_temp,GoodBet_temp(list_merge));
idx_g=find(idx_g==1);
idx_temp(idx_g)=GoodBet_temp(list_merge(1));
LinReg.(regionName).KmeansIdx_select=idx_temp;

Threshold=0.4;
i=2;
regionName=ItiaList{i};
temp_ZS=LinReg.(regionName).ZS_rsq;
idx_temp=LinReg.(regionName).KmeansIdx_select2;
GoodBet_temp=LinReg.(regionName).GoodBetas2;
list_merge=[2 3 4];
idx_g=ismember(idx_temp,GoodBet_temp(list_merge));
idx_g=find(idx_g==1);
idx_temp(idx_g)=GoodBet_temp(list_merge(1));
LinReg.(regionName).KmeansIdx_select=idx_temp;

Threshold=0.2;
i=3;
regionName=ItiaList{i};
temp_ZS=LinReg.(regionName).ZS_rsq;
idx_temp=LinReg.(regionName).KmeansIdx_select2;
GoodBet_temp=LinReg.(regionName).GoodBetas2;
list_merge=[2 4];
idx_g=ismember(idx_temp,GoodBet_temp(list_merge));
idx_g=find(idx_g==1);
idx_temp(idx_g)=GoodBet_temp(list_merge(1));
LinReg.(regionName).KmeansIdx_select=idx_temp;

Threshold=0.4;
i=4;
regionName=ItiaList{i};
temp_ZS=LinReg.(regionName).ZS_rsq;
idx_temp=LinReg.(regionName).KmeansIdx_select2;
GoodBet_temp=LinReg.(regionName).GoodBetas2;
list_merge=[2 3];
idx_g=ismember(idx_temp,GoodBet_temp(list_merge));
idx_g=find(idx_g==1);
idx_temp(idx_g)=GoodBet_temp(list_merge(1));
LinReg.(regionName).KmeansIdx_select=idx_temp;

Threshold=0.4;
i=5;
regionName=ItiaList{i};
temp_ZS=LinReg.(regionName).ZS_rsq;
idx_temp=LinReg.(regionName).KmeansIdx_select2;
GoodBet_temp=LinReg.(regionName).GoodBetas2;
list_merge=[1 4];
idx_g=ismember(idx_temp,GoodBet_temp(list_merge));
idx_g=find(idx_g==1);
idx_temp(idx_g)=GoodBet_temp(list_merge(1));
LinReg.(regionName).KmeansIdx_select=idx_temp;

Threshold=0.4;
i=5;
regionName=ItiaList{i};
temp_ZS=LinReg.(regionName).ZS_rsq;
idx_temp=LinReg.(regionName).KmeansIdx_select2;
GoodBet_temp=LinReg.(regionName).GoodBetas2;
list_merge=[1 4];
idx_g=ismember(idx_temp,GoodBet_temp(list_merge));
idx_g=find(idx_g==1);
idx_temp(idx_g)=GoodBet_temp(list_merge(1));
LinReg.(regionName).KmeansIdx_select=idx_temp;

Threshold=0.4;
i=6;
regionName=ItiaList{i};
temp_ZS=LinReg.(regionName).ZS_rsq;
idx_temp=LinReg.(regionName).KmeansIdx_select2;
GoodBet_temp=LinReg.(regionName).GoodBetas2;
list_merge=[1 4];
idx_g=ismember(idx_temp,GoodBet_temp(list_merge));
idx_g=find(idx_g==1);
idx_temp(idx_g)=GoodBet_temp(list_merge(1));
LinReg.(regionName).KmeansIdx_select=idx_temp;

Threshold=0.4;
i=7;
regionName=ItiaList{i};
temp_ZS=LinReg.(regionName).ZS_rsq;
idx_temp=LinReg.(regionName).KmeansIdx_select2;
GoodBet_temp=LinReg.(regionName).GoodBetas2;
list_merge=[1 3 4];
idx_g=ismember(idx_temp,GoodBet_temp(list_merge));
idx_g=find(idx_g==1);
idx_temp(idx_g)=GoodBet_temp(list_merge(1));
LinReg.(regionName).KmeansIdx_select=idx_temp;

Threshold=0.4;
i=9;
regionName=ItiaList{i};
temp_ZS=LinReg.(regionName).ZS_rsq;
idx_temp=LinReg.(regionName).KmeansIdx_select2;
GoodBet_temp=LinReg.(regionName).GoodBetas2;
list_merge=[2 3 5];
idx_g=ismember(idx_temp,GoodBet_temp(list_merge));
idx_g=find(idx_g==1);
idx_temp(idx_g)=GoodBet_temp(list_merge(1));
LinReg.(regionName).KmeansIdx_select=idx_temp;

Threshold=0.4;
i=11;
regionName=ItiaList{i};
temp_ZS=LinReg.(regionName).ZS_rsq;
idx_temp=LinReg.(regionName).KmeansIdx_select2;
GoodBet_temp=LinReg.(regionName).GoodBetas2;
list_merge=[4 2];
idx_g=ismember(idx_temp,GoodBet_temp(list_merge));
idx_g=find(idx_g==1)
idx_temp(idx_g)=GoodBet_temp(list_merge(1));
LinReg.(regionName).KmeansIdx_select=idx_temp;

figure;
counter=1;xplot=length(ItiaList);yplot=5;counter2=1;
for i=1:length(ItiaList)
regionName=ItiaList{i};
idx_temp=LinReg.(regionName).KmeansIdx_select;
GoodBet_temp=LinReg.(regionName).GoodBetas2;
counter=counter2;
for j=1:length(GoodBet_temp)
idx_temp2=find(idx_temp==GoodBet_temp(j));
subplot(xplot,yplot,counter+(j-1));plot(mean(LinReg.(regionName).ZS_AVG(idx_temp2,:),1));title(num2str(length(idx_temp2)));
end
counter2=counter2+yplot;
end

Clusters=[0 1  2 3 0; 0 1  2 5 0; 0 1  2 3 0;0 5  1 2 4; 2 0  5 3 1;3 1  2 5 0; 0 0  2 1 0; 0 0  1 0 0; 0 0  1 4 2;0 0  3 2 1;1 3  2 0 0];
ZS_clusters=cell(1,5);
for i=1:length(ItiaList)
    regionName=ItiaList{i};
    idx_temp=LinReg.(regionName).KmeansIdx_select;
    GoodBet_temp=LinReg.(regionName).GoodBetas2;
    for j=1:length(ZS_clusters)
        if Clusters(i,j)>0            
            idx_temp2=find(idx_temp==GoodBet_temp(Clusters(i,j)));
            if isempty(ZS_clusters{j})
                ZS_clusters{j}=LinReg.(regionName).ZS_AVG(idx_temp2,:);
            else
                ZS_clusters{j}=[ZS_clusters{j}; LinReg.(regionName).ZS_AVG(idx_temp2,:)];
            end
        end
    end
end

x = linspace(0.25,123/4,123);
for i=1:length(ZS_clusters)
    Fighandle=figure;
    set(Fighandle, 'Position', [100, 100, 250, 250]);
    temp=ZS_clusters{i};    
    for k=0:2           
            hold on;
            meanToPlot=mean(temp(:,4+(k*41):40+(k*41)),1);
            plot(x(4+(k*40):40+(k*40)),meanToPlot-mean(meanToPlot(1:4)),'color','k','LineWidth',3);axis([0 30 -3 3]);                       
    end
	hold off;
    print(Fighandle,strcat('_Cluster-',num2str(i)),'-dsvg','-r0');
end

ForCSVExport=cell(1,5);
for j=1:length(ZS_clusters)    
    for region_nb=1:length(ItiaList)
        regionName=ItiaList{region_nb};
        if Clusters(region_nb,j)>0
            idx_temp=LinReg.(regionName).KmeansIdx_select;
            GoodBet_temp=LinReg.(regionName).GoodBetas2;
            idx_temp2=find(idx_temp==GoodBet_temp(Clusters(region_nb,j)));
            ROI_temp=ROIsPerBrain.(regionName).ROIs;
            ROI_temp=ROI_temp(find(LinReg.(regionName).rsquared>0.1),:);
            if isempty(ForCSVExport{j})
                CSV_temp=zeros(length(idx_temp2),4);
                CSV_temp(:,1:3)=ROI_temp(idx_temp2,:);
                CSV_temp(:,4)=Clusters(region_nb,j);
                ForCSVExport{j}=CSV_temp;
            else
                CSV_temp=zeros(length(idx_temp2),4);
                CSV_temp(:,1:3)=ROI_temp(idx_temp2,:);
                CSV_temp(:,4)=Clusters(region_nb,j);
                ForCSVExport{j}=[ForCSVExport{j};CSV_temp];
            end
            
        end
    end
end

for j=1:length(ZS_clusters)
    CSV_temp=ForCSVExport{j};
    CSV_temp(:,3)=CSV_temp(:,3)*2;
    CSV_temp(:,4)=1;
    filename=strcat('ROIs_coord_clust_',num2str(j),'.csv');
    csvwrite(filename,CSV_temp);
end

%Pool everything

region_nb=1;regionName=ItiaList{region_nb};
ZS_pool=LinReg.(regionName).ZS_rsq;
for region_nb=2:length(ItiaList)
    regionName=ItiaList{region_nb};
    ZS_pool=[ZS_pool;LinReg.(regionName).ZS_rsq];    
end
options = statset('UseParallel',1);
[idxKmeans Cmap]=kmeans(ZS_pool,20,'Options',options,'Distance','cityblock','Replicates',5,'MaxIter',1000,'Display','final');
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
    subplot(xplot,yplot,j);plot(mean(ZS_pool(idx_temp,:),1));ylim([-3 3]);
end

% corr_temp=[];
GoodBetas_select=GoodBetas([2 4 5 6 7 14 15 16 19 20]);
% for j=1:length(GoodBetas_select)
%     idx_temp=find(KmeansIdx_select==GoodBetas_select(j));
%     Mean_temp=mean(ZS_pool(idx_temp,:),1);
%     for i=1:length(GoodBetas_select)
%         idx_temp2=find(KmeansIdx_select==GoodBetas_select(i));
%         %temp_corr=corrcoef(Mean_temp, mean(ZS_pool(idx_temp2,:),1));
%         %corr_temp(i,j)=temp_corr(1,2);
%         corr_temp(i,j)=pdist2(Mean_temp, mean(ZS_pool(idx_temp2,:),1),'cityblock');
%     end
% end

figure;
counter=1;counter2=1;xplot=floor(sqrt(length(GoodBetas_select)));yplot=ceil(length(GoodBetas_select)/xplot);coloring={'r','g','y'};
for j=1:length(GoodBetas_select)  
    idx_temp=find(KmeansIdx_select==GoodBetas_select(j));
    subplot(xplot,yplot,j);plot(mean(ZS_pool(idx_temp,:),1));ylim([-3 3]);xlim([1 420]);
end

KmeansIdx_merge=KmeansIdx_select;
idx_temp=ismember(KmeansIdx_select,GoodBetas_select([1 3 10]));
KmeansIdx_merge(idx_temp)=GoodBetas_select(1);
idx_temp=ismember(KmeansIdx_select,GoodBetas_select([2 8]));
KmeansIdx_merge(idx_temp)=GoodBetas_select(2);
idx_temp=ismember(KmeansIdx_select,GoodBetas_select([4 7 9]));
KmeansIdx_merge(idx_temp)=GoodBetas_select(4);

GoodBetas_merge=GoodBetas_select([1 2 4 5 6]);

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

colors=[0.14 1 0.14; 0.7 0.4 1];
for i=1:length(GoodBetas_merge)
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
    print(Fighandle,strcat('_Cluster-',num2str(i)),'-dsvg','-r0');
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


	 
print(Fighandle,strcat('_Cluster-','4_5'),'-dsvg','-r0');


%----------- AVG

region_nb=1;regionName=ItiaList{region_nb};
ZSAVG_pool=LinReg.(regionName).ZS_AVG;
for region_nb=2:length(ItiaList)
    regionName=ItiaList{region_nb};
    ZSAVG_pool=[ZSAVG_pool;LinReg.(regionName).ZS_AVG];    
end
options = statset('UseParallel',1);
[idxKmeans_AVG Cmap_AVG]=kmeans(ZSAVG_pool,20,'Options',options,'Distance','cityblock','Replicates',5,'MaxIter',1000,'Display','final');

options = statset('UseParallel',1);
[idxKmeans_AVG Cmap_AVG]=kmeans(ZSAVG_pool,20,'Options',options,'Distance','cityblock','Replicates',5,'MaxIter',1000,'Display','final');

[~,GoodBetas_AVG]=Test_Regress(Cmap_AVG,Stimuli_AVG,idxKmeans_AVG,0.3);
%Correlate back to Kmeans to remove crap
Threshold=0.4;
KmeansIdx_select_AVG=idxKmeans;
for j=1:length(GoodBetas_AVG)
    idx_g=find(idxKmeans_AVG==GoodBetas_AVG(j));    
    corr_temp=zeros(size(idx_g));
    Mean_temp=mean(ZSAVG_pool(idx_g,:),1);
    for k=1:length(idx_g)
        temp_corr=corrcoef(Mean_temp, ZSAVG_pool(idx_g(k),:));
        corr_temp(k)=temp_corr(1,2);
    end
    KmeansIdx_select_AVG(idx_g(find(corr_temp<=Threshold)))=0;
end
clearvars i j k temp ans coeff coefficients explained i Hindbrain_Mask idx_temp ish_nb

figure;
counter=1;counter2=1;xplot=floor(sqrt(length(GoodBetas_AVG)));yplot=ceil(length(GoodBetas_AVG)/xplot);coloring={'r','g','y'};
for j=1:length(GoodBetas_AVG)  
    idx_temp=find(KmeansIdx_select_AVG==GoodBetas_AVG(j));
    subplot(xplot,yplot,j);plot(mean(ZSAVG_pool(idx_temp,:),1));ylim([-3 3]);
end

figure;imagesc(ZS_pool(find(idxKmeans_AVG==4),:));

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

for i=1:length(GoodBetas_merge)    
    idx_temp=find(KmeansIdx_merge==GoodBetas_merge(i));
    CSV_temp=ROI_pool(idx_temp,:);
    CSV_temp(:,3)=CSV_temp(:,3)*2;
    CSV_temp(:,4)=1;
    filename=strcat('_ROIs_coord_clust_',num2str(i),'.csv');
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
print(Fighandle,strcat('_Cluster-','weird'),'-dsvg','-r0');

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




% figure;
% for i=1:size(Cmap_ZS_8,1)
%     plot(Cmap_ZS_8(i,:));pause
% end
% 
% 
% ModelBothSides=[];
% parfor i=1:length(ZS)
%     %mdl=stepwiselm(Regressor',DataSet(i,:),'linear','Criterion','adjrsquared','Intercept',false,'Upper','interactions','Verbose',0);
%     mdl=fitlm(Stimuli',ZS(i,:),'interactions');
%     ModelBothSides(i).coef=mdl.Coefficients;
%     %ModelResultsSeg_ZS(i).MSE=mdl.MSE;
%     %ModelResultsSeg_ZS(i).Fitted=mdl.Fitted;
%     ModelBothSides(i).rsquared=mdl.Rsquared.Adjusted;
% end
% idx_rsq=find([ModelBothSides.rsquared]>0.15);
% ZS_rsq=ZS(idx_rsq,:);
% 
% Fighandle=figure;
% set(Fighandle, 'Position', [100, 100, 1200, 900]);
% imagesc(x,y,ZS_rsq(randperm(size(ZS_rsq,1)),:),[-0.5 5]);colormap hot;set(gca,'YTickLabel',[]);
% 
% coefficients={};
% for idx=1:length(ModelBothSides)
%     coef=[ModelBothSides(idx).coef];
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
% 
% options = statset('UseParallel',1); [idxKmeans_ZS_rsq Cmap_ZS_rsq]=kmeans(ZS_rsq,30,'Options',options,'Distance','cityblock','Replicates',3,'MaxIter',1000,'Display','final');
% 
% GoodBetas=[5 7 17 19 23 24 26 28];
% 
% idxKmeans_final=zeros(size(ZS,1),1);
% idxKmeans_final(idx_rsq)=idxKmeans_ZS_rsq;
% 
% counter=1;
% x = linspace(0.5,size(ZS,2)/2,size(ZS,2));
% rows=length(GoodBetas);
% counter=1;
% Fighandle=figure;
% set(Fighandle, 'Position', [100, 100, 1300, 900]);
% for i=GoodBetas
%     idx_temp=find(idxKmeans_final==i);
%     subplot(rows,5,counter);plot(mean(ZS(idx_temp,:),1));
%     subplot(rows,5,counter+1);imagesc(ZS(idx_temp,:),[-0.5 4]);
%     subplot(rows,5,counter+2);histogram(idx_Plane(idx_temp));
%     subplot(rows,5,counter+3);histogram(idx_Position(idx_temp),[0.5:1:4.5]);h = gca;h.XTick=[1 2 3 4];h.XTickLabel={'Center','ELO-IRO','ERO-ILO','Out'};
%     subplot(rows,5,counter+4);histogram(idx_Fish(idx_temp));h = gca;h.XTickLabel={'1','2'};
%     counter=counter+5;
% end
% 
% GoodClustersData=[];
% for i=1:length(GoodBetas)
%     GoodClustersData(i).ZS=ZS(idxKmeans_final==GoodBetas(i),:);
%     GoodClustersData(i).Mean=mean(GoodClustersData(i).ZS,1);
%     GoodClustersData(i).STD=std(GoodClustersData(i).ZS,1,1);
% end
% 
% for i=1:numel(GoodClustersData)
%     corr_temp=zeros(size(GoodClustersData(i).ZS,1),1);
%     parfor j=1:size(GoodClustersData(i).ZS,1)
%         temp=corrcoef(GoodClustersData(i).Mean, GoodClustersData(i).ZS(j,:));
%         corr_temp(j)=temp(1,2);
%     end
%     GoodClustersData(i).CorrCoef=corr_temp;
% end
% 
% 
% GoodClusters_goodmembers=[];Threshold=0.5;
% idxKmeans_ZS_goodmembers=zeros(1,size(GoodCalcium,1));
% for i=1:length(GoodBetas)
% %GoodClusters_goodmembers(i).Spikes=GoodClustersData(i).Spikes(find(GoodClustersData(i).CorrCoef>=0.5),:);
% %GoodClusters_goodmembers(i).ZS=zscore(GoodClustersData(i).DF(find(GoodClustersData(i).CorrCoef>=0.5),:),1,2);
% GoodClusters_goodmembers(i).ZS=GoodClustersData(i).ZS(find(GoodClustersData(i).CorrCoef>=Threshold),:);
% temp=find(idxKmeans_final==GoodBetas(i));
% GoodClusters_goodmembers(i).idx=temp(find(GoodClustersData(i).CorrCoef>=0.5));
% GoodClusters_goodmembers(i).mean=mean(GoodClusters_goodmembers(i).ZS,1);
% GoodClusters_goodmembers(i).STD=std(GoodClusters_goodmembers(i).ZS,1,1);
% idx=find(idxKmeans_final==GoodBetas(i));
% idx=idx(find(GoodClustersData(i).CorrCoef>=Threshold));
% idxKmeans_ZS_goodmembers(idx)=GoodBetas(i);
% %GoodClusters_goodmembers(i).Fish=idx_Fish(idx);
% end
% 
% counter=1;
% x = linspace(0.5,size(ZS,2)/2,size(ZS,2));
% rows=length(GoodBetas);
% counter=1;
% Fighandle=figure;
% set(Fighandle, 'Position', [100, 100, 1300, 900]);
% for i=GoodBetas
%     idx_temp=find(idxKmeans_ZS_goodmembers==i);
%     subplot(rows,5,counter);plot(mean(ZS(idx_temp,:),1));
%     subplot(rows,5,counter+1);imagesc(ZS(idx_temp,:),[-0.5 4]);
%     subplot(rows,5,counter+2);histogram(idx_Plane(idx_temp));
%     subplot(rows,5,counter+3);histogram(idx_Position(idx_temp),[0.5:1:4.5]);h = gca;h.XTick=[1 2 3 4];h.XTickLabel={'Center','ELO-IRO','ERO-ILO','Out'};
%     subplot(rows,5,counter+4);histogram(idx_Fish(idx_temp));h = gca;h.XTickLabel={'1','2'};
%     counter=counter+5;
% end
% 
% All_ROIs=[];
% ROIs_idx=[];
% for i = 1:length(MatFiles)
%     name=strcat(MatFiles(i).name);
%     Rs=load(name, 'ROIs');
%     Rs=Rs.ROIs;
%     F=load(name, 'idx_components');
%     F=F.idx_components+1;
%     Rs=Rs(:,F);
%     All_ROIs{i}=Rs;
%     if i==1
%         ROIs_idx(i)=length(F);
%     else
%         ROIs_idx(i)=ROIs_idx(i-1)+length(F);
%     end
% end
% clearvars GC C S F N name i;
% 
% Numbers=[0 [ROIs_idx]];
% temp=[];
% counter=1;
% for i=GoodBetas
%     temp{counter}=find(idxKmeans_ZS_goodmembers==i);
%     %tempidx=find(idxKmeans==idx);
%     %temp{counter}=GoodClusters_goodmembers(counter).idx;
%     counter=counter+1;    
% end
% 
% colors = distinguishable_colors(length(GoodBetas),[1 1 1; 0 0 0]);
% colors = colors*256;
% for idx=1:length(MatFiles)
%     filename=MatFiles(idx).name;
%     ROIsNb=[];ClusterNb=[];
%     %for k = 1 : length(temp)
%     for k = 1 : length(temp)
%         tempROIsNb=find([temp{k}]<=Numbers(idx+1));
%         if tempROIsNb            
%             ROIsNb=[ROIsNb  temp{k}(tempROIsNb)];
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
% 
% Fighandle=figure;
% set(Fighandle, 'Position', [100, 100, 500, 1400]);x = linspace(0.5,size(ZS,2)/2,size(ZS,2));
% counter=1;counter2=1;xplot=1;yplot=length(GoodBetas);%yplot=ceil(length(GoodBetas)/xplot);
% coloring={'r','g','y'};
% for i=GoodBetas
%     idx_temp=find(idxKmeans_ZS_goodmembers==i);
%     subplot(xplot,yplot,counter);plot(x,mean(ZS(idx_temp,:),1),'color',colors(counter2,:)/256);axis([0 131 -1 4]);
%     start=30;
%     for k=1:3        
%         for j=1:3            
%             rectangle('FaceColor',coloring{k},'Position',[start+(j-1)*30 -1 2 0.25]);
%         end
%         start=10+start;
%     end    
%     counter=counter+1;
%     counter2=counter2+1;
% end
