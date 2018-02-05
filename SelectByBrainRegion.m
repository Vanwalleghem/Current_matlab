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
clearvars GC C S F N name i GS Fitness
ZS=zscore(Calcium+Noise,1,2);
ZS2=zscore(GoodCalcium+GoodNoise,1,2);
clearvars GoodCalcium GoodNoise Calcium Noise

i=15;
Coords=ROIs(i).coord;
fish_nb=1;

if iscell(Fish_list)
    Fish_name=Fish_list{fish_nb};
else
    Fish_name=num2str(Fish_list(fish_nb));
end
IndexC=strfind({MatFiles.name}, Fish_name);
MatFiles_fish = find(not(cellfun('isempty', IndexC)));
ROI_name=strsplit(Fish_name,'Fish2017');
if iscell(ROI_name)
    ROI_name=ROI_name{1};
end
IndexC=strfind({ROIs.name},ROI_name);
ROI_fish=find(not(cellfun('isempty', IndexC)));ROI_fish=ROI_fish(2);
ROI_fish_good=ROIs_good(ROI_fish).coord;ROI_fish_good(:,1:2)=round(ROI_fish_good(:,1:2));
ROI_fish_good(:,3)=round(((ROI_fish_good(:,3)-1)*2)+24);
ROI_fish=ROIs(ROI_fish).coord;ROI_fish(:,1:2)=round(ROI_fish(:,1:2));
ROI_fish(:,3)=round(((ROI_fish(:,3)-1)*2)+24);

if MatFiles_fish(1)==1
    numbersForROIs=[1 [MatFiles(MatFiles_fish).GoodNumber]];
else
    numbersForROIs=[MatFiles(MatFiles_fish(1)-1).GoodNumber+1 [MatFiles(MatFiles_fish).GoodNumber]];
end

%Example for Cereb
IsInCereb=[];
IndexC=strfind({Zbrain_Masks{:,2}}, 'Cerebellum');
IndexC=find(not(cellfun('isempty', IndexC)));
for i=IndexC
    if isempty(Mask)
        Mask=Zbrain_Masks{IndexC,3};
    else
        Mask=vertcat(Mask,Zbrain_Masks{IndexC,3});
    end
end
IsInCereb=ismember(ROI_fish,Mask,'rows');IsInCereb_good=ismember(ROI_fish_good,Mask,'rows');
Cereb=ZS(IsInCereb,:);
Cereb_good=ZS2(IsInCereb_good,:);

options = statset('UseParallel',1); [idxKmeans_cereb Cmap_cereb]=kmeans(Cereb,20,'Options',options,'Distance','cityblock','Replicates',5,'MaxIter',1000,'Display','final');
%---------------
figure;imagesc(Cmap_cereb,[-1 5]);colormap hot;

%Example for Thalamus
IsInMask={};

    IndexC=strfind({Zbrain_Masks{:,2}}, 'Thalamus');
    IndexC=find(not(cellfun('isempty', IndexC)));
    Mask=[];
    for i=IndexC
        if isempty(Mask)
            Mask=Zbrain_Masks{IndexC,3};
        else
            Mask=vertcat(Mask,Zbrain_Masks{IndexC,3});
        end        
    end    
    IsInThalamus=ismember(ROI_fish,Mask,'rows');IsInThalamus_good=ismember(ROI_fish_good,Mask,'rows');
    Thalamus=ZS(IsInThalamus,:);
Thalamus_good=ZS2(IsInThalamus_good,:);
    


Thalamus=ZS(IsInMask_temp,:);
Cereb_good=ZS2(IsInCereb_good,:);

options = statset('UseParallel',1); [idxKmeans_Thalamus Cmap_Thalamus]=kmeans(Thalamus,20,'Options',options,'Distance','cityblock','Replicates',5,'MaxIter',1000,'Display','final');
%---------------
figure;imagesc(Cmap_Thalamus,[-1 5]);colormap hot;

PerBrainRegions=struct();
ItiaList={'Thalamus','Cerebellum','NucMLF','Semicircularis','Pallium','Subpallium'};
progressbar;
for i=5:length(ItiaList)
    progressbar(i/length(ItiaList),[]);
    for fish_nb=1:13
        progressbar([],fish_nb/13);
        regionName=ItiaList{i};
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
        ROI_fish_good=ROIs_good(ROI_fish).coord;ROI_fish_good(:,1:2)=round(ROI_fish_good(:,1:2));
        ROI_fish_good(:,3)=round(((ROI_fish_good(:,3)-1)*2)+24);
        ROI_fish=ROIs(ROI_fish).coord;ROI_fish(:,1:2)=round(ROI_fish(:,1:2));
        ROI_fish(:,3)=round(((ROI_fish(:,3)-1)*2)+24);
        if MatFiles_fish(1)==1
            numbersForROIs_good=[1 [MatFiles(MatFiles_fish).GoodNumber]];
            numbersForROIs=[1 [MatFiles(MatFiles_fish).number]];
        else
            numbersForROIs_good=[MatFiles(MatFiles_fish(1)-1).GoodNumber+1 [MatFiles(MatFiles_fish).GoodNumber]];
            numbersForROIs=[MatFiles(MatFiles_fish(1)-1).number+1 [MatFiles(MatFiles_fish).number]];
        end
        IsInBrainRegion=[];
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
        IsInBrainRegion=ismember(ROI_fish,Mask,'rows');IsInBrainRegion_good=ismember(ROI_fish_good,Mask,'rows');
        PerBrainRegions(fish_nb).(regionName).ROIsCent=ROI_fish(IsInBrainRegion,:);
        PerBrainRegions(fish_nb).(regionName).ROIsCent_good=ROI_fish_good(IsInBrainRegion_good,:);
        temp=ZS(numbersForROIs(1):numbersForROIs(end),:);IsInBrainRegion=temp(IsInBrainRegion,:);
        temp=ZS2(numbersForROIs_good(1):numbersForROIs_good(end),:);IsInBrainRegion_good=temp(IsInBrainRegion_good,:);
        PerBrainRegions(fish_nb).(regionName).ZS=IsInBrainRegion;
        PerBrainRegions(fish_nb).(regionName).ZS_good=IsInBrainRegion_good;        
    end
end

%Now to make it into one big pool of data
ZS_Brain=struct();
ZS_Brain_good=struct();
for i=1:length(ItiaList)
    regionName=ItiaList{i};        
    for j=1:length(PerBrainRegions)
        if j==1;
            ZS_Brain.(regionName)=PerBrainRegions(j).(regionName).ZS;
            ZS_Brain_good.(regionName)=PerBrainRegions(j).(regionName).ZS_good;
        else
            ZS_Brain.(regionName)=vertcat(ZS_Brain.(regionName),PerBrainRegions(j).(regionName).ZS);
            ZS_Brain_good.(regionName)=vertcat(ZS_Brain_good.(regionName), PerBrainRegions(j).(regionName).ZS_good);        
        end
    end
end

%LinReg of all
LinReg=struct();
LinReg_good=struct();
for j=3:length(ItiaList)
    regionName=ItiaList{j};
    temp=ZS_Brain.(regionName);	
    coefficients=struct();
    rsquared=zeros(1,length(temp));
    parfor i=1:size(temp,1)
        mdl=fitlm(Stimuli',temp(i,:));%,'interactions');
        coefficients(i).coef=mdl.Coefficients;
        rsquared(i)=mdl.Rsquared.Adjusted;
    end    
    LinReg.(regionName).coef=coefficients;
    LinReg.(regionName).rsquared=rsquared;
    temp=ZS_Brain_good.(regionName);	
    coefficients=struct();
    rsquared=zeros(1,length(temp));
    parfor i=1:size(temp,1)
        mdl=fitlm(Stimuli',temp(i,:));%,'interactions');
        coefficients(i).coef=mdl.Coefficients;
        rsquared(i)=mdl.Rsquared.Adjusted;
    end    
    LinReg_good.(regionName).coef=coefficients;
    LinReg_good.(regionName).rsquared=rsquared;    
end

options = statset('UseParallel',1); 
for j=1:length(ItiaList)
    regionName=ItiaList{j};
    idx_temp=find(LinReg.(regionName).rsquared>0.1);
    LinReg.(regionName).ZS_rsq=ZS_Brain.(regionName)(idx_temp,:);
    [idxKmeans Cmap]=kmeans(LinReg.(regionName).ZS_rsq,10,'Options',options,'Distance','cityblock','Replicates',3,'MaxIter',1000,'Display','final');
    LinReg.(regionName).KmeansCenter=Cmap;
    LinReg.(regionName).KmeansIdx=idxKmeans;
    idx_temp=find(LinReg_good.(regionName).rsquared>0.1);
    LinReg_good.(regionName).ZS_rsq=ZS_Brain_good.(regionName)(idx_temp,:);
    [idxKmeans Cmap]=kmeans(LinReg_good.(regionName).ZS_rsq,5,'Options',options,'Distance','cityblock','Replicates',3,'MaxIter',1000,'Display','final');
    LinReg_good.(regionName).KmeansCenter=Cmap;
    LinReg_good.(regionName).KmeansIdx=idxKmeans;
end