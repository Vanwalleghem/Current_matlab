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
Calcium=Calcium+Noise;
clearvars Noise
Calcium=zscore(Calcium,1,2);

IndexC=strfind({MatFiles.name}, Fish_list{1});
MatFiles_fish = find(not(cellfun('isempty', IndexC)));
if MatFiles_fish(1)==1
        numbersForROIs=[1 [MatFiles(MatFiles_fish).number]];
    else
        numbersForROIs=[MatFiles(MatFiles_fish(1)-1).number [MatFiles(MatFiles_fish).number]];
end    
Calcium_temp=Calcium(numbersForROIs(1):numbersForROIs(end)-1,:);
Calcium_temp=Calcium_temp(:,21:end);
Calcium_temp(:,741:760)=Calcium_temp(:,end-19:end);
Calcium(numbersForROIs(1):numbersForROIs(end)-1,:)=Calcium_temp;

Zbrain_Masks={};
MaskFiles=dir('D:\Pictures\processed\Itia\ROI_ANTs\Masks\*.tif');
progressbar;
for i=1:length(MaskFiles)
    fname =strcat('D:\Pictures\processed\Itia\ROI_ANTs\Masks\',MaskFiles(i).name);
    temp=strsplit(fname,'Masks\');temp=temp{2};
    Brain_region=regexpi(temp,'(.+) -','tokens');
    if isempty(Brain_region)
        Brain_region=regexpi(temp,'(.+).tif','tokens');
    end
    while iscell(Brain_region)
        Brain_region=Brain_region{1};
    end
    SubRegion=regexpi(temp,' - (.+)\.tif','tokens');
    if isempty(SubRegion)
        SubRegion=' ';
    else
        while iscell(SubRegion)
            SubRegion=SubRegion{1};
        end
    end
    Zbrain_Masks{i,1}=Brain_region;
    Zbrain_Masks{i,2}=SubRegion;
    info = imfinfo(fname);
    num_images = numel(info);
    PixelList=[];
    for k = 1:num_images
        %transpose and rotated to align with Itia's stuff
        image_temp = imread(fname, k, 'Info', info)';image_temp=rot90(image_temp);image_temp=image_temp>0;
        [x,y]=find(image_temp>0);
        if x
            if isempty(PixelList)
                PixelList=[x y ones(length(x),1)*k];
            else
                PixelList=vertcat(PixelList,[x y ones(length(x),1)*k]);
            end
        end
    end
    Zbrain_Masks{i,3}=PixelList;
    progressbar(i/length(MaskFiles));
end
clearvars x y PixelList info num_images image_temp Brain_region SubRegion fname MaskFiles i k ans temp
for i=1:length(Zbrain_Masks)
    temp=strsplit(Zbrain_Masks{i,1},' ');temp=temp{1};
    Zbrain_Masks{i,1}=temp;
end

save('Zbrain_Masks.mat','Zbrain_Masks');

%Example for Cereb
IsInMask={};
for i=1:13
    IndexC=strfind({Zbrain_Masks{:,2}}, 'Cerebellum');
    IndexC=find(not(cellfun('isempty', IndexC)));
    Mask=Zbrain_Masks{131,3};%Mask(:,[1 2])=Mask(:,[2 1]);
    ROI_fish=ROIs(i).coord;
    ROI_fish(:,1:2)=round(ROI_fish(:,1:2));
    ROI_fish(:,3)=round(((ROI_fish(:,3)-1)*2)+24);
    IsInMask_temp=ismember(ROI_fish,Mask,'rows');IsInMask{i}=find(IsInMask_temp==1)';    
end
%---------------


temp=cell2mat(IsInMask);

figure;imagesc(Calcium(temp,:),[-1 4]);

load('_aMultipower.mat','idxKmeans_final_goodmember');
load('_aMultipower.mat','MatFiles');
load('_aMultipower.mat','GoodBetas');

Clusters_Brain=cell(length(Fish_list),length(Zbrain_Masks),length(GoodBetas));

for fish_nb=1:13
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
    ROI_fish=find(not(cellfun('isempty', IndexC)));
    ROI_fish=ROIs(ROI_fish).coord;ROI_fish(:,1:2)=round(ROI_fish(:,1:2));
    ROI_fish(:,3)=round(((ROI_fish(:,3)-1)*2)+24);%Get ROIs of fish
    if MatFiles_fish(1)==1
        numbersForROIs=[1 [MatFiles(MatFiles_fish).GoodNumber]];
    else
        numbersForROIs=[MatFiles(MatFiles_fish(1)-1).GoodNumber+1 [MatFiles(MatFiles_fish).GoodNumber]];
    end
    GoodROIs=idxKmeans_final_goodmember(numbersForROIs(1):numbersForROIs(end)-1);
    if find(GoodROIs>0)        
        for brain_nb=1:length(Zbrain_Masks)            
            Mask=Zbrain_Masks{brain_nb,3};%Mask(:,[1 2])=Mask(:,[2 1]);
            counter=1;
            for i=GoodBetas
                idx_temp=find(GoodROIs==i);
                if idx_temp
                    Coords=ROI_fish(idx_temp,:);
                    IsInMask_temp=ismember(Coords,Mask,'rows');Clusters_Brain{fish_nb,brain_nb,counter}=find(IsInMask_temp==1)';                        
                end
                counter=counter+1;
            end
        end
    end
end
clearvars GoodROIs Mask counter idx_temp Coords IsInMask IsInMask_temp


Mean_Clusters_brain=zeros(length(GoodBetas),length(Zbrain_Masks));
STD_Clusters_brain=zeros(length(GoodBetas),length(Zbrain_Masks));
for k=1:size(Clusters_Brain,3)
    for j=1:size(Clusters_Brain,2)
        temp=zeros(1,length(Fish_list));
        for i=1:size(Clusters_Brain,1)
            temp(i)=length(Clusters_Brain{i,j,k});
        end
        Mean_Clusters_brain(k,j)=mean(temp);
        STD_Clusters_brain(k,j)=std(temp);
    end
end

columnsWithAllZeros = all(Mean_Clusters_brain == 0);
columnsWithNoZeros=find(columnsWithAllZeros==0);
Brain_names={};
for i=1:length(MaskFiles)
    fname =strcat('D:\Pictures\processed\Itia\ROI_ANTs\Masks\',MaskFiles(i).name);
    temp=strsplit(fname,'Masks\');temp=temp{2};
    Brain_region=regexpi(temp,'(.+).tif','tokens');
    while iscell(Brain_region)
        Brain_region=Brain_region{1};
    end
    Brain_names{i}=Brain_region;
end

Mean_Clusters_brain_nozeros=Mean_Clusters_brain(:,columnsWithNoZeros);
STD_Clusters_brain_nozeros=STD_Clusters_brain(:,columnsWithNoZeros);
Brain_names_nozeros=Brain_names(columnsWithNoZeros);

temp=max(cumsum(Mean_Clusters_brain_nozeros,2),[],2);
for i=1:length(GoodBetas)
    Mean_Clusters_brain_nozeros_norm(i,:)=Mean_Clusters_brain_nozeros(i,:)/temp(i);
    STD_Clusters_brain_nozeros_norm(i,:)=STD_Clusters_brain_nozeros(i,:)/temp(i);
end
Mean_Clusters_brain_nozeros_norm=Mean_Clusters_brain_nozeros_norm*100;
STD_Clusters_brain_nozeros_norm=STD_Clusters_brain_nozeros_norm*100;

Sorted_Brain=zeros(size(Mean_Clusters_brain_nozeros_norm));
Idx_Brain=zeros(size(Mean_Clusters_brain_nozeros_norm));
for i=1:length(GoodBetas)
    [Sorted_Brain_perClust(i,:) Idx_Brain_perClust(i,:)]=sort(Mean_Clusters_brain_nozeros_norm(i,:));
end

[Sorted_Brain Idx_Brain]=sortrows(Mean_Clusters_brain_nozeros_norm');Sorted_Brain=Sorted_Brain';
Sorted_STD=STD_Clusters_brain_nozeros_norm(:,Idx_Brain);
Sorted_names=Brain_names_nozeros(Idx_Brain);


%Null Distribution
Null_Brain=cell(length(Fish_list),length(Zbrain_Masks));

for fish_nb=1:13
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
    ROI_fish=find(not(cellfun('isempty', IndexC)));
    ROI_fish=ROIs(ROI_fish).coord;ROI_fish(:,1:2)=round(ROI_fish(:,1:2));
    ROI_fish(:,3)=round(((ROI_fish(:,3)-1)*2)+24);%Get ROIs of fish
    if MatFiles_fish(1)==1
        numbersForROIs=[1 [MatFiles(MatFiles_fish).GoodNumber]];
    else
        numbersForROIs=[MatFiles(MatFiles_fish(1)-1).GoodNumber+1 [MatFiles(MatFiles_fish).GoodNumber]];
    end
    GoodROIs=idxKmeans_final_goodmember(numbersForROIs(1):numbersForROIs(end)-1);    
    if find(GoodROIs>0)        
        for brain_nb=1:length(Zbrain_Masks)            
            Mask=Zbrain_Masks{brain_nb,3};%Mask(:,[1 2])=Mask(:,[2 1]);
            counter=1;
            for i=GoodBetas
                idx_temp=find(GoodROIs==i);
                if idx_temp
                    idx_null=randi(length(GoodROIs),size(idx_temp));
                    Coords=ROI_fish(idx_null,:);
                    IsInMask_temp=ismember(Coords,Mask,'rows');Null_Brain{fish_nb,brain_nb,counter}=find(IsInMask_temp==1)';                        
                end
                counter=counter+1;
            end
        end
    end
end
clearvars GoodROIs Mask counter idx_temp Coords IsInMask IsInMask_temp idx_null

%Just the BigRegions
BigBrainRegions=[76 113 259 274 294];
Clusters_BigBrain=cell(length(Fish_list),length(BigBrainRegions),length(GoodBetas));

for fish_nb=1:13
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
    ROI_fish=find(not(cellfun('isempty', IndexC)));
    ROI_fish=ROIs(ROI_fish).coord;ROI_fish(:,1:2)=round(ROI_fish(:,1:2));
    ROI_fish(:,3)=round(((ROI_fish(:,3)-1)*2)+24);%Get ROIs of fish
    if MatFiles_fish(1)==1
        numbersForROIs=[1 [MatFiles(MatFiles_fish).GoodNumber]];
    else
        numbersForROIs=[MatFiles(MatFiles_fish(1)-1).GoodNumber+1 [MatFiles(MatFiles_fish).GoodNumber]];
    end
    GoodROIs=idxKmeans_final_goodmember(numbersForROIs(1):numbersForROIs(end)-1);
    if find(GoodROIs>0)        
        for brain_nb=1:length(BigBrainRegions)            
            Mask=Zbrain_Masks{BigBrainRegions(brain_nb),3};%Mask(:,[1 2])=Mask(:,[2 1]);
            counter=1;
            for i=GoodBetas
                idx_temp=find(GoodROIs==i);
                if idx_temp
                    Coords=ROI_fish(idx_temp,:);
                    IsInMask_temp=ismember(Coords,Mask,'rows');Clusters_BigBrain{fish_nb,brain_nb,counter}=find(IsInMask_temp==1)';                        
                end
                counter=counter+1;
            end
        end
    end
end
clearvars GoodROIs Mask counter idx_temp Coords IsInMask IsInMask_temp

Clusters_BigBrain=Clusters_Brain(:,BigBrainRegions,:);

Nb_Clusters_Bigbrain=zeros(size(Clusters_BigBrain));
Mean_Clusters_Bigbrain=zeros(length(GoodBetas),length(BigBrainRegions));
STD_Clusters_Bigbrain=zeros(length(GoodBetas),length(BigBrainRegions));
for k=1:size(Clusters_BigBrain,3)
    for j=1:size(Clusters_BigBrain,2)
        temp=zeros(1,length(Fish_list));
        for i=1:size(Clusters_BigBrain,1)
            temp(i)=length(Clusters_BigBrain{i,j,k});
            Nb_Clusters_Bigbrain(i,j,k)=length(Clusters_BigBrain{i,j,k});
        end
        Mean_Clusters_Bigbrain(k,j)=mean(temp);
        STD_Clusters_Bigbrain(k,j)=std(temp);
    end
end

temp=Nb_Clusters_Bigbrain;temp(temp==0)=nan;
temp2=temp(:,:,1);

% Results per class of brain region
Classes_brain=unique(Zbrain_Masks(:,1));
for i=1:length(Classes_brain);
    temp=Classes_brain(i);temp=temp{1};
    IndexC=strfind(Zbrain_Masks(:,1),temp);
    Idx_BigBrain = find(not(cellfun('isempty', IndexC)));Idx_BigBrain(ismember(Idx_BigBrain,BigBrainRegions))=[];    
    ClustersPerClass.(temp)=Mean_Clusters_brain(:,Idx_BigBrain);
    temp2=strcat(temp,'_sorted');
    [Sorted_Brain Idx_Brain]=sortrows(Mean_Clusters_brain(:,Idx_BigBrain)');Sorted_Brain=Sorted_Brain';
    columnsWithAllZeros = all(Sorted_Brain == 0);
    columnsWithNoZeros=find(columnsWithAllZeros==0);
    ClustersPerClass.(temp2)=Sorted_Brain(:,columnsWithNoZeros);
    temp2=strcat(temp,'_sortedIdx');
    ClustersPerClass.(temp2)=Idx_BigBrain(Idx_Brain(columnsWithNoZeros));
end

idx_temp=ClustersPerClass.Diencephalon_sortedIdx;
PrismReady=ones(length(GoodBetas)*3,length(idx_temp))*13;
for i=1:3:length(GoodBetas)*3
    PrismReady(i,:)=ClustersPerClass.Diencephalon_sorted(ceil(i/3),:);
    PrismReady(i+1,:)=STD_Clusters_brain(ceil(i/3),idx_temp);
end
Brain_prismNames=Brain_names(idx_temp);

idx_temp=ClustersPerClass.Mesencephalon_sortedIdx;
PrismReady=ones(length(GoodBetas)*3,length(idx_temp))*13;
for i=1:3:length(GoodBetas)*3
    PrismReady(i,:)=ClustersPerClass.Mesencephalon_sorted(ceil(i/3),:);
    PrismReady(i+1,:)=STD_Clusters_brain(ceil(i/3),idx_temp);
end
Brain_prismNames=Brain_names(idx_temp);

idx_temp=ClustersPerClass.Telencephalon_sortedIdx;
PrismReady=ones(length(GoodBetas)*3,length(idx_temp))*13;
for i=1:3:length(GoodBetas)*3
    PrismReady(i,:)=ClustersPerClass.Telencephalon_sorted(ceil(i/3),:);
    PrismReady(i+1,:)=STD_Clusters_brain(ceil(i/3),idx_temp);
end
Brain_prismNames=Brain_names(idx_temp);

idx_temp=ClustersPerClass.Rhombencephalon_sortedIdx;
PrismReady=ones(length(GoodBetas)*3,length(idx_temp))*13;
for i=1:3:length(GoodBetas)*3
    PrismReady(i,:)=ClustersPerClass.Rhombencephalon_sorted(ceil(i/3),:);
    PrismReady(i+1,:)=STD_Clusters_brain(ceil(i/3),idx_temp);
end
Brain_prismNames=Brain_names(idx_temp);



