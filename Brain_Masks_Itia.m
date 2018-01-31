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
