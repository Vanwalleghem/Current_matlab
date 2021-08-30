CSV_Files=dir('_ROIs*.csv');
ROIs2=struct();
for i=1:length(CSV_Files);
    temp=csvread(CSV_Files(i).name,1);
    Fishname=regexp(CSV_Files(i).name,'_ROIsFish_(\d+).csv','tokens');Fishname=Fishname{1}{1};
    ROIs2(i).name=Fishname;    
    ROIs2(i).coord=temp(:,1:3);
    ROIs2(i).idx=temp(:,5);
end
clearvars i temp CSV_Files Fishname

CSV_names={CSV_Files.name};
ROIs2=struct();
for fish_nb=1:length(Fish_list)
    IndexC=strfind(CSV_names,num2str(Fish_list(fish_nb)));
    i=find(not(cellfun('isempty', IndexC)));
    temp=csvread(CSV_Files(i).name,1);
    Fishname=regexp(CSV_Files(i).name,'_ROIsFish_(\d+).csv','tokens');Fishname=Fishname{1}{1};
    ROIs2(fish_nb).name=Fishname;    
    ROIs2(fish_nb).coord=temp(:,1:3);
    ROIs2(fish_nb).idx=temp(:,5);
end

i=1;ROI_pool2=ROIs2(i).coord;
for i=2:length(ROIs2);
    ROI_pool2=[ROI_pool2; ROIs2(i).coord];
end

Sort_ROIs2=[];temp_nb=0;
MatFiles_names={MatFiles.name};
for fish_nb=1:length(Fish_list)
    IndexC=find(idx_Fish==Fish_list(fish_nb));
    IndexC=ismember(ROIs_idx,IndexC);
    MatFiles_fish = find(IndexC>0);
    for file_nb=1:length(MatFiles_fish)
        if MatFiles_fish(file_nb)==1
            numbersForROIs=[1 MatFiles(1).GoodNumber];
        else
            numbersForROIs=[MatFiles(MatFiles_fish(file_nb)-1).GoodNumber+1 MatFiles(MatFiles_fish(file_nb)).GoodNumber];
        end
        if ismember(numbersForROIs,Sort_ROIs2)
            fish_name
            break
        end
        Sort_ROIs2=[Sort_ROIs2 numbersForROIs(1):1:numbersForROIs(2)];
    end
    if ~length(Sort_ROIs2)-temp_nb==sum(idx_Fish==(Fish_list(fish_nb)))
        fish_name
        break
    end
    temp_nb=length(Sort_ROIs2);
end
clearvars slice fish_nb roi_nb ROI Centroids IndexC file_nb

ROI_fish2(Sort_ROIs2,:)=ROI_pool2;

Fish_list=unique(idx_Fish);
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
        [slice,~]=regexp(filename,'_(\d+)um_','tokens','match');slice=str2num(slice{1}{1});
        idx_name=strcat(num2str(Fish_list(fish_nb)),'_',num2str(slice));
         if plane==1
            temp_roi=0;
        else
            temp_roi=temp_roi+size(ROI,3);
        end
        ROI=All_ROIs{MatFiles_fish(plane)};
        imagename=regexp(filename,'_output_analysis','split');
        imagename=strcat(imagename{1},'_mean.tiff');
        image=double(imread(imagename));image=image/max(max(image));image=image*128;
        ROI=reshape(full(ROI),size(image,1),size(image,2),size(ROI,2));                
        for roi_nb=1:size(ROI,3)
            progressbar([],[],roi_nb/size(ROI,3));
            temp=regionprops(uint16(squeeze(ROI(:,:,roi_nb)))==max(max(uint16(squeeze(ROI(:,:,roi_nb))))),'Centroid');
            Centroids(roi_nb+temp_roi,5)=roi_nb+temp_roi;
            temp=temp.Centroid;
            Centroids(roi_nb+temp_roi,1:2)=temp;
            if filename(2)=='9'
                Centroids(roi_nb+temp_roi,3)=((slice-12376)/20);
            else
                Centroids(roi_nb+temp_roi,3)=(slice/20);
            end
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

Fish_list_backup=Fish_list;

CSV_Files=dir('_4Warped*.csv');
CSV_names={CSV_Files.name};
ROIs=struct();
Fish_list(find(Fish_list>620000 & Fish_list<720000))=[];
for fish_nb=1:length(Fish_list)
    temp=num2str(Fish_list(fish_nb));    
    IndexC=strfind(CSV_names,temp);
    i=find(not(cellfun('isempty', IndexC)));
    Fishname=temp;
    temp=csvread(CSV_Files(i).name,1);
%     Fishname=regexp(CSV_Files(fish_nb).name,'_4Warped(\d+)1234(\d).+.csv','tokens');Fishname=Fishname{1}{1};
%     Fishname=strrep(Fishname,'_','1234');
    ROIs(fish_nb).name=Fishname;    
    ROIs(fish_nb).coord=temp(:,1:3);
    ROIs(fish_nb).idx=temp(:,5);
end
clearvars i temp CSV_Files Fishname

i=1;ROI_pool=ROIs(i).coord;
for i=2:length(ROIs);
    ROI_pool=[ROI_pool; ROIs(i).coord];
end

temp_nb=0;Sorted_ROIs=nan(length(ZS2),3);
MatFiles_names={MatFiles.name};
for fish_nb=1:length(Fish_list)
    temp=num2str(ROIs(fish_nb).name);
    temp=strsplit(temp,'1234');
    if temp{2}=='0'
        IndexC=regexp({MatFiles.name},strcat('[fFh]',temp{1},'\D.+\D_o'),'tokens');    
    else
        IndexC=regexp({MatFiles.name},strcat('[fFh]',temp{1},'\D.+_',temp{2},'_'),'tokens');    
    end
    MatFiles_fish = find(not(cellfun('isempty', IndexC)));
    ROI=ROIs(fish_nb).coord;
    Sort_ROIs=[];Counter_ROI_coord=1;
    for file_nb=1:length(MatFiles_fish)
        if file_nb>1
            Counter_ROI_coord=Counter_ROI_coord+length(Sorted_ROIs(numbersForROIs(1):numbersForROIs(2)));
        end
        if MatFiles_fish(file_nb)==1
            numbersForROIs=[1 MatFiles(1).GoodNumber];
        else
            numbersForROIs=[MatFiles(MatFiles_fish(file_nb)-1).GoodNumber+1 MatFiles(MatFiles_fish(file_nb)).GoodNumber];
        end
        if ismember(numbersForROIs,Sort_ROIs)
            temp
            break
        end        
        Sorted_ROIs(numbersForROIs(1):numbersForROIs(2),:)=ROI(Counter_ROI_coord:Counter_ROI_coord+length(Sorted_ROIs(numbersForROIs(1):numbersForROIs(2)))-1,:);
    end    
    if ~length(Sort_ROIs)-temp_nb==sum(idx_Fish==(Fish_list(fish_nb)))
        temp
        break
    end
    temp_nb=length(Sort_ROIs);
end
clearvars slice fish_nb roi_nb ROI Centroids IndexC file_nb
sum(~isnan(Sorted_ROIs(:,1))) %How many errors there are

IndexC=regexp({MatFiles.name},strcat('[fFh]',temp{1},'\D.+_+o'),'tokens');
MatFiles_fish = find(not(cellfun('isempty', IndexC)));
temp_nb=0;Sorted_ROIs=nan(length(ZS2),3);
MatFiles_names={MatFiles.name};
for fish_nb=1:length(Fish_list)
    temp=num2str(ROIs(fish_nb).name);
    temp=strsplit(temp,'1234');
    if temp{2}=='0'
        IndexC=regexp({MatFiles.name},strcat('[fFh]',temp{1},'\D.+\D_o'),'tokens');
    else
        IndexC=regexp({MatFiles.name},strcat('[fFh]',temp{1},'\D.+_',temp{2},'_'),'tokens');
    end
    MatFiles_fish = find(not(cellfun('isempty', IndexC)));
    ROI=ROIs(fish_nb).coord;
    Sort_ROIs=[];Counter_ROI_coord=1;
    for file_nb=1:length(MatFiles_fish)
        if file_nb>1
            Counter_ROI_coord=Counter_ROI_coord+length(Sorted_ROIs(numbersForROIs(1):numbersForROIs(2)));
        end
        if MatFiles_fish(file_nb)==1
            numbersForROIs=[1 MatFiles(1).GoodNumber];
        else
            numbersForROIs=[MatFiles(MatFiles_fish(file_nb)-1).GoodNumber+1 MatFiles(MatFiles_fish(file_nb)).GoodNumber];
        end
        if ismember(numbersForROIs,Sort_ROIs)
            temp
            break
        end
        Sorted_ROIs(numbersForROIs(1):numbersForROIs(2),:)=ROI(Counter_ROI_coord:Counter_ROI_coord+length(Sorted_ROIs(numbersForROIs(1):numbersForROIs(2)))-1,:);
    end
    if ~length(Sort_ROIs)-temp_nb==sum(idx_Fish==(Fish_list(fish_nb)))
        temp
        break
    end
    temp_nb=length(Sort_ROIs);
end
clearvars slice fish_nb roi_nb ROI Centroids IndexC file_nb
sum(~isnan(Sorted_ROIs(:,1))) %How many errors there are
unique(idx_Fish(isnan(Sorted_ROIs(:,1))))

ROI_temp=round(Sorted_ROIs);
dims_Zbrain=[1406 621];
ROI_rotated=ROI_temp;
ROI_rotated(:,2)=abs(ROI_temp(:,2)-dims_Zbrain(2));
ROI_rotated(:,1)=abs(ROI_temp(:,1)-dims_Zbrain(1));
ROI_rotated(:,3)=round(ROI_temp(:,3)/2);
clearvars ROI_temp

load('D:\Pictures\processed\Zbrain_Masks.mat')
Zbrain_AllMask=vertcat(Zbrain_Masks{[1:1:77 79:1:294],3});
Zbrain_AllMask=unique(Zbrain_AllMask,'rows');
idx_rand=randsample(length(Zbrain_AllMask),10000); %need to select a few members only
Fighandle=figure;
set(Fighandle, 'Position', [10,10, 600, 1400]);
scatter(Zbrain_AllMask(idx_rand,2),Zbrain_AllMask(idx_rand,1),'.')
hold on;
scatter(ROI_rotated(:,2),ROI_rotated(:,1),2,'r','.'); hold off;


ROI_temp=round(Sorted_ROIs);
dims_Zbrain=[1406 621];
ROI_rotated=ROI_temp;
ROI_rotated(:,2)=abs(ROI_temp(:,2)-dims_Zbrain(2));
ROI_rotated(:,1)=abs(ROI_temp(:,1)-dims_Zbrain(1));
ROI_rotated(:,3)=round(ROI_temp(:,3)/2);
clearvars ROI_temp

idx_NotF9orF7=find(~ismember(idx_Fish,[unique(idx_Fish(isnan(Sorted_ROIs(:,1))))' 912340 912341 912342]));

