MatFiles=dir('*analysis_matlab.mat');

MatFiles=dir('*ELO_IRO*analysis_matlab.mat');

idx_Fish_name={};
for i=1:length(MatFiles)
    name=strcat(MatFiles(i).name);
    %[Fish,~]=regexp(name,'Fish2017(\d+_E\DO)_','tokens');Fish=Fish{1};
    [Fish,~]=regexp(name,'Fish2018(\d+)_fish(\d)','tokens');Fish=str2num(strcat(Fish{1}{1},Fish{1}{2}));
    if iscell(Fish)
        Fish=Fish{1};
    end
    idx_Fish_name{i}=Fish;
end
clearvars i Fish Plane name counter
Fish_list=unique(cell2mat(idx_Fish_name));

%Creates ROI csv files
Errored_ROI={};
progressbar(0,0,0);
for fish_nb=1:length(Fish_list)
    if iscell(Fish_list)
        IndexC=strfind({MatFiles.name}, Fish_list{fish_nb});
    else
        IndexC=strfind({MatFiles.name}, num2str(Fish_list(fish_nb)));
    end       
    MatFiles_fish = find(not(cellfun('isempty', IndexC)));
    progressbar(fish_nb/length(Fish_list));
    Centroids=zeros(1,5);    
    for file_nb=1:length(MatFiles_fish)
        progressbar([],file_nb/length(MatFiles_fish));
        name=MatFiles(MatFiles_fish(file_nb)).name;
        Rs=load(name, 'ROIs');
        Rs=Rs.ROIs;
        %F=load(name, 'idx_components');%Include all the ROIs
        %F=F.idx_components+1;
        %Rs=Rs(:,F);
        cor_name=strrep(name,'analysis_matlab','correlation');
        cor_im=load(cor_name);cor_im=cor_im.Correlation_image;
        dims=size(cor_im);
        if file_nb==1
            temp_roi=0;
        else
            temp_roi=temp_roi+size(ROI,3);
        end
        ROI=reshape(full(Rs),dims(1),dims(2),size(Rs,2));counter=1;
        [slice,~]=regexp(name,'_(\d+)_','tokens','match');slice=str2num(slice{1}{1});
        for roi_nb=1:size(ROI,3)
            progressbar([],[],roi_nb/size(ROI,3));
            temp=regionprops(uint16(squeeze(ROI(:,:,roi_nb)))==max(max(uint16(squeeze(ROI(:,:,roi_nb))))),'Centroid');
            Centroids(roi_nb+temp_roi,5)=roi_nb+temp_roi;
            %Centroids(roi_nb+temp_roi,1:2)=temp.Centroid; With the resize
            %need to multiply coords by 1.5
            temp=temp.Centroid;
            Centroids(roi_nb+temp_roi,1:2)=temp;
            Centroids(roi_nb+temp_roi,3)=slice;
        end     
        
    end    
    if iscell(Fish_list)
        image_name=strcat('_ROIsFish',Fish_list{fish_nb},'.csv');
    else
        image_name=strcat('_ROIsFish',num2str(Fish_list(fish_nb)),'.csv');
    end
    csvwrite(image_name,Centroids);
end
clearvars slice fish_nb roi_nb ROI Centroids IndexC file_nb

%-------------------------------

idx_Fish=zeros(1,length(MatFiles));
for i=1:length(MatFiles)    
    name=strcat(MatFiles(i).name);
    [Fish,~]=regexp(name,'Fish2017(\d+)_','tokens','match');Fish=str2num(Fish{1}{1});
    idx_Fish(i)=Fish;
end
clearvars i Fish Plane name counter

Fish_list=unique(idx_Fish);

 idx_Fish={};
for i=1:length(MatFiles)    
    name=strcat(MatFiles(i).name);
    [Fish,~]=regexp(name,'Fish2017(\d+_E\DO)_','tokens');Fish=Fish{1};
    if iscell(Fish)
        Fish=Fish{1};
    end
    idx_Fish{i}=Fish;
end
clearvars i Fish Plane name counter

Fish_list=unique(idx_Fish);

% Errored_ROI={};
% progressbar(0,0,0);
% for fish_nb=1:length(Fish_list)
%     IndexC=strfind({MatFiles.name}, num2str(Fish_list(fish_nb)));
%     MatFiles_fish = find(not(cellfun('isempty', IndexC)));
%     image_name=strcat('_ROIsFish',num2str(Fish_list(fish_nb)),'.tif');
%     progressbar(fish_nb/length(Fish_list));
%     for slice=1:25
%         progressbar([],slice/25);
%         slice_nb=strcat('_',num2str(slice),'_');
%         IndexC=strfind({MatFiles(MatFiles_fish).name}, slice_nb);
%         file_nb=find(not(cellfun('isempty', IndexC)));
%         if isempty(file_nb)
%             image=zeros(size(cor_im));radius =1;
%         else
%             name=strcat(MatFiles(MatFiles_fish(file_nb)).name);
%             Rs=load(name, 'ROIs');
%             Rs=Rs.ROIs;
%             F=load(name, 'idx_components');
%             F=F.idx_components+1;
%             Rs=Rs(:,F);
%             cor_name=strrep(name,'analysis_matlab','correlation');
%             cor_im=load(cor_name);cor_im=cor_im.Correlation_image;
%             dims=size(cor_im);
%             if slice==1
%                 temp_roi=0;
%             else
%                 temp_roi=temp_roi+size(ROI,3);
%             end
%             ROI=reshape(full(Rs),dims(1),dims(2),size(Rs,2));counter=1;
%             Centroids={};
%             for roi_nb=1:size(ROI,3)
%                 temp=regionprops(uint16(squeeze(ROI(:,:,roi_nb)))==max(max(uint16(squeeze(ROI(:,:,roi_nb))))),'Centroid');
%                 Centroids{roi_nb}=temp.Centroid;
%             end
%             %imageSizeX = dims(1);
%             %imageSizeY = dims(2);
%             %[columnsInImage rowsInImage] = meshgrid(1:imageSizeX, 1:imageSizeY);
%             image=zeros(size(cor_im));radius =1;
%             for center=1:length(Centroids)
%                 progressbar([],[],center/length(Centroids));
%                 centerY=round(Centroids{center}(1));centerX=round(Centroids{center}(2));
%                 %circlePixels = (rowsInImage - centerY).^2 + (columnsInImage - centerX).^2 <= radius.^2;
%                 %image(circlePixels')=center;
%                 image(centerX,centerY)=center+temp_roi;
%                 if ~length(unique(image))-1==size(Rs,2)
%                     Errored_ROI{counter}=file_nb;
%                     counter=counter+1;
%                 end
%             end
%             %figure;imagesc(image)
%         end
%         imwrite(uint16(image),image_name,'WriteMode', 'append');
%     end
% end
% 
% idx_Fish={};
% for i=1:length(MatFiles)    
%     name=strcat(MatFiles(i).name);
%     [Fish,~]=regexp(name,'Fish2017(\d+_E\DO)_','tokens');Fish=Fish{1};
%     if iscell(Fish)
%         Fish=Fish{1};
%     end
%     idx_Fish{i}=Fish;
% end
% clearvars i Fish Plane name counter
% 
% Fish_list=unique(idx_Fish);
% 
% Errored_ROI={};
% progressbar(0,0,0);
% for fish_nb=1:length(Fish_list)
%     IndexC=strfind({MatFiles.name}, Fish_list{fish_nb});
%     MatFiles_fish = find(not(cellfun('isempty', IndexC)));
%     image_name=strcat('_ROIsFish',Fish_list{fish_nb},'.tif');
%     progressbar((fish_nb-1)/length(Fish_list));
%     for slice=1:25
%         progressbar([],slice/25);
%         slice_nb=strcat('_',num2str(slice),'_');
%         IndexC=strfind({MatFiles(MatFiles_fish).name}, slice_nb);
%         file_nb=find(not(cellfun('isempty', IndexC)));
%         if isempty(file_nb)
%             image=zeros(size(cor_im));radius =1;
%         else
%             name=strcat(MatFiles(MatFiles_fish(file_nb)).name);
%             Rs=load(name, 'ROIs');
%             Rs=Rs.ROIs;
%             F=load(name, 'idx_components');
%             F=F.idx_components+1;
%             Rs=Rs(:,F);
%             cor_name=strrep(name,'analysis_matlab','correlation');
%             cor_im=load(cor_name);cor_im=cor_im.Correlation_image;
%             dims=size(cor_im);
%             if slice==1
%                 temp_roi=0;
%             else
%                 temp_roi=temp_roi+size(ROI,3);
%             end
%             ROI=reshape(full(Rs),dims(1),dims(2),size(Rs,2));counter=1;
%             Centroids={};
%             for roi_nb=1:size(ROI,3)
%                 temp=regionprops(uint16(squeeze(ROI(:,:,roi_nb)))==max(max(uint16(squeeze(ROI(:,:,roi_nb))))),'Centroid');
%                 Centroids{roi_nb}=temp.Centroid;
%             end
%             %imageSizeX = dims(1);
%             %imageSizeY = dims(2);
%             %[columnsInImage rowsInImage] = meshgrid(1:imageSizeX, 1:imageSizeY);
%             image=zeros(size(cor_im));radius =1;
%             for center=1:length(Centroids)
%                 progressbar([],[],center/length(Centroids));
%                 centerY=round(Centroids{center}(1));centerX=round(Centroids{center}(2));                
%                 image(centerX,centerY)=center+temp_roi;
%                 if ~length(unique(image))-1==size(Rs,2)
%                     Errored_ROI{counter}=file_nb;
%                     counter=counter+1;
%                 end
%             end
%             %figure;imagesc(image)
%         end
%         imwrite(uint16(image),image_name,'WriteMode', 'append');
%     end
% end
% 
% Errored_ROI={};
% progressbar(0,0,0);
% for fish_nb=1:length(Fish_list)
%     if iscell(Fish_list)
%         IndexC=strfind({MatFiles.name}, Fish_list{fish_nb});
%     else
%         IndexC=strfind({MatFiles.name}, num2str(Fish_list(fish_nb)));
%     end   
%     
%     MatFiles_fish = find(not(cellfun('isempty', IndexC)));
%     progressbar(fish_nb/length(Fish_list));
%     Centroids=zeros(1,6,'uint16');
%     
%     for file_nb=1:length(MatFiles_fish)
%         progressbar([],file_nb/length(MatFiles_fish));
%         name=MatFiles(MatFiles_fish(file_nb)).name;
%         Rs=load(name, 'ROIs');
%         Rs=Rs.ROIs;
%         F=load(name, 'idx_components');
%         F=F.idx_components+1;
%         Rs=Rs(:,F);
%         cor_name=strrep(name,'analysis_matlab','correlation');
%         cor_im=load(cor_name);cor_im=cor_im.Correlation_image;
%         dims=size(cor_im);
%         if file_nb==1
%             temp_roi=0;
%         else
%             temp_roi=temp_roi+size(ROI,3);
%         end
%         ROI=reshape(full(Rs),dims(1),dims(2),size(Rs,2));counter=1;
%         [slice,~]=regexp(name,'_(\d+)_','tokens','match');slice=str2num(slice{1}{1});
%         for roi_nb=1:size(ROI,3)
%             progressbar([],[],roi_nb/size(ROI,3));
%             temp=regionprops(uint16(squeeze(ROI(:,:,roi_nb)))==max(max(uint16(squeeze(ROI(:,:,roi_nb))))),'Centroid');
%             Centroids(roi_nb+temp_roi,5)=roi_nb+temp_roi;
%             %Centroids(roi_nb+temp_roi,1:2)=temp.Centroid; With the resize
%             %need to multiply coords by 1.5
%             Centroids(roi_nb+temp_roi,1:2)=temp.Centroid*1.5;
%             Centroids(roi_nb+temp_roi,3)=slice;
%         end     
%         
%     end    
%     if iscell(Fish_list)
%         image_name=strcat('_ROIsFish',Fish_list{fish_nb},'.csv');
%     else
%         image_name=strcat('_ROIsFish',num2str(Fish_list(fish_nb)),'.csv');
%     end
%     csvwrite(image_name,Centroids);
% end
% clearvars slice fish_nb roi_nb ROI Centroids IndexC file_nb

load('_aMultipower.mat', 'All_ROIs')
load('_aMultipower.mat', 'idx_Plane')

%Creates ROI csv files
Errored_ROI={};
progressbar(0,0,0);
for fish_nb=1:length(Fish_list)
    Centroids=zeros(1,5);
    if iscell(Fish_list)
        IndexC=strfind({MatFiles.name}, Fish_list{fish_nb});
    else
        IndexC=strfind({MatFiles.name}, num2str(Fish_list(fish_nb)));
    end       
    MatFiles_fish = find(not(cellfun('isempty', IndexC)));
    if MatFiles_fish(1)==1
        numbersForROIs=[1 [MatFiles(MatFiles_fish).GoodNumber]];
    else
        numbersForROIs=[MatFiles(MatFiles_fish(1)-1).GoodNumber [MatFiles(MatFiles_fish).GoodNumber]];
    end
    GoodROIs=idxKmeans_final_goodmember(numbersForROIs(1):numbersForROIs(end)-1);
    GoodPlanes=idx_Plane(numbersForROIs(1):numbersForROIs(end)-1);
    ROIs_temp=All_ROIs{MatFiles_fish(1)};
        for i=2:length(MatFiles_fish)
        ROIs_temp=horzcat(ROIs_temp,All_ROIs{MatFiles_fish(i)});
    end    
    progressbar(fish_nb/length(Fish_list));        
    if find(GoodROIs>0)
        counter=1;
        for i=GoodBetas
            idx_temp=find(GoodROIs==i);            
            if idx_temp                
                ROIs_temp2=reshape(full(ROIs_temp(:,idx_temp)),540,640,length(idx_temp));%reshape cause of resize
                for roi=1:size(ROIs_temp2,3)
                    image_temp=squeeze(ROIs_temp2(:,:,roi));image_temp=image_temp/max(image_temp(:));
                    [x,y]=find(image_temp>0.5);
                    Centroids(counter:counter+length(x)-1,1)=y;
                    Centroids(counter:counter+length(x)-1,2)=x;
                    Centroids(counter:counter+length(x)-1,3)=GoodPlanes(idx_temp(roi));
                    Centroids(counter:counter+length(x)-1,5)=idx_temp(roi)+numbersForROIs(1)-1;
                    Centroids(counter:counter+length(x)-1,6)=i;
                    counter=counter+length(x);                    
                end
            end
        end
    end   
    if iscell(Fish_list)
        image_name=strcat('_ROIsFishFullGoodBetas',Fish_list{fish_nb},'.csv');
    else
        image_name=strcat('_ROIsFishFullGoodBetas',num2str(Fish_list(fish_nb)),'.csv');
    end
    csvwrite(image_name,Centroids);
end
clearvars slice fish_nb roi_nb ROI Centroids IndexC file_nb


%---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
% Read the output of ANTs to extract the ROIs centroid
%

CSV_Files=dir('_2Warped*.csv');
%CSV_Files=dir('_Warped*.csv');
ROIs=struct();
for i=1:length(CSV_Files);
    temp=csvread(CSV_Files(i).name,1);
    %Fishname=regexp(CSV_Files(i).name,'Fish(\d+_\w+).csv','tokens');Fishname=Fishname{1}{1};
    %Fishname=regexp(CSV_Files(i).name,'ROIsFish(.+)_Resized.csv','tokens');Fishname=Fishname{1}{1};
    Fishname=regexp(CSV_Files(i).name,'_2Warped(.+)_Resized.csv','tokens');Fishname=Fishname{1}{1};
    %Fishname=regexp(CSV_Files(i).name,'Betas(\d+_\w+).csv','tokens');Fishname=Fishname{1}{1};
    ROIs(i).name=Fishname;    
    ROIs(i).coord=temp(:,1:3);
    ROIs(i).idx=temp(:,5);
    %ROIs(i).cluster=temp(:,6);
end
clearvars i temp CSV_Files Fishname