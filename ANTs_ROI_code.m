MatFiles=dir('*analysis_matlab.mat');

idx_Fish=zeros(1,length(MatFiles));
for i=1:length(MatFiles)    
    name=strcat(MatFiles(i).name);
    [Fish,~]=regexp(name,'Fish2017(\d+)_','tokens','match');Fish=str2num(Fish{1}{1});
    idx_Fish(i)=Fish;
end
clearvars i Fish Plane name counter

Fish_list=unique(idx_Fish);

Errored_ROI={};
progressbar(0,0,0);
for fish_nb=1:length(Fish_list)
    IndexC=strfind({MatFiles.name}, num2str(Fish_list(fish_nb)));
    MatFiles_fish = find(not(cellfun('isempty', IndexC)));
    image_name=strcat('_ROIsFish',num2str(Fish_list(fish_nb)),'.tif');
    progressbar(fish_nb/length(Fish_list));
    for slice=1:25
        progressbar([],slice/25);
        slice_nb=strcat('_',num2str(slice),'_');
        IndexC=strfind({MatFiles(MatFiles_fish).name}, slice_nb);
        file_nb=find(not(cellfun('isempty', IndexC)));
        if isempty(file_nb)
            image=zeros(size(cor_im));radius =1;
        else
            name=strcat(MatFiles(MatFiles_fish(file_nb)).name);
            Rs=load(name, 'ROIs');
            Rs=Rs.ROIs;
            F=load(name, 'idx_components');
            F=F.idx_components+1;
            Rs=Rs(:,F);
            cor_name=strrep(name,'analysis_matlab','correlation');
            cor_im=load(cor_name);cor_im=cor_im.Correlation_image;
            dims=size(cor_im);
            ROI=reshape(full(Rs),dims(1),dims(2),size(Rs,2));counter=1;
            Centroids={};
            for roi_nb=1:size(ROI,3)
                temp=regionprops(uint16(squeeze(ROI(:,:,roi_nb)))==max(max(uint16(squeeze(ROI(:,:,roi_nb))))),'Centroid');
                Centroids{roi_nb}=temp.Centroid;
            end
            %imageSizeX = dims(1);
            %imageSizeY = dims(2);
            %[columnsInImage rowsInImage] = meshgrid(1:imageSizeX, 1:imageSizeY);
            image=zeros(size(cor_im));radius =1;
            for center=1:length(Centroids)
                progressbar([],[],center/length(Centroids));
                centerY=round(Centroids{center}(1));centerX=round(Centroids{center}(2));
                %circlePixels = (rowsInImage - centerY).^2 + (columnsInImage - centerX).^2 <= radius.^2;
                %image(circlePixels')=center;
                image(centerX,centerY)=center;
                if ~length(unique(image))-1==size(Rs,2)
                    Errored_ROI{counter}=file_nb;
                    counter=counter+1;
                end
            end
            %figure;imagesc(image)
        end
        imwrite(uint16(image),image_name,'WriteMode', 'append');
    end
end

idx_Fish={};
for i=1:length(MatFiles)    
    name=strcat(MatFiles(i).name);
    [Fish,~]=regexp(name,'Fish2017(\d+_E\DO)_','tokens','match');Fish=Fish{1}{1});
    idx_Fish{i}=Fish;
end
clearvars i Fish Plane name counter

Fish_list=unique(idx_Fish);

Errored_ROI={};
progressbar(0,0,0);
for fish_nb=1:length(Fish_list)
    IndexC=strfind({MatFiles.name}, Fish_list{fish_nb});
    MatFiles_fish = find(not(cellfun('isempty', IndexC)));
    image_name=strcat('_ROIsFish',Fish_list{fish_nb},'.tif');
    progressbar(fish_nb/length(Fish_list));
    for slice=1:25
        progressbar([],slice/25);
        slice_nb=strcat('_',num2str(slice),'_');
        IndexC=strfind({MatFiles(MatFiles_fish).name}, slice_nb);
        file_nb=find(not(cellfun('isempty', IndexC)));
        if isempty(file_nb)
            image=zeros(size(cor_im));radius =1;
        else
            name=strcat(MatFiles(MatFiles_fish(file_nb)).name);
            Rs=load(name, 'ROIs');
            Rs=Rs.ROIs;
            F=load(name, 'idx_components');
            F=F.idx_components+1;
            Rs=Rs(:,F);
            cor_name=strrep(name,'analysis_matlab','correlation');
            cor_im=load(cor_name);cor_im=cor_im.Correlation_image;
            dims=size(cor_im);
            ROI=reshape(full(Rs),dims(1),dims(2),size(Rs,2));counter=1;
            Centroids={};
            for roi_nb=1:size(ROI,3)
                temp=regionprops(uint16(squeeze(ROI(:,:,roi_nb)))==max(max(uint16(squeeze(ROI(:,:,roi_nb))))),'Centroid');
                Centroids{roi_nb}=temp.Centroid;
            end
            %imageSizeX = dims(1);
            %imageSizeY = dims(2);
            %[columnsInImage rowsInImage] = meshgrid(1:imageSizeX, 1:imageSizeY);
            image=zeros(size(cor_im));radius =1;
            for center=1:length(Centroids)
                progressbar([],[],center/length(Centroids));
                centerY=round(Centroids{center}(1));centerX=round(Centroids{center}(2));
                %circlePixels = (rowsInImage - centerY).^2 + (columnsInImage - centerX).^2 <= radius.^2;
                %image(circlePixels')=center;
                image(centerX,centerY)=center;
                if ~length(unique(image))-1==size(Rs,2)
                    Errored_ROI{counter}=file_nb;
                    counter=counter+1;
                end
            end
            %figure;imagesc(image)
        end
        imwrite(uint16(image),image_name,'WriteMode', 'append');
    end
end




% Create a logical image of a circle with specified
% diameter, center, and image size.
% First create the image.
imageSizeX = dims(1);
imageSizeY = dims(2);
[columnsInImage rowsInImage] = meshgrid(1:imageSizeX, 1:imageSizeY);
% Next create the circle in the image.
centerX = 320;
centerY = 240;
radius = 100;
circlePixels = (rowsInImage - centerY).^2 ...
    + (columnsInImage - centerX).^2 <= radius.^2;


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