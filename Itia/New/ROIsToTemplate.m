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

fname = 'Template_Itia_16012018.tif';
info = imfinfo(fname);
num_images = numel(info);
Template=zeros(info(1).Width,info(1).Height,length(info));
for k = 1:num_images
    Template(:,:,k) = imread(fname, k, 'Info', info)';
    % ... Do something with image A ...
end

for fish_nb=1:length(Fish_list)
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
    ROI_fish=ROIs(ROI_fish).coord;ROI_fish=round(ROI_fish);    
    if MatFiles_fish(1)==1
        numbersForROIs=[1 [MatFiles(MatFiles_fish).GoodNumber]];
    else
        numbersForROIs=[MatFiles(MatFiles_fish(1)-1).GoodNumber [MatFiles(MatFiles_fish).GoodNumber]];
    end
    PlanesForROIs=idx_Plane(numbersForROIs(1):numbersForROIs(end)-1);
    GoodROIs=idxKmeans_final_goodmember(numbersForROIs(1):numbersForROIs(end)-1);
    if find(GoodROIs>0)        
        for i=GoodBetas
            idx_temp=find(GoodROIs==i);
            
        end
        
        for slice=1:25
            
            
            
        end
    end
    
    for slice=1:25
        idx_Planes=find(PlanesForROIs==slice);
        if MatFiles_fish(file_nb)==1
            base_number=0;
            number=MatFiles(MatFiles_fish(file_nb)).GoodNumber;
        end
        if ~isempty(file_nb)
            number=MatFiles(MatFiles_fish(file_nb)).GoodNumber;
        end
    end
end





progressbar(0,0,0);
image=
for fish_nb=1:length(Fish_list)
     if iscell(Fish_list)
        Fish_name=Fish_list{fish_nb};        
     else
        Fish_name=num2str(Fish_list(fish_nb));
                
     end
    IndexC=strfind({MatFiles.name}, Fish_name);
    MatFiles_fish = find(not(cellfun('isempty', IndexC)));
    
    
    coord_name=strcat('_ROIsFish',num2str(Fish_list(fish_nb)),'.tif');
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
            if slice==1
                temp_roi=0;
            else
                temp_roi=temp_roi+size(ROI,3);
            end
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
                image(centerX,centerY)=center+temp_roi;
                if ~length(unique(image))-1==size(Rs,2)
                    Errored_ROI{counter}=file_nb;
                    counter=counter+1;
                end
            end
            %figure;imagesc(image)
        end
        imwrite(uint16(image),coord_name,'WriteMode', 'append');
    end
end
