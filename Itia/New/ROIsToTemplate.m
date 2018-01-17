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

Fish_list_backup=Fish_list;
IndexC=strfind(Fish_list_backup, 'ELO');
ELO_fish=find(not(cellfun('isempty', IndexC)));
ERO_fish=find(cellfun('isempty', IndexC));
Fish_list=Fish_list_backup;

%fname = 'Template_Itia_16012018.tif';
fname = 'Itia_templatetemplate0Fish201712141_ELO_IRO_HS_mean20WarpedToTemplate.tif';
info = imfinfo(fname);
num_images = numel(info);
Template=zeros(info(1).Width,info(1).Height,length(info),'uint8');
for k = 1:num_images
    image_temp = imread(fname, k, 'Info', info)';
    image_temp=image_temp/max(max(image_temp));image_temp=image_temp*128;
    Template(:,:,k) = image_temp;
    % ... Do something with image A ...
end


Template2=repmat(Template,1,1,1,3);
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
    GoodROIs=idxKmeans_final_goodmember(numbersForROIs(1):numbersForROIs(end)-1);
    if find(GoodROIs>0)
        counter=1;
        for i=GoodBetas
            idx_temp=find(GoodROIs==i);
            if idx_temp
                for slice=1:size(Template2,3)
                    idx_slice=find(ROI_fish(idx_temp,3)==slice);
                    if idx_slice
                        for roi_nb=1:length(idx_slice)
                            xcoord=ROI_fish(idx_temp(idx_slice(roi_nb)),1);
                            ycoord=ROI_fish(idx_temp(idx_slice(roi_nb)),2);
                            if xcoord>0 & ycoord>0
                                for col=1:3
                                    Template2(xcoord,ycoord,slice,col)=colors(counter,col);
                                end
                            end
                        end
                    end
                end
            end
            counter=counter+1;
        end        
    end  
end

 for slice=1:size(Template2,3)
        image_temp=squeeze(Template2(:, :,slice,:));%image_temp=uint8(image_temp);
        %image=(image-min_score)/(max_score-min_score);image=image*256;image=uint16(image);
        OutputName='Template_Multipower.tif';
        imwrite(image_temp, OutputName, 'WriteMode', 'append');
 end

 for slice=1:size(Template2,3)
        image_temp=squeeze(Template2(:, :,slice,:));
        imshow(image_temp);
        pause;
end

Template2=repmat(Template,1,1,1,3);
%for fish_nb=1:length(Fish_list)
for fish_nb=ERO_fish    
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
    GoodROIs=idxKmeans_final_goodmember(numbersForROIs(1):numbersForROIs(end)-1);
    ROIs_temp=All_ROIs{MatFiles_fish(1)};
    for i=2:length(MatFiles_fish)
        ROIs_temp=horzcat(ROIs_temp,All_ROIs{MatFiles_fish(i)});
    end     
    if find(GoodROIs>0)
        counter=1;
        for i=GoodBetas
            idx_temp=find(GoodROIs==i);            
            if idx_temp
                ROIs_temp2=reshape(full(ROIs_temp(:,idx_temp)),imageSizeY,imageSizeX,length(idx_temp));
                for slice=1:size(Template2,3)
                    idx_slice=find(ROI_fish(idx_temp,3)==slice);
                    if idx_slice
                        for roi_nb=1:length(idx_slice)
                            ROI_temp=uint16(squeeze(ROIs_temp2(:,:,idx_slice(roi_nb))));
                            temp=regionprops(uint16(squeeze(ROI_temp))==max(max(ROI_temp)),'Centroid');
                            Centroids=temp.Centroid;
                            xcoord=ROI_fish(idx_temp(idx_slice(roi_nb)),1);
                            ycoord=ROI_fish(idx_temp(idx_slice(roi_nb)),2);
                            shift_x=xcoord-Centroids(1);
                            shift_y=ycoord-Centroids(2);
                            ROI_shifted=circshift(ROI_temp,[round(shift_y),round(shift_x)]);
                            ROI_shifted=double(ROI_shifted')/double(max(max(ROI_shifted)));
                            if xcoord>1 & ycoord>1                                
                                for col=1:3
                                    image=squeeze(squeeze(Template2(:,:,slice,col)));
                                    image(ROI_shifted>0.7)=((ROI_shifted(ROI_shifted>0.7)-0.2)/0.8)*colors(counter,col);                                    
                                    Template2(:,:,slice,col)=image;
                                end
                            end
                        end
                    end
                end
            end
            counter=counter+1;
        end        
    end  
    fish_counter=fish_counter+xplot;
end

 for slice=1:size(Template2,3)
        image_temp=squeeze(Template2(:, :,slice,:));%image_temp=uint8(image_temp);
        %image=(image-min_score)/(max_score-min_score);image=image*256;image=uint16(image);
        OutputName='Template_Multipower_ERO_highThreshROI.tif';
        imwrite(image_temp, OutputName, 'WriteMode', 'append');
 end
