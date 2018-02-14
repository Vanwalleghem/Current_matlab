%All brain regions
Zbrain_AllMask=Zbrain_Masks{1,3};
for brain_nb=2:length(MaskFiles)   
    Mask=Zbrain_Masks{brain_nb,3};
    temp=unique([Zbrain_AllMask; Mask],'rows');
    Zbrain_AllMask=temp;
end
clearvars temp brain_nb
%Removing the eyes
Mask=Zbrain_Masks{78,3};
IsInEyes_temp=ismember(Zbrain_AllMask,Mask,'rows');IsInEyes_temp=find(IsInEyes_temp==1);
Zbrain_AllMask(IsInEyes_temp,:)=[];

idxKmeans_final_goodmemberInBrain=idxKmeans_final_goodmember;
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
        for i=GoodBetas
            idx_temp=find(GoodROIs==i);
            if idx_temp
                Coords=ROI_fish(idx_temp,:);
                IsInMask_temp=ismember(Coords,Zbrain_AllMask,'rows');
                GoodROIs(idx_temp(find(IsInMask_temp==0)))=0;
            end           
        end
        idxKmeans_final_goodmemberInBrain(numbersForROIs(1):numbersForROIs(end)-1)=GoodROIs;
    end
end
clearvars GoodROIs Mask counter idx_temp Coords IsInMask IsInMask_temp

fname = 'Huc_H2B_RFP_8bit.tif';
info = imfinfo(fname);
num_images = numel(info);
Template=zeros(info(1).Height,info(1).Width,length(info),'uint8');
for k = 1:num_images
    image_temp = imread(fname, k, 'Info', info)';image_temp=double(image_temp);
    image_temp=image_temp/1.5;%prctile(prctile(image_temp,90),90);image_temp=image_temp*100;
    Template(:,:,k) = image_temp';    
end
clearvars info fname num_images k i j 

imageSizeY = size(Template,1);
imageSizeX = size(Template,2);
[columnsInImage rowsInImage] = meshgrid(1:imageSizeX, 1:imageSizeY);radius =3;

Template2=zeros([size(Template),3]);
Template2=repmat(Template,1,1,1,3);
colors=[256 50 0; 240 159 0; 0 158 115; 0 114 250];
for fish_nb=ELO_fish
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
    ROI_fish(:,3)=round(((ROI_fish(:,3)-1)*2)+24);
    if MatFiles_fish(1)==1
        numbersForROIs=[1 [MatFiles(MatFiles_fish).GoodNumber]];
    else
        numbersForROIs=[MatFiles(MatFiles_fish(1)-1).GoodNumber+1 [MatFiles(MatFiles_fish).GoodNumber]];
    end
    GoodROIs=idxKmeans_final_goodmemberInBrain(numbersForROIs(1):numbersForROIs(end)-1);
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
                            if xcoord>1 & ycoord>1
                                circlePixels = (rowsInImage - ycoord).^2 + (columnsInImage - xcoord).^2 <= radius.^2;
                                for col=1:3
                                    image_temp=squeeze(squeeze(Template2(:,:,slice,col)));
                                    image_temp(circlePixels)=colors(counter,col);
                                    Template2(:,:,slice,col)=image_temp;
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

OutputName='Zbrain_Multipower_ELO_highThreshROI_noBrain.tif';
delete(OutputName);
 for slice=1:size(Template2,3)
        image_temp=uint8(squeeze(Template2(:, :,slice,:)));%image_temp=uint8(image_temp);
        %image=(image-min_score)/(max_score-min_score);image=image*256;image=uint16(image);        
        imwrite(image_temp, OutputName, 'WriteMode', 'append');
 end
 
 Template2=zeros([size(Template),3]);
Template2=repmat(Template,1,1,1,3);
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
    ROI_fish=ROIs(ROI_fish).coord;ROI_fish(:,1:2)=round(ROI_fish(:,1:2));
    ROI_fish(:,3)=round(((ROI_fish(:,3)-1)*2)+24);
    if MatFiles_fish(1)==1
        numbersForROIs=[1 [MatFiles(MatFiles_fish).GoodNumber]];
    else
        numbersForROIs=[MatFiles(MatFiles_fish(1)-1).GoodNumber+1 [MatFiles(MatFiles_fish).GoodNumber]];
    end
    GoodROIs=idxKmeans_final_goodmemberInBrain(numbersForROIs(1):numbersForROIs(end)-1);
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
                            if xcoord>1 & ycoord>1
                                circlePixels = (rowsInImage - ycoord).^2 + (columnsInImage - xcoord).^2 <= radius.^2;
                                for col=1:3
                                    image_temp=squeeze(squeeze(Template2(:,:,slice,col)));
                                    image_temp(circlePixels)=colors(counter,col);
                                    Template2(:,:,slice,col)=image_temp;
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

OutputName='Zbrain_Multipower_ERO_highThreshROI_noBrain.tif';
delete(OutputName);
 for slice=1:size(Template2,3)
        image_temp=uint8(squeeze(Template2(:, :,slice,:)));%image_temp=uint8(image_temp);
        %image=(image-min_score)/(max_score-min_score);image=image*256;image=uint16(image);        
        imwrite(image_temp, OutputName, 'WriteMode', 'append');
 end