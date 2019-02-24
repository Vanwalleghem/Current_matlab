Template2=zeros([size(Template),3]);
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
    ROI_cluster=ROIs(ROI_fish).cluster;
    ROI_fish=ROIs(ROI_fish).coord;ROI_fish(:,1:2)=round(ROI_fish(:,1:2));
    ROI_fish(:,3)=round(((ROI_fish(:,3)-1)*2)+24);
    counter=1;
    for i=GoodBetas
        idx_temp=find(ROI_cluster==i);
        if idx_temp
            for slice=1:size(Template2,3)
                idx_slice=find(ROI_fish(idx_temp,3)==slice);
                if idx_slice
                    for roi_nb=1:length(idx_slice)
                        ycoord=ROI_fish(idx_temp(idx_slice(roi_nb)),1);
                        xcoord=ROI_fish(idx_temp(idx_slice(roi_nb)),2);
                        if (( 0 < xcoord & xcoord <= size(Template,1) ) & (0 < ycoord & ycoord <= size(Template,2)))
                            for col=1:3
                                    image_temp=squeeze(squeeze(Template2(:,:,slice,col)));
                                    image_temp(xcoord,ycoord)=colors(counter,col);
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


OutputName='Zbrain_Multipower_ERO_highThreshFullROI_nogray.tif';
delete(OutputName);
 for slice=1:size(Template2,3)
        image_temp=uint8(squeeze(Template2(:, :,slice,:)));%image_temp=uint8(image_temp);
        %image=(image-min_score)/(max_score-min_score);image=image*256;image=uint16(image);        
        imwrite(image_temp, OutputName, 'WriteMode', 'append');
 end
 