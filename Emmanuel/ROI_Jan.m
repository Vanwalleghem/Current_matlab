%Creates ROI csv files
progressbar();
Centroids=zeros(1,3);counter=1;
for file_nb=1:length(MatFiles)
    progressbar(file_nb/length(MatFiles));
    name=MatFiles(file_nb).name;
    Rs=load(name, 'ROIs');
    Rs=Rs.ROIs;
    F=load(name, 'idx_components');%Include all the ROIs
    F=F.idx_components+1;
    Rs=Rs(:,F);
    cor_name=strrep(name,'analysis_matlab','correlation');
    cor_im=load(cor_name);cor_im=cor_im.Correlation_image;
    dims=size(cor_im);
    if file_nb==1
        temp_roi=0;
    else
        temp_roi=temp_roi+size(ROI,3);
    end
    ROI=reshape(full(Rs),dims(1),dims(2),size(Rs,2));
    [slice,~]=regexp(name,'Slice(\d+)_','tokens','match');slice=str2num(slice{1}{1});
    for roi_nb=1:size(ROI,3)
        temp=regionprops(uint16(squeeze(ROI(:,:,roi_nb)))==max(max(uint16(squeeze(ROI(:,:,roi_nb))))),'Centroid');        
        temp=temp.Centroid;
        Centroids(counter,1:2)=temp;
        Centroids(counter,3)=slice;
        counter=counter+1;
    end    
end
