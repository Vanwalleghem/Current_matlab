Twoplanes_transfo=load('D:\Pictures\processed\Flow\Basic\transformation_crop_to_full_2planes.mat');
Errored_ROI={};
progressbar(0,0,0);
for fish_nb=1:2%length(Fish_list)
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
        if findstr(filename,'2planes')
            [slice,~]=regexp(filename,'f\d-(\d+)um_','tokens','match');slice=str2num(slice{1}{1});
            idx_name=strcat(num2str(Fish_list(fish_nb)),'-',num2str(slice));
            Transfo_idx = ~cellfun('isempty',strfind(cellstr(Twoplanes_transfo.FileName),idx_name));Transfo_idx=find(Transfo_idx>0);
        else
            [slice,~]=regexp(filename,'\d+_(\d+)_','tokens','match');slice=str2num(slice{1}{1});
            idx_name=strcat(num2str(Fish_list(fish_nb)),'_',num2str(slice));
            Transfo_idx = ~cellfun('isempty',strfind(cellstr(FileName),strcat(idx_name,'.tif')));Transfo_idx=find(Transfo_idx>0);
        end
        ROI=All_ROIs{MatFiles_fish(plane)};
        imagename=regexp(filename,'_output_analysis','split');
        imagename=strcat('AVG_',imagename{1},'.tif');
        image=double(imread(imagename));image=image/max(max(image));image=image*128;
        ROI=reshape(full(ROI),size(image,1),size(image,2),size(ROI,2));        
        if isempty(Transfo_idx)
            FileName
            break
        end
        if findstr(filename,'2planes')
            Transfo_tmp=Twoplanes_transfo.TransfoMatrix(Transfo_idx,:);
             if Transfo_tmp(1)<-5
                Transfo_tmp(1)=830+Transfo_tmp(1);
            end
            if Transfo_tmp(2)<-5
                Transfo_tmp(2)=1120+Transfo_tmp(2);
            end
        else
            Transfo_tmp=TransfoMatrix(Transfo_idx,:);
            if Transfo_tmp(1)<-5
                Transfo_tmp(1)=1080+Transfo_tmp(1);
            end
            if Transfo_tmp(2)<-5
                Transfo_tmp(2)=1280+Transfo_tmp(2);
            end
        end
        for roi_nb=1:size(ROI,3)
            test=ROI(:,:,roi_nb);
            if fish_nb<3                
                test=padarray(test,[830-size(test,1),1120-size(test,2)],0,'post');
                test=circshift(test,Transfo_tmp);
                test=imrotate(test,90);
            end
            [M I]=max(test(:));
            [I_row I_col]=ind2sub(size(test),I);
            Centroids(counter,5)=counter;
            temp=[I_col I_row];
            if fish_nb>2
                Centroids(counter,1:2)=temp+Transfo_tmp([2 1]);
            else
                Centroids(counter,1:2)=temp;
            end
            if fish_nb==5
                Centroids(counter,3)=((280-slice)/20)+1;
            else
                Centroids(counter,3)=(slice/20)+1;
            end
            progressbar([],[],roi_nb/size(ROI,3));
            counter=counter+1;
        end
        progressbar([],fish_nb/length(MatFiles_fish));
    end
    if iscell(Fish_list)
        image_name=strcat('_ROIsFish_',Fish_list{fish_nb},'b.csv');
    else
        image_name=strcat('_ROIsFish_',num2str(Fish_list(fish_nb)),'b.csv');
    end
    csvwrite(image_name,Centroids);
    progressbar(fish_nb/length(Fish_list));
end
