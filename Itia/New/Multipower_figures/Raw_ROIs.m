PerBrainRegions_raw=struct();
ItiaList={'Thalamus','Cerebellum','NucMLF','Semicircularis','Telencephalon','Tectum','Longitudinalis','Tegmentum','Habenula'};
progressbar;
for i=1:length(ItiaList)
    progressbar(i/length(ItiaList),[]);
    regionName=ItiaList{i};
    Mask=[];
    IndexC=strfind({Zbrain_Masks{:,2}}, regionName);
    IndexC=find(not(cellfun('isempty', IndexC)));
    for j=IndexC
        if isempty(Mask)
            Mask=Zbrain_Masks{j,3};
        else
            Mask=vertcat(Mask,Zbrain_Masks{j,3});
        end
    end
    for fish_nb=1:13
        progressbar([],fish_nb/13);
        if iscell(Fish_list)
            Fish_name=Fish_list{fish_nb};
        else
            Fish_name=num2str(Fish_list(fish_nb));
        end
        IndexC=strfind({MatFiles.name}, Fish_name);
        MatFiles_fish = find(not(cellfun('isempty', IndexC)));
        ROI_name=strsplit(Fish_name,'Fish2017');
        while iscell(ROI_name)
            ROI_name=ROI_name{1};
        end
        IndexC=strfind({ROIs.name},ROI_name);
        ROI_fish=find(not(cellfun('isempty', IndexC)));        
        ROI_fish=ROIs(ROI_fish).coord;ROI_fish(:,1:2)=round(ROI_fish(:,1:2));
        ROI_fish(:,3)=round(((ROI_fish(:,3)-1)*2)+24);        
        if MatFiles_fish(1)==1
            numbersForROIs=[1 [MatFiles(MatFiles_fish).number]];
        else
            numbersForROIs=[MatFiles(MatFiles_fish(1)-1).number+1 [MatFiles(MatFiles_fish).number]];
        end
        IsInBrainRegion=[];        
        IsInBrainRegion=ismember(ROI_fish,Mask,'rows');%IsInBrainRegion_good=ismember(ROI_fish_good,Mask,'rows');
        ROI_fish=find(not(cellfun('isempty', IndexC)));        
        ROI_fish=ROIs(ROI_fish).coord;ROI_fish(:,1:2)=ROI_fish(:,1:2);
        ROI_fish(:,3)=((ROI_fish(:,3)-1)*2)+24;
        PerBrainRegions_raw(fish_nb).(regionName).ROIsCent=ROI_fish(IsInBrainRegion,:);
    end
end

Hindbrain_Mask=Zbrain_Masks{259,3};
Mask=Zbrain_Masks{131,3};
IsInEyes_temp=ismember(Hindbrain_Mask,Mask,'rows');IsInEyes_temp=find(IsInEyes_temp==1);
Hindbrain_Mask(IsInEyes_temp,:)=[];
clearvars i j fish_nb Mask

regionName='Hindbrain';
progressbar;
 for fish_nb=1:13
        progressbar(fish_nb/13);        
        if iscell(Fish_list)
            Fish_name=Fish_list{fish_nb};
        else
            Fish_name=num2str(Fish_list(fish_nb));
        end
        IndexC=strfind({MatFiles.name}, Fish_name);
        MatFiles_fish = find(not(cellfun('isempty', IndexC)));
        ROI_name=strsplit(Fish_name,'Fish2017');
        while iscell(ROI_name)
            ROI_name=ROI_name{1};
        end
        IndexC=strfind({ROIs.name},ROI_name);
        ROI_fish=find(not(cellfun('isempty', IndexC)));        
        ROI_fish=ROIs(ROI_fish).coord;ROI_fish(:,1:2)=round(ROI_fish(:,1:2));
        ROI_fish(:,3)=round(((ROI_fish(:,3)-1)*2)+24);        
        if MatFiles_fish(1)==1           
            numbersForROIs=[1 [MatFiles(MatFiles_fish).number]];
        else            
            numbersForROIs=[MatFiles(MatFiles_fish(1)-1).number+1 [MatFiles(MatFiles_fish).number]];
        end
        IsInBrainRegion=[];
        IsInBrainRegion=ismember(ROI_fish,Hindbrain_Mask,'rows');
        ROI_fish=find(not(cellfun('isempty', IndexC)));        
        ROI_fish=ROIs(ROI_fish).coord;ROI_fish(:,1:2)=ROI_fish(:,1:2);
        ROI_fish(:,3)=((ROI_fish(:,3)-1)*2)+24;
        PerBrainRegions_raw(fish_nb).(regionName).ROIsCent=ROI_fish(IsInBrainRegion,:);
 end     
 
 ROIsPerBrain_raw=struct();
for i=1:length(ItiaList)
    regionName=ItiaList{i};
    for j=1:length(PerBrainRegions)
        if j==1;
            ROIsPerBrain_raw.(regionName).ROIs=PerBrainRegions_raw(j).(regionName).ROIsCent;
            ROIsPerBrain_raw.(regionName).Numbers=length(PerBrainRegions_raw(j).(regionName).ROIsCent);
            temp=length(PerBrainRegions(j).(regionName).ROIsCent);
        else
            ROIsPerBrain_raw.(regionName).ROIs=vertcat(ROIsPerBrain_raw.(regionName).ROIs, PerBrainRegions_raw(j).(regionName).ROIsCent); 
            temp=temp+length(PerBrainRegions_raw(j).(regionName).ROIsCent);
            ROIsPerBrain_raw.(regionName).Numbers=horzcat(ROIsPerBrain_raw.(regionName).Numbers,temp);
        end
    end    
end
regionName='Hindbrain';
for j=1:length(PerBrainRegions)
    if j==1;
        ROIsPerBrain_raw.(regionName).ROIs=PerBrainRegions_raw(j).(regionName).ROIsCent;
        ROIsPerBrain_raw.(regionName).Numbers=length(PerBrainRegions_raw(j).(regionName).ROIsCent);
        temp=length(PerBrainRegions(j).(regionName).ROIsCent);
    else
        ROIsPerBrain_raw.(regionName).ROIs=vertcat(ROIsPerBrain_raw.(regionName).ROIs, PerBrainRegions_raw(j).(regionName).ROIsCent);
        temp=temp+length(PerBrainRegions_raw(j).(regionName).ROIsCent);
        ROIsPerBrain_raw.(regionName).Numbers=horzcat(ROIsPerBrain_raw.(regionName).Numbers,temp);
    end
end

clearvars i j k ans Cmap coef coefficients counter Fish_name fish_nb GCaMP6 GoodBetas_nuc idxStart idx_Fish idx_temp idxKmeans IndexC image_temp IsInBrainRegion
clearvars IsInBrainRegion_good Mask MatFiles_fish mdl numbersForROIs_good numbersForROIs regionName ROI_fish ROI_fish_good ROI_name rsquared  temp

 ItiaList{10}='Hindbrain';
 for i=1:length(ItiaList)
    regionName=ItiaList{i};
    idx_temp=find(LinReg.(regionName).rsquared>0.1);
    ROIs_temp=ROIsPerBrain_raw.(regionName).ROIs(idx_temp,:);
    Numbers=ROIsPerBrain_raw.(regionName).Numbers;Numbers=[0 Numbers];
    for j=1:length(LinReg.(regionName).GoodBetas_AVG_final)
        idx_goodROIs=find(LinReg.(regionName).KmeansIdx_select_AVG==LinReg.(regionName).GoodBetas_AVG_final(j));
        idx_goodROIs_correct=idx_temp(idx_goodROIs);
        idx_ELO=[];
        idx_ERO=[];
        for fishNb=1:length(Fish_list)
            if Index_ELO{fishNb}                
                idx_ELOtemp=find(Numbers(fishNb) < idx_goodROIs_correct & idx_goodROIs_correct <= Numbers(fishNb+1));
                if isempty(idx_ELO)
                    idx_ELO=idx_ELOtemp;
                else
                    idx_ELO=horzcat(idx_ELO, idx_ELOtemp);
                end
            else
                idx_EROtemp=find(Numbers(fishNb) < idx_goodROIs_correct & idx_goodROIs_correct <= Numbers(fishNb+1));
                if isempty(idx_ERO)
                    idx_ERO=idx_EROtemp;
                else
                    idx_ERO=horzcat(idx_ERO, idx_EROtemp);
                end
            end
        end
        goodROIs=ROIs_temp(idx_goodROIs,:); 
        ROI_coord_ELO_raw{i,j}=goodROIs(idx_ELO,:);
        ROI_coord_ERO_raw{i,j}=goodROIs(idx_ERO,:);    
    end
end
clearvars goodROIs ROIs_temp regionName idx_temp Template2 idx_slice roi_nb xcoord ycoord circlePixels image_temp columnsInImage rowsInImage col corr_temp Fighdandle framerate GoodBet_temp H idx_g Model_cereb


for region_nb=1:length(ItiaList)
    ForCSVExport=[];
    for cluster=1:3
        if ROI_coord_ERO{region_nb,cluster}
            if ForCSVExport
                CSV_temp=zeros(length(ROI_coord_ERO_raw{region_nb,cluster}),4);
                CSV_temp(:,1:3)=ROI_coord_ERO_raw{region_nb,cluster};
                CSV_temp(:,4)=cluster;
                ForCSVExport=[ForCSVExport;CSV_temp];
            else
                ForCSVExport=zeros(length(ROI_coord_ERO_raw{region_nb,cluster}),4);
                ForCSVExport(:,1:3)=ROI_coord_ERO_raw{region_nb,cluster};
                ForCSVExport(:,4)=cluster;
            end
        end
    end
    ForCSVExport(:,3)=ForCSVExport(:,3)*2;
    filename=strcat(ItiaList{region_nb},'.csv');    
    csvwrite(filename,ForCSVExport);
end

%DT vs VT
region_nb=1;
ForCSVExport=[];
ForCSVExportb=[];
for cluster=1:3
    if ROI_coord_ERO{region_nb,cluster}
        if ForCSVExport
            CSV_temp=zeros(length(ROI_coord_ERO_raw{region_nb,cluster}),4);
            CSV_temp(:,1:3)=ROI_coord_ERO_raw{region_nb,cluster};
            CSV_temp(:,4)=cluster;
            temp=ROI_coord_ERO{region_nb,cluster};
            ForCSVExport=[ForCSVExport;CSV_temp];
            ForCSVExportb=[ForCSVExportb;temp];
        else
            ForCSVExport=zeros(length(ROI_coord_ERO_raw{region_nb,cluster}),4);
            ForCSVExport(:,1:3)=ROI_coord_ERO_raw{region_nb,cluster};
            ForCSVExport(:,4)=cluster;
            ForCSVExportb=ROI_coord_ERO{region_nb,cluster};
        end
    end
end
regionName=ItiaList{1};
IndexC=strfind({Zbrain_Masks{:,2}}, regionName);
IndexC=find(not(cellfun('isempty', IndexC)));
dors_ROIs=Zbrain_Masks{IndexC(1),3};
IsInDorsal=ismember(ForCSVExportb,dors_ROIs,'rows');
idx_dorsal=find(IsInDorsal==1);
idx_ventral=find(IsInDorsal==0);
ForCSVExport(:,3)=ForCSVExport(:,3)*2;
filename=strcat('ROIs_raw_dorsal',ItiaList{region_nb},'.csv');
csvwrite(filename,ForCSVExport(idx_dorsal,:));
filename=strcat('ROIs_raw_ventral',ItiaList{region_nb},'.csv');
csvwrite(filename,ForCSVExport(idx_ventral,:));
clearvars goodROIs ForCSVExportb temp ROIs_temp regionName idx_temp Template2 idx_slice roi_nb xcoord ycoord circlePixels image_temp columnsInImage rowsInImage col corr_temp Fighdandle framerate GoodBet_temp H idx_g Model_cereb


%Separate PVL neuropil
region_nb=6;
ForCSVExport=[];
ForCSVExportb=[];
for cluster=1:3
    if ROI_coord_ERO{region_nb,cluster}
        if ForCSVExport
            CSV_temp=zeros(length(ROI_coord_ERO_raw{region_nb,cluster}),4);
            CSV_temp(:,1:3)=ROI_coord_ERO_raw{region_nb,cluster};
            CSV_temp(:,4)=cluster;
            temp=ROI_coord_ERO{region_nb,cluster};
            ForCSVExport=[ForCSVExport;CSV_temp];
            ForCSVExportb=[ForCSVExportb;temp];
        else
            ForCSVExport=zeros(length(ROI_coord_ERO_raw{region_nb,cluster}),4);
            ForCSVExport(:,1:3)=ROI_coord_ERO_raw{region_nb,cluster};
            ForCSVExport(:,4)=cluster;
            ForCSVExportb=ROI_coord_ERO{region_nb,cluster};
        end
    end
end
regionName=ItiaList{region_nb};
IndexC=strfind({Zbrain_Masks{:,2}}, regionName);
IndexC=find(not(cellfun('isempty', IndexC)));
PVL=Zbrain_Masks{IndexC(1),3};
IsInPVL=ismember(ForCSVExportb(:,1:3),PVL,'rows');
idx_PVL=find(IsInPVL==1);
idx_neuropil=find(IsInPVL==0);
ForCSVExport(:,3)=ForCSVExport(:,3)*2;
filename=strcat('ROIs_raw_PVL',ItiaList{region_nb},'.csv');
csvwrite(filename,ForCSVExport(idx_PVL,:));
filename=strcat('ROIs_raw_neuropil',ItiaList{region_nb},'.csv');
csvwrite(filename,ForCSVExport(idx_neuropil,:));
clearvars goodROIs ForCSVExportb temp ROIs_temp regionName idx_temp Template2 idx_slice roi_nb xcoord ycoord circlePixels image_temp columnsInImage rowsInImage col corr_temp Fighdandle framerate GoodBet_temp H idx_g Model_cereb

