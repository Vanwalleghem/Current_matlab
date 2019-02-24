ItiaList={'Thalamus','Cerebellum','NucMLF','Semicircularis','Telencephalon','Tectum','Longitudinalis','Tegmentum','Habenula','Hindbrain','MON','Pretectum'};
i=12;
regionName=ItiaList{i};
progressbar;
for i=12
    progressbar(i/length(ItiaList),[]);
    regionName=ItiaList{i};
    Mask=[];
    IndexC=strfind({Zbrain_Masks{:,2}}, regionName);
    IndexC=find(not(cellfun('isempty', IndexC)));
    if i==5
        IndexC=294;
    end
    for j=IndexC
        if isempty(Mask)
            Mask=Zbrain_Masks{j,3};
        else
            Mask=vertcat(Mask,Zbrain_Masks{j,3});
        end
    end
    for fish_nb=1:length(Fish_list)
        progressbar([],fish_nb/length(Fish_list));
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
        IsInBrainRegion=ismember(ROI_fish,Mask,'rows');
        PerBrainRegions(fish_nb).(regionName).ROIsCent=ROI_fish(IsInBrainRegion,:);        
        temp=ZS(numbersForROIs(1):numbersForROIs(end),:);IsInBrainRegion=temp(IsInBrainRegion,:);        
        PerBrainRegions(fish_nb).(regionName).ZS=IsInBrainRegion;        
    end
end


i=12;
regionName=ItiaList{i};
for j=1:length(PerBrainRegions)
    if j==1;
        ZS_Brain.(regionName)=PerBrainRegions(j).(regionName).ZS;
    else
        ZS_Brain.(regionName)=vertcat(ZS_Brain.(regionName),PerBrainRegions(j).(regionName).ZS);
    end
end

for i=12
    regionName=ItiaList{i};
    for j=1:length(PerBrainRegions)
        if j==1;
            ROIsPerBrain.(regionName).ROIs=PerBrainRegions(j).(regionName).ROIsCent;
            ROIsPerBrain.(regionName).Numbers=length(PerBrainRegions(j).(regionName).ROIsCent);
            temp=length(PerBrainRegions(j).(regionName).ROIsCent);
        else
            ROIsPerBrain.(regionName).ROIs=vertcat(ROIsPerBrain.(regionName).ROIs, PerBrainRegions(j).(regionName).ROIsCent); 
            temp=temp+length(PerBrainRegions(j).(regionName).ROIsCent);
            ROIsPerBrain.(regionName).Numbers=horzcat(ROIsPerBrain.(regionName).Numbers,temp);
        end
    end    
end
clearvars i j k ans Cmap coef coefficients counter Fish_name fish_nb GCaMP6 GoodBetas_nuc idxStart idx_Fish idx_temp idxKmeans IndexC image_temp IsInBrainRegion
clearvars IsInBrainRegion_good Mask MatFiles_fish mdl numbersForROIs_good numbersForROIs regionName ROI_fish ROI_fish_good ROI_name rsquared  temp



for j=12
    regionName=ItiaList{j};
    temp=ZS_Brain.(regionName);	
    coefficients=struct();
    rsquared=zeros(1,length(temp));
    parfor i=1:size(temp,1)
        mdl=fitlm(Stimuli',temp(i,:));%,'interactions');
        coefficients(i).coef=mdl.Coefficients;
        rsquared(i)=mdl.Rsquared.Adjusted;
    end    
    LinReg.(regionName).coef=coefficients;
    LinReg.(regionName).rsquared=rsquared;    
end

%Kmeans of all
options = statset('UseParallel',1); 
for j=12
    regionName=ItiaList{j};
    idx_temp=find(LinReg.(regionName).rsquared>0.1);
    LinReg.(regionName).ZS_rsq=ZS_Brain.(regionName)(idx_temp,:);
    [idxKmeans Cmap]=kmeans(LinReg.(regionName).ZS_rsq,10,'Options',options,'Distance','cityblock','Replicates',3,'MaxIter',1000,'Display','final');
    LinReg.(regionName).KmeansCenter=Cmap;
    LinReg.(regionName).KmeansIdx=idxKmeans;
end

%Get the good ones
for j=12
    regionName=ItiaList{j};
    [~,LinReg.(regionName).GoodBetas]=Test_Regress( LinReg.(regionName).KmeansCenter,Stimuli,LinReg.(regionName).KmeansIdx,0.5);
end

%Correlate back to Kmeans to remove crap
Threshold=0.4;
for i=12
    regionName=ItiaList{i};
    temp=LinReg.(regionName).ZS_rsq;
    idx_temp=LinReg.(regionName).KmeansIdx;
    GoodBet_temp=LinReg.(regionName).GoodBetas;
    for j=1:length(GoodBet_temp)
        idx_g=find(idx_temp==GoodBet_temp(j));
        temp_g=temp(idx_g,:);
        corr_temp=zeros(size(idx_g));
        parfor k=1:length(idx_g)        
            temp=corrcoef(LinReg.(regionName).KmeansCenter(GoodBet_temp(j),:), temp_g(k,:));
            corr_temp(k)=temp(1,2);
        end        
        idx_temp(idx_g(find(corr_temp<=Threshold)))=0;
    end
    LinReg.(regionName).KmeansIdx_select=idx_temp;    
end
clearvars i j k temp ans Cmap coeff coefficients explained i Hindbrain_Mask idx_temp idxKmeans ish_nb 

figure;
counter=1;xplot=length(ItiaList);yplot=4;counter2=1;
for i=12
    regionName=ItiaList{i};
    idx_temp=LinReg.(regionName).KmeansIdx_select;
    GoodBet_temp=LinReg.(regionName).GoodBetas;
    counter=counter2;
    for j=1:length(GoodBet_temp)
        idx_temp2=find(idx_temp==GoodBet_temp(j));
        subplot(xplot,yplot,counter+(j-1));plot(mean(LinReg.(regionName).ZS_rsq(idx_temp2,:),1));
    end
    counter2=counter2+yplot;
end

%Building the AVG across 3 presentations
for j=12
    regionName=ItiaList{j};
    ZS2=LinReg.(regionName).ZS_rsq;
    ZS_AVG2=zeros(size(ZS2,1),246);
    parfor idx_ZS=1:size(ZS2,1)
        start=30;
        AVG=[];
        for i=1:3
            AVG(i,:)=ZS2(idx_ZS,start:start+40);
            start=start+40;
        end
        AVG=mean(AVG,1);
        AVG=AVG-mean(AVG(1:10));%-min(AVG);
        temp=[];
        for j=2:6
            for i=1:3
                temp(i,:)=ZS2(idx_ZS,start:start+40);
                start=start+40;
            end
            temp=mean(temp,1);
            temp=temp-mean(temp(1:10));%-min(temp);
            AVG=[AVG temp];
        end
        ZS_AVG2(idx_ZS,:)=AVG;
    end
    LinReg.(regionName).ZS_AVG=ZS_AVG2;
end

GoodBet_temp=LinReg.(regionName).GoodBetas;
KmeansIdx_merge=LinReg.(regionName).KmeansIdx_select;
idx_temp=ismember(KmeansIdx_merge,GoodBet_temp([2 3]));
KmeansIdx_merge(idx_temp)=GoodBet_temp(3);
LinReg.(regionName).GoodBeta_merge=GoodBet_temp([3 1]);
LinReg.(regionName).KmeansIdx_merge=KmeansIdx_merge;

i=12;colors{i}=[0.14 1 0.14; 0.7 0.4 1];regionName=ItiaList{i};
%Graph figure with thick lines and no SD
counter=1;xplot=3;yplot=1;
x = linspace(0.25,246/4,246);
region_nb=12;
Fighandle=figure;
set(Fighandle, 'Position', [100, 100, 600, 250]);
counter=1;
regionName=ItiaList{region_nb};
idx_temp=LinReg.(regionName).KmeansIdx_merge;
GoodBet_temp=LinReg.(regionName).GoodBeta_merge;
if region_nb==5
    counter=counter+1;
end
for j=1:length(GoodBet_temp)
    idx_temp2=find(idx_temp==GoodBet_temp(j));
    for k=0:5
        meanToPlot=mean(LinReg.(regionName).ZS_AVG(idx_temp2,4+(k*41):40+(k*41)),1);
        hold on;
        plot(x(4+(k*40):40+(k*40)),meanToPlot-mean(meanToPlot(1:4)),'color',colors{region_nb}(j,:),'LineWidth',3);axis([0 60 -3 3]);
        if (j==1 & k<5)
            rectangle('EdgeColor','none','FaceColor',[0.25, 0.25, 0.25, 0.4-k*0.06],'Position',[x(9+(k*40)) -3 1 7]);
        end
    end
    hold on;
end
hold off;

print(Fighandle,regionName,'-dsvg','-r0');
close all;


for i=12
    regionName=ItiaList{i};
    idx_temp=find(LinReg.(regionName).rsquared>0.1);
    ROIs_temp=ROIsPerBrain.(regionName).ROIs(idx_temp,:);
    Numbers=ROIsPerBrain.(regionName).Numbers;Numbers=[0 Numbers];
    for j=1:length(LinReg.(regionName).GoodBeta_merge)
        idx_goodROIs=find(LinReg.(regionName).KmeansIdx_merge==LinReg.(regionName).GoodBeta_merge(j));
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
        ROI_coord_ELO{i,j}=goodROIs(idx_ELO,:);
        ROI_coord_ERO{i,j}=goodROIs(idx_ERO,:);    
    end
end
clearvars goodROIs ROIs_temp regionName idx_temp Template2 idx_slice roi_nb xcoord ycoord circlePixels image_temp columnsInImage rowsInImage col corr_temp Fighdandle framerate GoodBet_temp H idx_g Model_cereb


for region_nb=12
    ForCSVExport=[];
    for cluster=1:3
        if ROI_coord_ERO{region_nb,cluster}
            if ForCSVExport
                CSV_temp=zeros(length(ROI_coord_ERO{region_nb,cluster}),4);
                CSV_temp(:,1:3)=ROI_coord_ERO{region_nb,cluster};
                CSV_temp(:,4)=cluster;
                ForCSVExport=[ForCSVExport;CSV_temp];
            else
                ForCSVExport=zeros(length(ROI_coord_ERO{region_nb,cluster}),4);
                ForCSVExport(:,1:3)=ROI_coord_ERO{region_nb,cluster};
                ForCSVExport(:,4)=cluster;
            end
        end
    end
    ForCSVExport(:,3)=ForCSVExport(:,3)*2;
    filename=strcat('__ROIs_coord',ItiaList{region_nb},'.csv');    
    csvwrite(filename,ForCSVExport);
end

