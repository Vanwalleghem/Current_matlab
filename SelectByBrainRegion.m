MatFiles=dir('*analysis_matlab.mat');
name=strcat(MatFiles(1).name);
Calcium=load(name, 'DenoisedTraces');
Calcium=Calcium.DenoisedTraces;
MatFiles(1).number=size(Calcium,1);
%Spikes=load(name, 'Spikes');
%Spikes=Spikes.Spikes;
Noise=load(name, 'Noise');
Noise=Noise.Noise;
Fitness=load(name, 'idx_components');
Fitness=Fitness.idx_components+1;
GoodCalcium=Calcium(Fitness,:);
GoodNoise=Noise(Fitness,:);
MatFiles(1).GoodNumber=length(Fitness);
%MatFiles(1).GC=GoodCalcium;
for i = 2:length(MatFiles)
    name=strcat(MatFiles(i).name);
    C=load(name, 'DenoisedTraces');
    C=C.DenoisedTraces;
%     if i==3
%         C=[C(:,1) C(:,1) C(:,1:58)];
%     end
    %S=load(name, 'Spikes');
    %S=S.Spikes;
    N=load(name, 'Noise');
    N=N.Noise;
    F=load(name, 'idx_components');
    F=F.idx_components+1;
    GC=C(F,:);
    GN=N(F,:);
    %GS=S(F,:);
    Noise=vertcat(Noise,N);
    Calcium=vertcat(Calcium,C);
    %Spikes=vertcat(Spikes,S);
    Fitness=horzcat(Fitness,F);
    GoodCalcium=vertcat(GoodCalcium,GC);
    GoodNoise=vertcat(GoodNoise,GN);
    %GoodSpikes=vertcat(GoodSpikes,GS);
    MatFiles(i).number=size(Calcium,1);
    MatFiles(i).GoodNumber=MatFiles(i-1).GoodNumber+length(F);
    %MatFiles(i).GC=GC;
end
clearvars GC C S F N name i GS Fitness
ZS_all=zscore(Calcium+Noise,1,2);
ZS2=zscore(GoodCalcium+GoodNoise,1,2);
clearvars GoodCalcium GoodNoise Calcium Noise

IndexC=strfind({MatFiles.name}, Fish_list{1});
MatFiles_fish = find(not(cellfun('isempty', IndexC)));
if MatFiles_fish(1)==1
        numbersForROIs=[1 [MatFiles(MatFiles_fish).GoodNumber]];
    else
        numbersForROIs=[MatFiles(MatFiles_fish(1)-1).GoodNumber [MatFiles(MatFiles_fish).GoodNumber]];
end    
ZS_temp=ZS2(numbersForROIs(1):numbersForROIs(end)-1,:);
ZS_temp=ZS_temp(:,21:end);
ZS_temp(:,741:760)=ZS_temp(:,end-19:end);
ZS2(numbersForROIs(1):numbersForROIs(end)-1,:)=ZS_temp;

IndexC=strfind({MatFiles.name}, Fish_list{2});
MatFiles_fish = find(not(cellfun('isempty', IndexC)));
if MatFiles_fish(1)==1
        numbersForROIs=[1 [MatFiles(MatFiles_fish).GoodNumber]];
    else
        numbersForROIs=[MatFiles(MatFiles_fish(1)-1).GoodNumber [MatFiles(MatFiles_fish).GoodNumber]];
end    
ZS_temp=ZS2(numbersForROIs(1):numbersForROIs(end)-1,:);
ZS_temp=ZS_temp(:,21:end);
ZS_temp(:,741:760)=ZS_temp(:,end-19:end);
ZS2(numbersForROIs(1):numbersForROIs(end)-1,:)=ZS_temp;


PerBrainRegions=struct();
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
        %ROI_fish_good=ROIs_good(ROI_fish).coord;ROI_fish_good(:,1:2)=round(ROI_fish_good(:,1:2));
        %ROI_fish_good(:,3)=round(((ROI_fish_good(:,3)-1)*2)+24);
        ROI_fish=ROIs(ROI_fish).coord;ROI_fish(:,1:2)=round(ROI_fish(:,1:2));
        ROI_fish(:,3)=round(((ROI_fish(:,3)-1)*2)+24);
        if MatFiles_fish(1)==1
            % numbersForROIs_good=[1 [MatFiles(MatFiles_fish).GoodNumber]];
            numbersForROIs=[1 [MatFiles(MatFiles_fish).number]];
        else
            %numbersForROIs_good=[MatFiles(MatFiles_fish(1)-1).GoodNumber+1 [MatFiles(MatFiles_fish).GoodNumber]];
            numbersForROIs=[MatFiles(MatFiles_fish(1)-1).number+1 [MatFiles(MatFiles_fish).number]];
        end
        IsInBrainRegion=[];        
        IsInBrainRegion=ismember(ROI_fish,Mask,'rows');%IsInBrainRegion_good=ismember(ROI_fish_good,Mask,'rows');
        PerBrainRegions(fish_nb).(regionName).ROIsCent=ROI_fish(IsInBrainRegion,:);
        %PerBrainRegions(fish_nb).(regionName).ROIsCent_good=ROI_fish_good(IsInBrainRegion_good,:);
        temp=ZS(numbersForROIs(1):numbersForROIs(end),:);IsInBrainRegion=temp(IsInBrainRegion,:);
        %temp=ZS2(numbersForROIs_good(1):numbersForROIs_good(end),:);IsInBrainRegion_good=temp(IsInBrainRegion_good,:);
        PerBrainRegions(fish_nb).(regionName).ZS=IsInBrainRegion;
        %PerBrainRegions(fish_nb).(regionName).ZS_good=IsInBrainRegion_good;
    end
end

Hindbrain_Mask=Zbrain_Masks{259,3};
Mask=Zbrain_Masks{131,3};
IsInEyes_temp=ismember(Hindbrain_Mask,Mask,'rows');IsInEyes_temp=find(IsInEyes_temp==1);
Hindbrain_Mask(IsInEyes_temp,:)=[];
clearvars i j fish_nb Mask

regionName='Hindbrain';
progressbar;
 for fish_nb=1:length(ItiaList)
        progressbar(fish_nb/length(ItiaList));        
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
        PerBrainRegions(fish_nb).(regionName).ROIsCent=ROI_fish(IsInBrainRegion,:);
        temp=ZS(numbersForROIs(1):numbersForROIs(end),:);IsInBrainRegion=temp(IsInBrainRegion,:);
        PerBrainRegions(fish_nb).(regionName).ZS=IsInBrainRegion;
 end     
        
ItiaList={'Thalamus','Cerebellum','NucMLF','Semicircularis','Telencephalon','Tectum','Longitudinalis','Tegmentum','Habenula','Hindbrain'}; 
     
%Now to make it into one big pool of data
ZS_Brain=struct();
%ZS_Brain_good=struct();
for i=1:length(ItiaList)
    regionName=ItiaList{i};        
    for j=1:length(PerBrainRegions)
        if j==1;
            ZS_Brain.(regionName)=PerBrainRegions(j).(regionName).ZS;
            %ZS_Brain_good.(regionName)=PerBrainRegions(j).(regionName).ZS_good;
        else
            ZS_Brain.(regionName)=vertcat(ZS_Brain.(regionName),PerBrainRegions(j).(regionName).ZS);
            %ZS_Brain_good.(regionName)=vertcat(ZS_Brain_good.(regionName), PerBrainRegions(j).(regionName).ZS_good);        
        end
    end
end

%LinReg of all
LinReg=struct();
%LinReg_good=struct();
for j=1:length(ItiaList)
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
    %temp=ZS_Brain_good.(regionName);	
    %coefficients=struct();
    %rsquared=zeros(1,length(temp));
    %parfor i=1:size(temp,1)
       % mdl=fitlm(Stimuli',temp(i,:));%,'interactions');
        %coefficients(i).coef=mdl.Coefficients;
        %rsquared(i)=mdl.Rsquared.Adjusted;
    %end    
    %LinReg_good.(regionName).coef=coefficients;
    %LinReg_good.(regionName).rsquared=rsquared;    
end

%Kmeans of all
options = statset('UseParallel',1); 
for j=1:length(ItiaList)
    regionName=ItiaList{j};
    idx_temp=find(LinReg.(regionName).rsquared>0.1);
    LinReg.(regionName).ZS_rsq=ZS_Brain.(regionName)(idx_temp,:);
    [idxKmeans Cmap]=kmeans(LinReg.(regionName).ZS_rsq,5,'Options',options,'Distance','cityblock','Replicates',3,'MaxIter',1000,'Display','final');
    LinReg.(regionName).KmeansCenter=Cmap;
    LinReg.(regionName).KmeansIdx=idxKmeans;
%     idx_temp=find(LinReg_good.(regionName).rsquared>0.1);
%     LinReg_good.(regionName).ZS_rsq=ZS_Brain_good.(regionName)(idx_temp,:);
%     [idxKmeans Cmap]=kmeans(LinReg_good.(regionName).ZS_rsq,5,'Options',options,'Distance','cityblock','Replicates',3,'MaxIter',1000,'Display','final');
%     LinReg_good.(regionName).KmeansCenter=Cmap;
%     LinReg_good.(regionName).KmeansIdx=idxKmeans;
end

%Get the good ones
for j=1:length(ItiaList)
    regionName=ItiaList{j};
    [~,LinReg.(regionName).GoodBetas]=Test_Regress( LinReg.(regionName).KmeansCenter,Stimuli,LinReg.(regionName).KmeansIdx,0.5);
end

%NucMLF is noisy
regionName=ItiaList{3};
[idxKmeans Cmap]=kmeans(LinReg.(regionName).ZS_rsq,5,'Options',options,'Distance','cityblock','Replicates',3,'MaxIter',1000,'Display','final');
LinReg.(regionName).KmeansCenter=Cmap;
LinReg.(regionName).KmeansIdx=idxKmeans;
[temp,LinReg.(regionName).GoodBetas]=Test_Regress( LinReg.(regionName).KmeansCenter,Stimuli,LinReg.(regionName).KmeansIdx,0.335);

%Correlate back to Kmeans to remove crap
Threshold=0.4;
for i=1:length(ItiaList)
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
for i=1:length(ItiaList)
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

colors={};
colors{1}=[1 0 0; 0 0 1; 0 1 0];
colors{2}=[0.7 0 0.7; 0 1 0; 0 0 1; 1 0 0];
colors{3}=[1 0 0; 0 1 0];
colors{4}=[0 0 1; 0 1 0];
colors{5}=[0 0 1; 0 1 0; 1 0 0];
colors{6}=[0 1 0; 1 0 0];


ROIsPerBrain=struct();
for i=1:length(ItiaList)
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


%Make pretty images
fname = 'Huc_H2B_RFP.tif';
info = imfinfo(fname);
num_images = numel(info);
Template=zeros(info(1).Height,info(1).Width,length(info),'uint8');
for k = 1:num_images
    image_temp = imread(fname, k, 'Info', info)';image_temp=double(image_temp);
    image_temp=image_temp/max(max(image_temp));image_temp=image_temp*128;
    Template(:,:,k) = image_temp';    
end
clearvars info fname num_images k i j 


imageSizeX = size(Template,1);
imageSizeY = size(Template,2);
[columnsInImage rowsInImage] = meshgrid(1:imageSizeX, 1:imageSizeY);radius =3;
for i=1:length(ItiaList)
    Template2=zeros([size(Template),3]);
    regionName=ItiaList{i};
    idx_temp=find(LinReg.(regionName).rsquared>0.1);
    ROIs_temp=ROIsPerBrain.(regionName).ROIs(idx_temp,:);
    for j=1:length(LinReg.(regionName).GoodBetas)
        goodROIs=find(LinReg.(regionName).KmeansIdx==LinReg.(regionName).GoodBetas(j));
        goodROIs=ROIs_temp(goodROIs,:);
        counter=1;
        for slice=1:size(Template2,3)
            idx_slice=find(goodROIs(:,3)==slice);
            if idx_slice
                for roi_nb=1:length(idx_slice)
                    xcoord=goodROIs(idx_slice(roi_nb),1);
                    ycoord=goodROIs(idx_slice(roi_nb),2);
                    if xcoord>1 & ycoord>1
                        circlePixels = (rowsInImage - ycoord).^2 + (columnsInImage - xcoord).^2 <= radius.^2;
                        for col=1:3
                            image_temp=squeeze(squeeze(Template2(:,:,slice,col)));
                            image_temp(circlePixels)=colors{i}(counter,col);
                            Template2(:,:,slice,col)=image_temp;
                        end
                    end
                end
            end
        end
        counter=counter+1;
    end
    ROIsInBrain{i}=Template2;
end

colors={};
colors{1}=[1 0 0; 0 0 1; 0 1 0];
colors{2}=[0.7 0 0.7; 0 1 0; 0 0 1; 1 0 0];
colors{3}=[1 0 0; 0 1 0];
colors{4}=[0 0 1; 0 1 0];
colors{5}=[0 0 1; 0 1 0; 1 0 0];
colors{6}=[0 1 0; 1 0 0];
for i=1:length(ItiaList)
    Template2=zeros([size(Template),3]);
    regionName=ItiaList{i};
    idx_temp=find(LinReg.(regionName).rsquared>0.1);
    ROIs_temp=ROIsPerBrain.(regionName).ROIs(idx_temp,:);
    for j=1:length(LinReg.(regionName).GoodBetas)
        goodROIs=find(LinReg.(regionName).KmeansIdx_select==LinReg.(regionName).GoodBetas(j));
        goodROIs=ROIs_temp(goodROIs,:);        
        for slice=1:size(Template2,3)
            idx_slice=find(goodROIs(:,3)==slice);
            if idx_slice
                for roi_nb=1:length(idx_slice)
                    xcoord=goodROIs(idx_slice(roi_nb),1);
                    ycoord=goodROIs(idx_slice(roi_nb),2);
                    if xcoord>1 & ycoord>1
                        circlePixels = (rowsInImage - ycoord).^2 + (columnsInImage - xcoord).^2 <= radius.^2;
                        for col=1:3
                            image_temp=squeeze(squeeze(Template2(:,:,slice,col)));
                            image_temp(circlePixels)=colors{i}(j,col)*256;
                            Template2(:,:,slice,col)=image_temp;
                        end
                    end
                end
            end
        end        
    end
    ROIsInBrain_select{i}=Template2;
end
clearvars goodROIs ROIs_temp regionName idx_temp Template2 idx_slice roi_nb xcoord ycoord circlePixels image_temp columnsInImage rowsInImage col corr_temp Fighdandle framerate GoodBet_temp H idx_g Model_cereb

imageSizeY = size(Template,1);
imageSizeX = size(Template,2);
[columnsInImage rowsInImage] = meshgrid(1:imageSizeX, 1:imageSizeY);radius =3;


framerate=4;
x = linspace(1/framerate,size(ZS,2)/framerate,size(ZS,2));
for i=4:length(ItiaList)
    regionName=ItiaList{i};
    OutputName=strcat('Zbrain_Multipower_Both_HighThresh',regionName,'.tif');    
    delete(OutputName);
    Template2=ROIsInBrain_select{i};
    for slice=1:size(Template2,3)
        image_temp=uint8(squeeze(Template2(:, :,slice,:)));
        imwrite(image_temp, OutputName, 'WriteMode', 'append');
    end
%     Fighandle=figure;
%     set(Fighandle, 'Position', [100, 100, 1300, 900]);
%     xplot=length(LinReg.(regionName).GoodBetas);
%     for j=1:length(LinReg.(regionName).GoodBetas)
%         goodROIs=find(LinReg.(regionName).KmeansIdx_select==LinReg.(regionName).GoodBetas(j));
%         ZS_temp=LinReg.(regionName).ZS_rsq;
%         subplot(xplot,1,j);
%         H=shadedErrorBar(x, mean(ZS_temp(goodROIs,:),1), std(ZS_temp(goodROIs,:),1,1));axis([0 190 -3 6]);
%         H.mainLine.Color=colors{i}(j,:);
%         H.patch.FaceColor=colors{i}(j,:);
%         H.edge(1).Color=colors{i}(j,:);
%         H.edge(2).Color=colors{i}(j,:);
%     end
%     print(Fighandle,'-dpng','-r300',strcat('ColoredClusters_',regionName))
end
clearvars OutputName radius slice Template Threshold xplot ZS_temp
clearvars goodROIs ROIs_temp regionName idx_temp Template2 idx_slice roi_nb xcoord ycoord circlePixels image_temp columnsInImage rowsInImage col corr_temp Fighdandle framerate GoodBet_temp H idx_g Model_cereb

%Now to do ELO vs ERO

Index_ELO=strfind(Fish_list, 'ELO');
ELO_fish=find(not(cellfun('isempty', Index_ELO)));
ERO_fish=find(cellfun('isempty', Index_ELO));
imageSizeY = size(Template,1);
imageSizeX = size(Template,2);
[columnsInImage rowsInImage] = meshgrid(1:imageSizeX, 1:imageSizeY);radius =3;

for i=1:length(ItiaList)
    Template_ELO=zeros([size(Template),3]);
    Template_ERO=zeros([size(Template),3]);
    regionName=ItiaList{i};
    idx_temp=find(LinReg.(regionName).rsquared>0.1);
    ROIs_temp=ROIsPerBrain.(regionName).ROIs(idx_temp,:);
    Numbers=ROIsPerBrain.(regionName).Numbers;Numbers=[0 Numbers];
    for j=1:length(LinReg.(regionName).GoodBetas)
        idx_goodROIs=find(LinReg.(regionName).KmeansIdx_select==LinReg.(regionName).GoodBetas(j));
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
        for slice=1:size(Template_ELO,3)
            idx_slice=find(goodROIs(:,3)==slice);
            if idx_slice
                for roi_nb=1:length(idx_slice)
                    xcoord=goodROIs(idx_slice(roi_nb),1);
                    ycoord=goodROIs(idx_slice(roi_nb),2);
                    if xcoord>1 & ycoord>1
                        circlePixels = (rowsInImage - ycoord).^2 + (columnsInImage - xcoord).^2 <= radius.^2;
                        for col=1:3                           
                            if ismember(idx_slice(roi_nb),idx_ELO)
                                image_temp=squeeze(squeeze(Template_ELO(:,:,slice,col)));
                                image_temp(circlePixels)=colors{i}(j,col)*256;
                                Template_ELO(:,:,slice,col)=image_temp;
                            elseif ismember(idx_slice(roi_nb),idx_ERO)
                                image_temp=squeeze(squeeze(Template_ERO(:,:,slice,col)));
                                image_temp(circlePixels)=colors{i}(j,col)*256;
                                Template_ERO(:,:,slice,col)=image_temp;
                            else
                                error_roi=[idx_slice roi_nb];
                                break;
                            end
                        end
                    end
                end
            end
        end        
    end
    InBrain_select_ELO{i}=Template_ELO;
    InBrain_select_ERO{i}=Template_ERO;
end
clearvars goodROIs ROIs_temp regionName idx_temp Template2 idx_slice roi_nb xcoord ycoord circlePixels image_temp columnsInImage rowsInImage col corr_temp Fighdandle framerate GoodBet_temp H idx_g Model_cereb

for i=1:length(ItiaList)
    regionName=ItiaList{i};
    OutputName=strcat('Zbrain_Multipower_ELO_HighThresh',regionName,'.tif');    
    delete(OutputName);
    Template2=InBrain_select_ELO{i};
    for slice=1:size(Template2,3)
        image_temp=uint8(squeeze(Template2(:, :,slice,:)));
        imwrite(image_temp, OutputName, 'WriteMode', 'append');
    end
end
clearvars OutputName radius slice Template Threshold xplot ZS_temp
clearvars goodROIs ROIs_temp regionName idx_temp Template2 idx_slice roi_nb xcoord ycoord circlePixels image_temp columnsInImage rowsInImage col corr_temp Fighdandle framerate GoodBet_temp H idx_g Model_cereb

for i=1:length(ItiaList)
    regionName=ItiaList{i};
    OutputName=strcat('Zbrain_Multipower_ERO_HighThresh',regionName,'.tif');    
    delete(OutputName);
    Template2=InBrain_select_ERO{i};
    for slice=1:size(Template2,3)
        image_temp=uint8(squeeze(Template2(:, :,slice,:)));
        imwrite(image_temp, OutputName, 'WriteMode', 'append');
    end
end
clearvars OutputName radius slice Template Threshold xplot ZS_temp
clearvars goodROIs ROIs_temp regionName idx_temp Template2 idx_slice roi_nb xcoord ycoord circlePixels image_temp columnsInImage rowsInImage col corr_temp Fighdandle framerate GoodBet_temp H idx_g Model_cereb


%Defining laterality
fname = 'Huc_H2B_RFP_8bit.tif';
info = imfinfo(fname);
num_images = numel(info);
Template=zeros(info(1).Height,info(1).Width,length(info),'uint8');
for k = 1:num_images
    image_temp = imread(fname, k, 'Info', info)';image_temp=double(image_temp);
    image_temp=image_temp/max(max(image_temp));image_temp=image_temp*128;
    Template(:,:,k) = image_temp';    
end
clearvars info fname num_images k i j 

PrismLaterality={};
Midline=306;imageSizeY = size(Template,1);
Index_ELO=strfind(Fish_list, 'ELO');
ELO_fish=find(not(cellfun('isempty', Index_ELO)));
ERO_fish=find(cellfun('isempty', Index_ELO));counter=1;
for i=1:length(ItiaList)    
    regionName=ItiaList{i};
    idx_temp=find(LinReg.(regionName).rsquared>0.1);
    ROIs_temp=ROIsPerBrain.(regionName).ROIs(idx_temp,:);
    Numbers=ROIsPerBrain.(regionName).Numbers;Numbers=[0 Numbers];
    counter=1;
    Prism_temp=nan(2,length(LinReg.(regionName).GoodBetas)*3);
    for j=1:length(LinReg.(regionName).GoodBetas)
        idx_goodROIs=find(LinReg.(regionName).KmeansIdx_select==LinReg.(regionName).GoodBetas(j));
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
        if length(idx_goodROIs)==(length(idx_ELO)+length(idx_ERO))
            goodROIs=ROIs_temp(idx_goodROIs,:);
            goodROIs_ELO_y{i,j}=(goodROIs(idx_ELO,2)-Midline)/(imageSizeY-Midline);
            goodROIs_ERO_y{i,j}=(goodROIs(idx_ERO,2)-Midline)/(imageSizeY-Midline);
            if ~isempty(goodROIs_ELO_y{i,j})
                Prism_temp(1,counter)=nanmean(goodROIs_ELO_y{i,j});
                Prism_temp(1,counter+1)=nanstd(goodROIs_ELO_y{i,j});
                Prism_temp(1,counter+2)=length(ELO_fish);
            end
            if ~isempty(goodROIs_ERO_y{i,j})
                Prism_temp(2,counter)=nanmean(goodROIs_ERO_y{i,j});
                Prism_temp(2,counter+1)=nanstd(goodROIs_ERO_y{i,j});
                Prism_temp(2,counter+2)=length(ERO_fish);
            end
            counter=counter+3;
        else
            break
        end
        PrismLaterality{i}=Prism_temp;
    end
end
clearvars goodROIs ROIs_temp regionName idx_temp Template2 idx_slice roi_nb xcoord ycoord circlePixels image_temp columnsInImage rowsInImage col corr_temp Fighdandle framerate GoodBet_temp H idx_g Model_cereb
