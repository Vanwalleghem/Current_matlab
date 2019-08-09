%Get all the files in the folder
MatFiles=dir('*analysis_matlab.mat');
MatFiles_names={MatFiles.name};

%Design a regular expression that match your naming scheme (https://regexr.com/)
fin = cellfun(@(x)regexp(x,'fish(\d+)_(\d)','tokens'), MatFiles_names, 'UniformOutput', false);
names=[];
for i=1:length(fin)
    Fish=str2num(strcat(fin{i}{1}{1},fin{i}{1}{2}));%Concatenate all the matched groups from the regex
    names(i)=Fish;
end
FishList=unique(names);
clearvars i Fish names fin

MatFiles_fish=[];
for IndividualFish=FishList
    IndividualFish=num2str(IndividualFish);
    IndexC=strfind({MatFiles.name}, strcat(IndividualFish(1),'_',IndividualFish(2))); %Make the string match the pattern of the naming scheme you used
    MatFiles_fish = [MatFiles_fish find(not(cellfun('isempty', IndexC)))];
end
MatFiles_ordered=MatFiles(MatFiles_fish); %This orders the Matlab files per fish, it helps with the indexing of the ROIs after ANTs
MatFiles_ordered = rmfield(MatFiles_ordered,'isdir');MatFiles_ordered = rmfield(MatFiles_ordered,'folder');MatFiles_ordered = rmfield(MatFiles_ordered,'bytes');MatFiles_ordered = rmfield(MatFiles_ordered,'date');
MatFiles_ordered = rmfield(MatFiles_ordered,'datenum');
clearvars MatFiles_fish IndexC IndividualFish

%Load up the data you want from the CaImAn output (DenoisedTraces, Noise, Baseline, ROIs, Spikes and idx_components(the "good" components)) 
name=MatFiles_ordered(1).name;
DenoisedCalcium=load(name, 'DenoisedTraces');
DenoisedCalcium=DenoisedCalcium.DenoisedTraces;
Noise=load(name, 'Noise'); %The noise may contain the weak responses and/or the inhibition since CaImAn convolves with an exponential decay to clean up
Noise=Noise.Noise;
GoodComponents=load(name, 'idx_components');
GoodComponents=GoodComponents.idx_components+1;
GoodCalcium=DenoisedCalcium(GoodComponents,:);
GoodNoise=Noise(GoodComponents,:);
MatFiles_ordered(1).GoodNumber=length(GoodComponents);

%Loads up the ROIs and compute their centroid
Rs=load(name, 'ROIs');
Rs=Rs.ROIs;Rs=Rs(:,GoodComponents);
cor_name=strrep(name,'analysis_matlab','correlation');
cor_im=load(cor_name);cor_im=cor_im.Correlation_image;
dims=size(cor_im);
ROI=reshape(full(Rs),dims(1),dims(2),size(Rs,2));
Centroids=zeros(size(Rs,2),2);
for roi_nb=1:size(ROI,3)
    progressbar([],roi_nb/size(ROI,3));
    temp=regionprops(uint16(squeeze(ROI(:,:,roi_nb)))==max(max(uint16(squeeze(ROI(:,:,roi_nb))))),'Centroid');    
    temp=temp.Centroid;
    Centroids(roi_nb,1:2)=temp;
end
MatFiles_ordered(1).ROIs=Centroids;

%Now that the initial variables are created you iterate over all the files
clearvars Noise Calcium
for i = 2:length(MatFiles_ordered)
    progressbar(i/length(MatFiles_ordered),[]);
    name=MatFiles_ordered(i).name;
    C=load(name, 'DenoisedTraces');
    C=C.DenoisedTraces;
    N=load(name, 'Noise');
    N=N.Noise;
    Rs=load(name, 'ROIs');
    Rs=Rs.ROIs;
    F=load(name, 'idx_components');
    F=F.idx_components+1;
    Rs=Rs(:,F);
    cor_name=strrep(name,'analysis_matlab','correlation');
    cor_im=load(cor_name);cor_im=cor_im.Correlation_image;
    dims=size(cor_im);
    ROI=reshape(full(Rs),dims(1),dims(2),size(Rs,2));
    Centroids=zeros(size(Rs,2),2);
    for roi_nb=1:size(ROI,3)
        progressbar([],roi_nb/size(ROI,3));
        temp=regionprops(uint16(squeeze(ROI(:,:,roi_nb)))==max(max(uint16(squeeze(ROI(:,:,roi_nb))))),'Centroid');
        temp=temp.Centroid;
        Centroids(roi_nb,1:2)=temp;
    end
    GC=C(F,:);
    GN=N(F,:);
    GoodComponents=horzcat(GoodComponents,F);
    GoodCalcium=vertcat(GoodCalcium,GC);
    GoodNoise=vertcat(GoodNoise,GN);
    MatFiles_ordered(i).GoodNumber=MatFiles_ordered(i-1).GoodNumber+length(F);
    MatFiles_ordered(i).ROIs=Centroids;
end
%Clear up all the temporary variables so it doesn't clutter the workspace
clearvars temp GC C S F N name i GS Fitness GN Rs GoodComponents cor_im cor_name MatFiles_names MatFiles ROI roi_nb Centroids DenoisedCalcium dims

%We Z-score the data, you can choose to include the noise or not
%Your data should be in a matrix of NbNeurons x TimePoints
ZS=zscore(GoodCalcium+GoodNoise,1,2);

%You may want to save the "raw" data at this stage and remove the non
%z-scored data
save('Raw_data','-v7.3');
clearvars GoodCalcium GoodNoise

%You can also build an index of the Fish and the plane of each ROI
Numbers=[0 [MatFiles_ordered.GoodNumber]];
idx_Plane=nan(size(ZS,2),1);
idx_Fish=nan(size(ZS,2),1);
name=strcat(MatFiles_ordered(1).name);
for i=1:length(MatFiles_ordered)	
    name=strcat(MatFiles_ordered(i).name);
    [Plane,~]=regexp(name,'Slice(\d+)_','tokens','match');Plane=str2num(Plane{1}{1}); %Again double check the regex
    idx_Plane(Numbers(i)+1:Numbers(i+1))=Plane;    
    [Fish,~]=regexp(name,'fish(\d)_(\d)','tokens');Fish=str2num(strcat(Fish{1}{1},Fish{1}{2})); %Again double check the regex
    idx_Fish(Numbers(i)+1:Numbers(i+1))=Fish;    
end
clearvars i Fish Plane name counter

%If you're working with evoked activity, now is the time to create the
%regressor. You need :
%framerate, the GCaMP spike is the one I created for GCaMP6s at a framerate
%of 4 (may need to interpolate it or make a new one.
spike=[0,1.69644104899772,5.13796058542217,8.27886020152244,10.3756715204800,11.8173714529814,12.2425184714093,10.8571417354877,8.80831829681196,6.91339112244670,5.46959264663869,4.30868766622567,3.42533619066766,2.75378443486879,2.18017250852183,1.72816235135824,1.32732537295463,1.00684435500268,0.730210038304555,0.530242444093118,0.362253250339685,0.227668255288566,0.0869242416152502,0.000718266708050853,-0.0828334873368325]';
spike=spike/max(spike);
framerate=4;
%How many stimuli do you have ? Here I use an example from Lena
Stimuli=zeros(5,size(ZS,2));
OMR=[53	151	186];
for i=OMR
    idx=i*framerate;idx=idx;
    Stimuli(1,idx-1:idx-1+length(spike)-1)=spike';
end
for i = [73	108	213];%Loom	303	331	383	439	526	578
    idx=(i+4)*framerate;idx=idx;
    Stimuli(2,idx-1:idx-1+length(spike)-1)=spike';
end
Loom_audio=[308 444 583];
for i=Loom_audio
    idx=i*framerate;idx=idx;
    Stimuli(3,idx-1:idx-1+length(spike)-1)=spike';
end
for i = [99	134	169]%Audio
    idx=i*framerate;idx=idx;
    Stimuli(4,idx-1:idx-1+length(spike)-1)=spike';
end
for i = [92	127	204]%Vestib
    idx=i*framerate;idx=idx;
    Stimuli(5,idx-1:idx-1+length(spike)-1)=spike';
end
%You may choose to include a linear increase to remove decay or what not
temp=[0:1:size(ZS,2)-1];
Stimuli=Stimuli(:,1:size(ZS,2));
Stimuli(13,:)=temp;
clearvars temp OMR OMR_Audio Audio_vestib Fade_out_Loom i j k Loom_vestib Loom_audio Fade_out Fade_in OMR

%A simpler example from Itia
Stimuli=zeros(3,size(ZS,2));
start=40;
for i = 1:3
    for j = 0:2
        Stimuli(i,start+j*120:start+j*120+length(spike)-1)=spike';
    end
    start=start+40;
end

%If you have the fish with concatenated data run this
load('Mite_InBrainAndTiming.mat');
Stimuli=RebRegressorAudHab;

%Your stimulus should look like a multicolored set of GCaMP waves
figure;plot(Stimuli');

%Our basic analysis starts with a Linear Regression
%You can choose to only use part of the timeseries
LinearModelStimuli=[];
parfor i=1:size(ZS,2)    
    mdl=fitlm(Stimuli',ZS(i,:));
    LinearModelStimuli(i).coef=mdl.Coefficients;    
    LinearModelStimuli(i).rsquared=mdl.Rsquared.Adjusted;
end

%You can pick a Threshold to select the high fidelity responses
ThresholdR2=0.1
idx_rsq=find([LinearModelStimuli.rsquared]>ThresholdR2);
ZS_rsq=ZS(idx_rsq,:);

%You can also get a matrix of all the regression coefficients
coefficients={}; %%%to make the coefficients variable that we will use. Regression coefficients represent the mean change in the response variable for one unit of change in the predictor variable while holding other predictors in the model constant.
for idx=1:length(LinearModelStimuli)%%% to make a variable the size of ModelResults
    coef=[LinearModelStimuli(idx).coef];%%%% and then put in another variable the coef field from ModelResults
    temp=coef.Properties.RowNames;temp=regexp(temp,'x(\d+)','tokens');%%%to take the name of the rows of the coef variable
    if ~isempty(temp)%%% if temp is not empty...
        %temp=[temp{:}];temp=[temp{:}];temp=[temp{:}];%temp=str2num(temp);
        for coef_idx=2:height(coef)%%%take the number of rows from coef, except the first one(i think because is the intercept)
            %if coef.pValue(coef_idx)<0.05%%%to select the coef that are bellow the p value we want, in this case 0.05
                coefficients{idx,str2num(temp{coef_idx}{1}{1})}=coef.Estimate(coef_idx); %%%to make an array the size of idx,10 with the coefficient values that were significant
            %end
        end
    end
end
idxempty=cellfun('isempty',coefficients); %%%to make a variable with where we will aply in every cell the isempty function wich will help us find the empty places
coefficients(idxempty)={0}; %%% and put a 0 in the places where we found that there were empty cells
clearvars idxempty idx coef_idx coef  %%%clear variables
coefficients=cell2mat(coefficients); %%%to make a matrix of the coefficients array

%The alternative/follow up on linear regression is clustering, we usually
%use k-means, before or after filtering, as you prefer
options = statset('UseParallel',1); %parallelize the replicates
[idxKmeans_ZS Cmap_ZS]=kmeans(ZS,20,'Options',options,'Distance','cityblock','Replicates',3,'MaxIter',1000,'Display','final');
figure;
imagesc(Cmap_ZS);

%Build the ROI csv files for ANTs (not needed here)
progressbar;
for fish_nb=1:length(FishList)
    temp=num2str(FishList(fish_nb));
    IndexC=strfind({MatFiles_ordered.name}, strcat(temp(1),'_',temp(2))); %Again double check the pattern
    MatFiles_fish = find(not(cellfun('isempty', IndexC)));
    progressbar(fish_nb/length(FishList));
    for file_nb=1:length(MatFiles_fish)
        progressbar([],file_nb/length(MatFiles_fish));
        [Plane,~]=regexp(MatFiles_ordered(MatFiles_fish(file_nb)).name,'Slice(\d+)_','tokens','match');Plane=str2num(Plane{1}{1});
        if file_nb==1
            Centroids=MatFiles_ordered(MatFiles_fish(file_nb)).ROIs;
            Centroids(:,3)=Plane;
        else
            CentTemp=MatFiles_ordered(MatFiles_fish(file_nb)).ROIs;CentTemp(:,3)=Plane;
            Centroids=vertcat(Centroids,CentTemp);
        end
    end
    if length(Centroids)~= sum(idx_Fish==FishList(fish_nb));
        temp
        break
    end
    csvwrite(strcat('_ROIs',temp,'.csv'),Centroids);
end
clearvars slice fish_nb roi_nb ROI Centroids IndexC file_nb

%Once ANTs is processed you should load back the csv files
CSV_Files=dir('_2Warped*.csv');
MatFiles_fish=[];
for IndividualFish=FishList
    IndividualFish=num2str(IndividualFish);
    IndexC=strfind({CSV_Files.name}, IndividualFish);    
    MatFiles_fish = [MatFiles_fish find(not(cellfun('isempty', IndexC)))];
end
CSV=CSV_Files(MatFiles_fish);%To order them in the same way as the Z-scored data

ROIs_warped=struct();truth=[];
for i=1:length(CSV)
    temp=csvread(CSV(i).name,1);   
    Fishname=regexp(CSV(i).name,'_2Warped(\d+).csv','tokens');Fishname=Fishname{1}{1};    %Regex yadda yadda
    ROIs_warped(i).name=Fishname;    
    ROIs_warped(i).coord=temp(:,1:3);    
    truth(i)=size(temp,1)==sum(idx_Fish==str2num(Fishname)) %checks the numbers of ROI match the expected number
end
clearvars i temp Fishname

%Sorting the data and checks everything is fine
temp_nb=0;truth=[];
MatFiles_names={MatFiles_ordered.name};
Sorted_ROIs=nan(length(idx_Fish),3);
for fish_nb=1:length(FishList)
    temp=num2str(ROIs_warped(fish_nb).name);
	IndexC=strfind({MatFiles_ordered.name}, strcat(temp(1:3),'_fish0',temp(end)));    
    MatFiles_fish = find(not(cellfun('isempty', IndexC)));
    ROI=ROIs_warped(fish_nb).coord;
    Sort_ROIs=[];Counter_ROI_coord=1;
    for file_nb=1:length(MatFiles_fish)
        if file_nb>1
            Counter_ROI_coord=Counter_ROI_coord+length(Sorted_ROIs(numbersForROIs(1):numbersForROIs(2)));
        end
        if MatFiles_fish(file_nb)==1
            numbersForROIs=[1 MatFiles_ordered(1).GoodNumber];
        else
            numbersForROIs=[MatFiles_ordered(MatFiles_fish(file_nb)-1).GoodNumber+1 MatFiles_ordered(MatFiles_fish(file_nb)).GoodNumber];
        end
        if ismember(numbersForROIs,Sort_ROIs)
            temp
            break
        end        
        Sorted_ROIs(numbersForROIs(1):numbersForROIs(2),:)=ROI(Counter_ROI_coord:Counter_ROI_coord+length(Sorted_ROIs(numbersForROIs(1):numbersForROIs(2)))-1,:);
    end    
    if ~length(Sort_ROIs)-temp_nb==sum(idx_Fish==(FishList(fish_nb)))
        temp
        break
    end
    temp_nb=length(Sort_ROIs);
end
clearvars slice fish_nb roi_nb ROI Centroids IndexC file_nb
sum(~isnan(Sorted_ROIs(:,1))) %How many errors there are

%Now you need to round up the ROI coordinates to match with Z-brain and
%ensure that the orientation matches the Z-brain index
ROI_fish=round(Sorted_ROIs);
ROI_correct(:,1)=round(ROI_fish(:,2));
ROI_correct(:,2)=round(ROI_fish(:,1));
ROI_correct(:,3)=round(ROI_fish(:,3)/2);

load('Zbrain_Masks.mat')
Zbrain_AllMask=vertcat(Zbrain_Masks{[1:1:77 79:1:294],3}); %Makes a massive matrix of all the Zbrain regions except the eyes
Zbrain_AllMask=unique(Zbrain_AllMask,'rows');

%You can check the alignments with :
idx_rand=randsample(length(Zbrain_AllMask),10000); %need to select a few members only
figure;
subplot(1,2,1);
scatter(Zbrain_AllMask(idx_rand,1),Zbrain_AllMask(idx_rand,2),'.')
hold on;scatter(ROI_correct(:,1),ROI_correct(:,2),'.');
subplot(1,2,1);
scatter(Zbrain_AllMask(idx_rand,1),Zbrain_AllMask(idx_rand,3),'.')
hold on;scatter(ROI_correct(:,1),ROI_correct(:,3),'.');

%Then we subdivide the responses for each subregion selected in RegionList
PerBrainRegions=struct();
RegionList={'Thalamus','Cerebellum','Semicircularis','Telencephalon','Tectum','Tegmentum','Habenula','Pretectum','MON','Hindbrain','Stratum'};
progressbar;
for i=1:length(RegionList)
    progressbar(i/length(RegionList));
    regionName=RegionList{i};
    if strcmp(regionName,'Telencephalon')
        Mask=Zbrain_Masks{294,3};
    elseif strcmp(regionName,'Hindbrain')
        Hindbrain_Mask=Zbrain_Masks{259,3};
        Mask=Zbrain_Masks{131,3};
        IsInEyes_temp=ismember(Hindbrain_Mask,Mask,'rows');IsInEyes_temp=find(IsInEyes_temp==1);%remove cerebellum
        Hindbrain_Mask(IsInEyes_temp,:)=[];
        Mask=Zbrain_Masks{295,3};
        IsInEyes_temp=ismember(Hindbrain_Mask,Mask,'rows');IsInEyes_temp=find(IsInEyes_temp==1);%remove MON
        Hindbrain_Mask(IsInEyes_temp,:)=[];
        Mask=Hindbrain_Mask;
    else
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
    end
    Mask=unique(Mask,'rows');
    IsInBrainRegion=ismember(ROI_correct,Mask,'rows');
    PerBrainRegions.(regionName).idx=find(IsInBrainRegion==1);    
end
clearvars i j k l m temp Mask Hindbrain_Mask IsInBrainRegion IsInEyes_temp IndexC

%And now you decide what to do with it, like histograms of the brain
%region contribution to something
progressbar;
h=zeros(length(RegionList),length(FishList));
for i=1:length(RegionList)
    progressbar(i/length(RegionList));
    regionName=RegionList{i};
    idx_temp=PerBrainRegions.(regionName).idx;    
    for fish=1:length(FishList)
        temp_WT_fish=length(find(idx_Fish(idx_temp)==FishList(fish)));
        h(i,fish)=temp_WT_fish;
    end    
end
figure;bar(sum(h,2));

%Or K-means per brain region or whatever you like