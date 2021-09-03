%Get all the files in the folder
MatFiles=dir('*analysis_matlab.mat');
MatFiles_names={MatFiles.name};

%Design a regular expression that match your naming scheme (https://regexr.com/)
fin = cellfun(@(x)regexp(x,'fish(\d+)_ENS(.*)_(\d)DPF','tokens'), MatFiles_names, 'UniformOutput', false);
names=zeros(3,length(fin));
for i=1:length(fin)
    %Fish=str2num(strcat(fin{i}{1}{1},fin{i}{1}{2}));%Concatenate all the matched groups from the regex
    names(1,i)=str2num(fin{i}{1}{1});names(2,i)=str2num(fin{i}{1}{3});
    if strcmp(fin{i}{1}{2},'GF')
        names(3,i)=2;
    elseif strcmp(fin{i}{1}{2},'Fed')
        names(3,i)=3;
    else
        names(3,i)=1;
    end      
end
FishList=unique(names);
clearvars i Fish fin

% MatFiles_fish=[];
% for IndividualFish=FishList
%     IndividualFish=num2str(IndividualFish);
%     IndexC=strfind({MatFiles.name}, strcat(IndividualFish(1),'_',IndividualFish(2))); %Make the string match the pattern of the naming scheme you used
%     MatFiles_fish = [MatFiles_fish find(not(cellfun('isempty', IndexC)))];
% end
% MatFiles_ordered=MatFiles(MatFiles_fish); %This orders the Matlab files per fish, it helps with the indexing of the ROIs after ANTs
[~,idx_sort]=sortrows(names',[2 3]);
MatFiles_ordered=MatFiles(idx_sort'); %This orders the Matlab files per fish, it helps with the indexing of the ROIs after ANTs
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
ROI_temp=reshape(full(Rs),dims(1),dims(2),size(Rs,2));
ROI_GPU=gpuArray(ROI_temp);
ROI_centroid=zeros(size(ROI_temp,3),2);
for i=1:size(ROI_GPU,3)
    temp=squeeze(ROI_GPU(:,:,i));   
    s=regionprops(temp>0,temp,'WeightedCentroid');
    ROI_centroid(i,:)=s.WeightedCentroid;    
end
MatFiles_ordered(1).ROIs=ROI_centroid;


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
    ROI_temp=reshape(full(Rs),dims(1),dims(2),size(Rs,2));
    ROI_GPU=ROI_temp;%gpuArray(ROI_temp);
    ROI_centroid=zeros(size(ROI_temp,3),2);
    for ij=1:size(ROI_GPU,3)
        temp=squeeze(ROI_GPU(:,:,ij));
        s=regionprops(temp>0,temp,'WeightedCentroid');
        ROI_centroid(ij,:)=s.WeightedCentroid;
    end
    GC=C(F,:);
    GN=N(F,:);
    if size(GC,2)<1200
        GC=padarray(GC,[0 1200-size(GC,2)],'replicate','post');
        GN=padarray(GN,[0 1200-size(GN,2)],'replicate','post');
    end
    GoodComponents=horzcat(GoodComponents,F);
    GoodCalcium=vertcat(GoodCalcium,GC);
    GoodNoise=vertcat(GoodNoise,GN);
    MatFiles_ordered(i).GoodNumber=MatFiles_ordered(i-1).GoodNumber+length(F);
    MatFiles_ordered(i).ROIs=ROI_centroid;
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
idx_DPF=nan(size(ZS,2),1);
idx_Fish=nan(size(ZS,2),1);
idx_Treatment=nan(size(ZS,2),1);
names_ordered=names(:,idx_sort);
for i=1:length(MatFiles_ordered)	    
    idx_DPF(Numbers(i)+1:Numbers(i+1))=names(2,i);
    idx_Fish(Numbers(i)+1:Numbers(i+1))=names(1,i);
    idx_Treatment(Numbers(i)+1:Numbers(i+1))=names(3,i);
end
clearvars i Fish Plane name counter

figure;
for i=1:length(MatFiles_ordered)
    im_temp=imread(strrep(MatFiles_ordered(i).name,'output_analysis_matlab.mat','mean.tiff'));
    imagesc(im_temp);colormap gray;
    title(num2str(i));
    hold on;
    scatter(MatFiles_ordered(i).ROIs(:,1),MatFiles_ordered(i).ROIs(:,2),40,'.','r');
    hold off;
    pause
end

fid = fopen('Gut limits');
Limits=inf(2,length(MatFiles_ordered));
counter=1;
while ~feof(fid) % feof(fid) is true when the file ends
      textLineEntry = fgetl(fid); % read one line
      if regexp(textLineEntry,'x')
        Limits(1,counter)=str2num(textLineEntry(3:end));
      elseif regexp(textLineEntry,'y')
        Limits(2,counter)=str2num(textLineEntry(3:end));
      else
          warning('error');
      end
      counter=counter+1;
end
fclose(fid); % close the file


%You can also build an index of the Fish and the plane of each ROI
Numbers=[0 [MatFiles_ordered.GoodNumber]];
idx_Gut=nan(size(ZS,2),1);
for i=1:length(MatFiles_ordered)	    
    ROI_temp=MatFiles_ordered(i).ROIs;
    idx_ROI=(ROI_temp(:,1)<Limits(1,i) & ROI_temp(:,1)>20) & (ROI_temp(:,2)<Limits(2,i) & ROI_temp(:,2)>20);
    idx_Gut(Numbers(i)+1:Numbers(i+1))=idx_ROI;    
end
clearvars i Fish Plane name counter
idx_Gut=boolean(idx_Gut);

Fighandle=figure;
set(Fighandle, 'Position', [10,10, 1200, 800]);
xplot=5;
yplot=3;
ha = tight_subplot(yplot,xplot,[.01 .01],[.01 .01],[.01 .01]);
counter=1;
for treat=1:3
    for dpf=3:7    
        axes(ha(counter));
        idx_temp=find(idx_DPF==dpf & idx_Treatment==treat);
        idx_temp=idx_temp(idx_Gut(idx_temp));
        imagesc(ZS(idx_temp(randperm(length(idx_temp))),:),[-0.5 4]);colormap(hot);
        counter=counter+1;
        set(gca,'XTickLabel',[]);set(gca,'YTickLabel',[]);
    end    
end
print(Fighandle,'ZS_DPF_Treatment','-dsvg','-r0');    

ZS_guts={};
for treat=1:3
    for dpf=3:7  
        idx_temp=find(idx_DPF==dpf & idx_Treatment==treat);
        idx_temp=idx_temp(idx_Gut(idx_temp));
        ZS_guts{treat,dpf-2,1}=ZS(idx_temp,:); 
        ZS_guts{treat,dpf-2,2}=idx_temp;
    end    
end

ops.nCall=[3 100];
ops.useGPU=0;
smooth_ZS_fish={};
results_ENS=struct();
for treat=1:3
    for dpf=3:7
        idx_temp=ZS_guts{treat,dpf-2,2};
        for fish_nb=1:max(unique(idx_Fish(idx_temp)))
            ZS_temp=ZS(idx_temp(idx_Fish(idx_temp)==fish_nb),:);
            if ~isempty(ZS_temp)
                ops.iPC=1:size(ZS_temp,1);
                [isort1, isort2, Sm] = mapTmap(ZS_temp,ops);
                smooth_ZS_fish{treat,dpf-2,fish_nb,1}=isort1;
                smooth_ZS_fish{treat,dpf-2,fish_nb,2}=isort2;
                smooth_ZS_fish{treat,dpf-2,fish_nb,3}=Sm;
                pks_smooth=cell(1,size(Sm,1));
                locs_smooth=cell(1,size(Sm,1));
                Periodicity_smooth=nan(size(Sm,1),1);
                for roi_nb=1:size(ZS_temp,1)
                    [pks_smooth{roi_nb},locs_smooth{roi_nb}] = findpeaks(Sm(roi_nb,:),'MinPeakProminence',1.5,'MinPeakDistance',5);
                    if length(locs_smooth{roi_nb})>5
                        locs_temp=locs_smooth{roi_nb}(2:end)-locs_smooth{roi_nb}(1:end-1);
                        locs_temp=locs_temp/2;
                        Periodicity_smooth(roi_nb)=mean(locs_temp);
                    end
                end
                results_ENS(treat,dpf-2,fish_nb).peaks=pks_smooth;
                results_ENS(treat,dpf-2,fish_nb).locs=locs_smooth;
                results_ENS(treat,dpf-2,fish_nb).Periodicity=Periodicity_smooth;
            end
        end
    end
end

edge=[0:5:120];
treat=1;
h=nan(length(edge)-1,5*8);
for dpf=3:7
    for fish_nb=1:8
        h(:,fish_nb+(dpf-3)*8)=histcounts(results_ENS(treat,dpf-2,fish_nb).Periodicity,edge,'Normalization','probability')';
    end
end


PrismTemp=nan(3,5*8);
for treat=1:3
    for dpf=3:7
        for fish_nb=1:8
            PrismTemp(treat,fish_nb+(dpf-3)*8)=nanmean(results_ENS(treat,dpf-2,fish_nb).Periodicity);
        end
    end
end
PrismTemp(PrismTemp==0)=nan;

edge=[0:5:120];
treat=3;
h=nan(length(edge)-1,5*8);
for dpf=3:7
    for fish_nb=1:8
        h(:,fish_nb+(dpf-3)*8)=histcounts(results_ENS(treat,dpf-2,fish_nb).Periodicity,edge,'Normalization','probability')';
    end
end


%% Testing cross_correlation vs correlation
Cross_corr_temp=nan(size(ZS_temp,1));
Cross_corr_lag=nan(size(ZS_temp,1));
for rows=1:size(ZS_temp,1)
    for columns=1:size(ZS_temp,1)
        if columns>=rows            
            [Cross_corr_temp(rows,columns) Cross_corr_lag(rows,columns)]=max(xcorr(ZS_temp(rows,:)',ZS_temp(columns,:)',10,'coeff'));
        end
    end
end
Cross_corr_temp(isnan(Cross_corr_temp))=0;
Cross_corr_temp(1:1+size(Cross_corr_temp,1):end) = 1;
temp=Cross_corr_temp+tril(Cross_corr_temp',-1);
Cross_corr_temp_vector = squareform(1-temp,'tovector');
corr_temp=pdist(ZS_temp,'correlation');
CrossCorrLink=linkage(Cross_corr_temp_vector,'single');
CorrLink=linkage(corr_temp,'single');

figure;
subplot(1,2,1);dendrogram(CrossCorrLink,0);
subplot(1,2,2);dendrogram(CorrLink,0);

CrossCorrClust=cluster(CrossCorrLink,'cutoff',1);
CorrClust=cluster(CorrLink,'cutoff',1);

CrossCorrClust=cluster(CrossCorrLink,'maxclust',10);
CorrClust=cluster(CorrLink,'maxclust',10);


CrossClustMean=zeros(max(CrossCorrClust),size(ZS_temp,2));
for i=1:max(CrossCorrClust)
    CrossClustMean(i,:)=mean(ZS_temp(CrossCorrClust==i,:),1);
end


ClustMean=zeros(max(CrossCorrClust),size(ZS_temp,2));
for i=1:max(CrossCorrClust)
    ClustMean(i,:)=mean(ZS_temp(CorrClust==i,:),1);
end

figure;
subplot(1,2,1);imagesc(CrossClustMean,[0 1]);colormap magma
subplot(1,2,2);imagesc(ClustMean,[0 1]);colormap magma

figure;
subplot(1,2,1);histogram(CrossCorrClust)
subplot(1,2,2);histogram(CorrClust)