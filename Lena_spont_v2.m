%Need to make the stacks for 1112_fish6 and 1122_fish2
MatFiles=dir('*analysis_matlab.mat'); %%to get the files

MatFiles_names={MatFiles.name};
fin = cellfun(@(x)regexp(x,'2018(\d+)_fish(\d+)','tokens'), MatFiles_names, 'UniformOutput', false);
names=[];
for i=1:length(fin)
    Fish=str2num(strcat(fin{i}{1}{1},fin{i}{1}{2}));
    names(i)=Fish;
end
WT_fishNB=unique(names);

MatFiles_fish=[];
for WT_fish=WT_fishNB
    WT_fish=num2str(WT_fish);
    IndexC=strfind({MatFiles.name}, strcat(WT_fish(1:4),'_fish',WT_fish(end)));    
    MatFiles_fish = [MatFiles_fish find(not(cellfun('isempty', IndexC)))];
end
MatFiles=MatFiles(MatFiles_fish);

name=strcat(MatFiles(1).name); %%%to get the name of the files
Calcium=load(name, 'DenoisedTraces'); %%to load only the DenoisedTraces from the file, the raw data was denoised by the CNMF (The Cluster Analysis tool calculates clusters based on a Constrained non-negative matrix factorization (NMF) clustering method.)
Calcium=Calcium.DenoisedTraces; %%%% <-- take the field called DenoisedTraces from the Calcium structure and make it the new Calcium
%MatFiles(1).number=size(Calcium,1);
%Spikes=load(name, 'Spikes');
%Spikes=Spikes.Spikes;
Noise=load(name, 'Noise');
Noise=Noise.Noise;
%DF=load(name, 'dFonF');
%DF=DF.dFonF;
Fitness=load(name, 'idx_components');%%to load only the idx_components from the file, they are based on what a Gcamp spike should be and they will filter the true spikes in our data
Fitness=Fitness.idx_components+1; %%%% <-- take the field called idx_components from the Fitness structure and make it the new Fitness but why +1?? Because python indexing starts at 0 ant matlab at 1
GoodCalcium=Calcium(Fitness,:);  %%%to combine the Calcium and Fitness variables (need to ask Gilles what Fitness is). Fitness here is the variable were we take the good calcium responses from the HPC analysis and pairthem with their index number.
%GoodSpikes=Spikes(Fitness,:);
GoodNoise=Noise(Fitness,:);
%GoodDF=DF(Fitness,:);
Rs=load(name, 'ROIs');
Rs=Rs.ROIs;Rs=Rs(:,Fitness);
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
MatFiles(1).ROIs=Centroids;
MatFiles(1).GoodNumber=length(Fitness); %%%% <-- Create a field inside MatFilesCalcium called GoodNumber the size of Fitness.
for i = 2:length(MatFiles) %%%%to take the slices one by one starting by the second one cause we already did this with the first one
    %%%% we are going to do the same thing that before but for all the
    %%%% slices
    progressbar(i/length(MatFiles),[]);
    name=strcat(MatFiles(i).name);%%%%to take the name of the slice in turn
    C=load(name, 'DenoisedTraces');%%to load only the DenoisedTraces from the file
    C=C.DenoisedTraces;%%%% <-- take the field called DenoisedTraces from the C structure and make it the new C
    N=load(name, 'Noise');
    N=N.Noise;
    F=load(name, 'idx_components');
    F=F.idx_components+1;%%%because indexing in python is from 0 and matlab is at 1
    GC=C(F,:);
    Noise=vertcat(Noise,N);
    GN=N(F,:);
    Fitness=horzcat(Fitness,F);
    GoodCalcium=vertcat(GoodCalcium,GC); %The fish 20+ are longer
    GoodNoise=vertcat(GoodNoise,GN);
    Rs=load(name, 'ROIs');
    Rs=Rs.ROIs;Rs=Rs(:,F);
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
    MatFiles(i).ROIs=Centroids;
    MatFiles(i).GoodNumber=MatFiles(i-1).GoodNumber+length(F);
    MatFiles(i).GC=GC;
end
clearvars GC C S F N name i GS GN N;%%%to get rid of vairables we will not use anymore
%%

save('Lena_Spont.mat','GoodCalcium', 'GoodNoise','MatFiles','-v7.3'); clearvars GoodCalcium GoodNoise GoodCalNoise Fitness Calcium Noise

Numbers=[0 [MatFiles.GoodNumber]];
idx_Fish=nan(length(GoodCalcium),1);
name=strcat(MatFiles(1).name);%%%to get the name of the files (is actually to create the variable name before the loop)
for i=1:length(MatFiles) %%%%to take slices one by one
    name=strcat(MatFiles(i).name);
    [name2,~]=regexp(name,'10minsNoStim_(\d+)_','tokens','match'); %%%to get the number of the fish
    [name3,~]=regexp(name,'fish(\d+)_','tokens','match'); %%%to get the number of the fish
    Fish=strcat(name2{1}{1},name3{1}{1});Fish=str2double(Fish); %%%to get the number of the fish
    idx_Fish(Numbers(i)+1:Numbers(i+1))=Fish; %%%to put the number of the plane on the correspondent goodnumbers in the idx_Fish
end
clearvars i Fish Plane name counter name2 name3 rows %%%to get rid of vairables we will not use anymore
FishList=unique(idx_Fish);

Corr_matrices=cell(length(FishList),1);
progressbar;
for i=1:length(FishList)
    fish_nb=FishList(i);
    idx_temp=idx_Fish==fish_nb;
    Corr_matrices{i}=pdist(ZS_CN(idx_temp,:),'correlation');
    progressbar(i/length(FishList));
end

save('Lena_Spont_CorrMatricesAll.mat','Corr_matrices','-v7.3'); clearvars Corr_matrices

MatFiles_names={MatFiles.name};
fin = cellfun(@(x)regexp(x,'2018(\d+)_fish(\d+)','tokens'), MatFiles_names, 'UniformOutput', false);
names=[];
for i=1:length(fin)
    Fish=str2num(strcat(fin{i}{1}{1},fin{i}{1}{2}));
    names(i)=Fish;
end
WT_fishNB=unique(names);

counter=1;
ROI_coords=cell(length(WT_fishNB),1);
for WT_fish=WT_fishNB
    WT_fish=num2str(WT_fish);
    IndexC=strfind({MatFiles.name}, strcat(WT_fish(1:4),'_fish',WT_fish(end)));    
    IndexC=find(not(cellfun('isempty', IndexC)));
    ROI_pool=[];
    for slice_nb=1:length(IndexC)
        Rois_temp=MatFiles(IndexC(slice_nb)).ROIs;
        [slice,~]=regexp(MatFiles(IndexC(slice_nb)).name,'Slice(\d+)_','tokens','match');slice=str2num(slice{1}{1});
        Rois_temp(:,3)=slice;
        ROI_pool=[ROI_pool ;Rois_temp];
    end
    ROI_coords{counter}=ROI_pool;
    counter=counter+1;
end

for i=1:length(WT_fishNB)
    WT_fish=num2str(WT_fishNB(i));
    csvwrite(strcat('_ROIs',WT_fish,'.csv'),ROI_coords{i});
end

ROI_coords_warped=cell(length(WT_fishNB),1);
for i=46:length(WT_fishNB)
    WT_fish=num2str(WT_fishNB(i));
    FileName=strcat('_2Warped',WT_fish(1:4),'_fish',WT_fish(5),'.csv');
    ROI_coords_warped{i}=csvread(FileName,1);
end

ROI_fish=cell(length(WT_fishNB),1);
for i=1:length(WT_fishNB)
    if ~isempty(ROI_coords_warped{i})
        ROI_WT=[];
        ROI_temp=round(ROI_coords_warped{i}(:,1:3));        
        ROI_WT(:,1)=round(ROI_temp(:,2));
        ROI_WT(:,2)=round(ROI_temp(:,1));
        ROI_WT(:,3)=round(ROI_temp(:,3)/2);
        ROI_fish{i}=ROI_WT;
    end
end

load('Zbrain_Masks.mat')
Zbrain_AllMask=vertcat(Zbrain_Masks{[1:1:77 79:1:294],3});
Zbrain_AllMask=unique(Zbrain_AllMask,'rows');

for i=1:length(WT_fishNB)
    if ~isempty(ROI_coords_warped{i})
        IsInBrain=ismember(ROI_fish{i,1},Zbrain_AllMask,'rows');%remove cerebellum
        ROI_fish{i,2}=IsInBrain;
    end
end

FishList=WT_fishNB;

%Excluding crappy hets
Genome{1}=[201810222 201810223 201810224 201810232 201810234 201810292 201811054 201811056 201811062 201811127 201811194 201811222];
Genome{2}=[201810225 201810226 201810291 201810293 201810294 201810295 201810296 201810301 201810304 201811052 201811053];
Genome{3}=[201810231 201810233 201810297 201810302 201811055 201811061 201811064 201811065 201811124 201811195];
Genome{4}=[201811063 201811125 201811126 201811191 201811192 201811196 201811223 201811224];

%Exclude blurry fish
% Genome{1}=[201810222 201810223 201810224 201810232 201810234 201810292 201811054 201811062 201811127];%mutant
Genome{1}=[201810223, 201810224, 201810234, 201811054, 201811056, 201811062, 201811127];%mutant
Genome{2}=[201810225 201810291 201810294 201810301 201811053 201811224];
% Genome{3}=[201810233 201810297 201810302 201811055 201811064 201811124];
Genome{3}=[ 201810233 201810297 201811055 201811064 201811124];
Genome{4}=[201810226 201810293 201810295  201811052 201811063 201811125 201811191 201811196];%Hets


for idx_gen=[3]
    idx_temp=find(ismember(FishList,Genome{idx_gen}));        
    Correlation_WT_all=[];
    for fish_nb=1:length(idx_temp)
        if ~isempty(ROI_fish{idx_temp(fish_nb),2})     
            Correlation_WT_all=[Correlation_WT_all Corr_matrices_noisy{idx_temp(fish_nb)}];
        end
    end
end

for idx_gen=[3]
    idx_temp=find(ismember(FishList,Genome{idx_gen}));        
    Correlation_WT_all_denoised=[];
    for fish_nb=1:length(idx_temp)
        if ~isempty(ROI_fish{idx_temp(fish_nb),2})     
            Correlation_WT_all_denoised=[Correlation_WT_all_denoised Corr_matrices_denoised{idx_temp(fish_nb)}];
        end
    end
end

Threshold_corr_WT=mean(Correlation_WT_all)-2*std(Correlation_WT_all);

MeanSTD_corr_WT_denoised(1)=nanmean(Correlation_WT_all_denoised);
MeanSTD_corr_WT_denoised(2)=nanstd(Correlation_WT_all_denoised);
Threshold_corr_WT_denoised=MeanSTD_corr_WT_denoised(1)-2.5*MeanSTD_corr_WT_denoised(2);

Graphs_correlation_WTthr=cell(length(WT_fishNB),1);
for Fish_nb=1:length(WT_fishNB)
    if ~isempty(ROI_fish{Fish_nb})
        Correlation_matrix=squareform(Corr_matrices{Fish_nb});
        Correlation_matrix_binary=Correlation_matrix(ROI_fish{Fish_nb,2},ROI_fish{Fish_nb,2})<Threshold_corr_WT;
        Graphs_correlation_WTthr{Fish_nb}=graph(Correlation_matrix_binary,'omitselfloops');
    end
end


%corr_matrix on denoised data
Corr_matrices_denoised=cell(length(FishList),1);
progressbar;
for Fish_nb=1:length(FishList)
    fish_nb=FishList(Fish_nb);
    Correlation_matrix=squareform(Corr_matrices{Fish_nb});
    Correlation_matrix=Correlation_matrix(ROI_fish{Fish_nb,2},ROI_fish{Fish_nb,2});    
    Corr_matrices_denoised{Fish_nb}=squareform(Correlation_matrix);
    progressbar(Fish_nb/length(FishList));
end

Graphs_correlation=cell(length(WT_fishNB),1);
for Fish_nb=1:length(WT_fishNB)
    if ~isempty(ROI_fish{Fish_nb})
        Correlation_matrix=squareform(Corr_matrices_denoised{Fish_nb});
        Correlation_matrix_binary=Correlation_matrix<0.25;
        Graphs_correlation{Fish_nb}=graph(Correlation_matrix_binary,'omitselfloops');
    end
end

%corr_matrix on denoised data
Dist_matrices=cell(length(FishList),1);
progressbar;
for i=1:length(FishList)
    ROI_temp=ROI_fish{i,1};
    Dist_matrices{i}=pdist(ROI_temp(ROI_fish{i,2},:));
    progressbar(i/length(FishList));
end
Max_distance_ROIs=0;
for i=1:length(FishList)
    Max_temp=max(Dist_matrices{i});
    if Max_distance_ROIs<Max_temp
        Max_distance_ROIs=Max_temp;
    end
end

edges={[0:0.005:2],[0:0.005:2]};
CorrVsDist=cell(2,length(FishList));
for i=1:length(FishList)   
    progressbar(i/length(FishList));
    if ~isempty(Corr_matrices_denoised{i})
        [CorrVsDist{1,i}, CorrVsDist{2,i}] = hist3([Corr_matrices_denoised{i}' Dist_matrices{i}'/(Max_distance_ROIs/2)],'Edges',edges);    
    end
end

rows=4;yplot=3;counter=1;
for idx_gen=[3 1 2]
    Fighandle=figure;
    set(Fighandle, 'Position', [200, 200, 1500, 1500]);
    idx_temp=find(ismember(FishList,Genome{idx_gen}));
    ha = tight_subplot(rows,yplot,[.02 .02],[.02 .02],[.02 .02]);
    h=[];
    for fish_nb=1:length(idx_temp)        
        axes(ha(fish_nb));
        imagesc(CorrVsDist{1,idx_temp(fish_nb)})
    end
    print(Fighandle,strcat('\\shares.rdm.uq.edu.au\scottlab-q0291\Lena\Multi_results\Spontaneous\','WTvsFMR_CorrelationVsDistance_',num2str(idx_gen),'_NoEyes'),'-dpng','-r0');
end
clearvars rows counter i idx_temp idx_temp2 idx_gen Fighandle yplot
idx_gen=3;
idx_temp=find(ismember(FishList,Genome{idx_gen}));
CorrVDist_temp=zeros(size(CorrVsDist{1,1},1),size(CorrVsDist{1,1},1),length(idx_temp));
for fish_nb=1:length(idx_temp)
    temp=CorrVsDist{1,idx_temp(fish_nb)};
    if ~isempty(temp)
        CorrVDist_temp(:,:,fish_nb)=temp/sum(temp(:));    
    end
end
AVG_WT=squeeze(mean(CorrVDist_temp,3));
idx_gen=1;
idx_temp=find(ismember(FishList,Genome{idx_gen}));
CorrVDist_temp=zeros(size(CorrVsDist{1,1},1),size(CorrVsDist{1,1},1),length(idx_temp));
for fish_nb=1:length(idx_temp)
    temp=CorrVsDist{1,idx_temp(fish_nb)};
    if ~isempty(temp)
        CorrVDist_temp(:,:,fish_nb)=temp/sum(temp(:));        
    end
end
AVG_FMR=squeeze(mean(CorrVDist_temp,3));
idx_gen=2;
idx_temp=find(ismember(FishList,Genome{idx_gen}));
idx_gen=4;
idx_temp2=find(ismember(FishList,Genome{idx_gen}));
idx_temp=union(idx_temp,idx_temp2);
CorrVDist_temp=zeros(size(CorrVsDist{1,1},1),size(CorrVsDist{1,1},1),length(idx_temp));
for fish_nb=1:length(idx_temp)
    temp=CorrVsDist{1,idx_temp(fish_nb)};
    if ~isempty(temp)
        CorrVDist_temp(:,:,fish_nb)=temp/sum(temp(:));
    end
end
AVG_Het=squeeze(mean(CorrVDist_temp,3));
% idx_gen=4;
% idx_temp=find(ismember(FishList,Genome{idx_gen}));
% CorrVDist_temp=zeros(size(CorrVsDist{1,1},1),size(CorrVsDist{1,1},1),length(idx_temp));
% for fish_nb=1:length(idx_temp)
%     temp=CorrVsDist{1,idx_temp(fish_nb)};
%     CorrVDist_temp(:,:,fish_nb)=temp/sum(temp(:));
% end
% AVG_Het2=squeeze(mean(CorrVDist_temp,3));

Fighandle=figure;
set(Fighandle, 'Position', [200, 200, 3200, 800]);
ha = tight_subplot(1,3,[.02 .02],[.02 .02],[.02 .02]);
axes(ha(1));imagesc(AVG_WT,[0 8e-5]);
axes(ha(2));imagesc(AVG_FMR,[0 8e-5]);
axes(ha(3));imagesc(AVG_Het,[0 8e-5]);
%axes(ha(4));imagesc(AVG_Het2,[0 8e-5]);
print(Fighandle,strcat('\\shares.rdm.uq.edu.au\scottlab-q0291\Lena\Multi_results\Spontaneous\','WTvsFMRvsHetmerged_CorrelationVsDistance_norm','_NoEyes'),'-dpng','-r0');

ZS=zscore(GoodCalcium,1,2);
%FishList=FishList+201800000;
%Firing rate
Firing_rate=cell(length(FishList),1);
Firing_rate_mean_std=nan(length(FishList),2);
progressbar(0,0);
for i=1:length(FishList)
    fish_nb=FishList(i);
    if ~isempty(ROI_fish{i,2})
        idx_temp=find(idx_Fish==fish_nb);        
        ZS_temp=ZS(idx_temp(ROI_fish{i,2}),:);
        temp=nan(sum(ROI_fish{i,2}),1);
        for roi_nb=1:length(temp)
            [~,loc]=findpeaks(ZS_temp(roi_nb,:),'MinPeakDistance',8,'MinPeakProminence',1);
            temp(roi_nb)=length(loc);
            progressbar(roi_nb/length(temp),[]);
        end
        Firing_rate{i}=temp;
        Firing_rate_mean_std(i,1)=mean(Firing_rate{i});
        Firing_rate_mean_std(i,2)=std(Firing_rate{i});
    end
    progressbar([],i/length(FishList));
end

for idx_gen=[3 1 2]
    idx_temp=find(ismember(FishList,Genome{idx_gen}));
    if idx_gen==2;
        idx_gen=4;
        idx_temp2=find(ismember(FishList,Genome{idx_gen}));
        idx_temp=union(idx_temp,idx_temp2);
        idx_gen=2;
    end
    h=[Firing_rate_mean_std(idx_temp,1) Firing_rate_mean_std(idx_temp,2)];
    csvwrite(strcat('\\shares.rdm.uq.edu.au\scottlab-q0291\Lena\Multi_results\Spontaneous\Spont_Firing_rate_',num2str(idx_gen),'_meanSTD_NoEyes.csv'),h');
end

Fighandle=figure;
set(Fighandle, 'Position', [200, 200, 1500, 750]);
edges=[0:1:100];
color=[1 0 0; 0 1 0; 0 0 1];
for idx_gen=[3 1 2]
    idx_temp=find(ismember(FishList,Genome{idx_gen}));
    if idx_gen==2;
        idx_gen=4;
        idx_temp2=find(ismember(FishList,Genome{idx_gen}));
        idx_temp=union(idx_temp,idx_temp2);
        idx_gen=2;
    end
    h=[];
    for fish_nb=1:length(idx_temp)
        [h(fish_nb,:),~]=histcounts(Firing_rate{idx_temp(fish_nb)},edges);
    end
    bar(nanmean(h,1),'FaceColor','none','EdgeColor',color(idx_gen,:));hold on;
    csvwrite(strcat('\\shares.rdm.uq.edu.au\scottlab-q0291\Lena\Multi_results\Spontaneous\Spont_Firing_rate_',num2str(idx_gen),'_NoEyes.csv'),h');
    
end
clearvars rows counter i idx_temp idx_temp2 idx_gen Fighandle yplot

rows=4;yplot=1;counter=1;edges=[-1:0.02:1];
Fighandle=figure;
set(Fighandle, 'Position', [200, 200, 1200, 800]);
for idx_gen=1:length(Genome)
    idx_temp=find(ismember(FishList,Genome{idx_gen}));
    %subplot(rows,yplot,idx_gen);
    h=[];
    for fish_nb=1:length(idx_temp)        
        histogram((1-Corr_matrices_noisy{idx_temp(fish_nb)}),edges,'normalization','probability');hold on;
        [h(fish_nb,:),~]=histcounts((1-Corr_matrices_noisy{idx_temp(fish_nb)}),edges);   
    end
    csvwrite(strcat('\\shares.rdm.uq.edu.au\scottlab-q0291\Lena\Multi_results\Spontaneous\BinnedCorr_Spontanous_',num2str(idx_gen),'_noEyes.csv'),h');
    print(Fighandle,strcat('\\shares.rdm.uq.edu.au\scottlab-q0291\Lena\Multi_results\Spontaneous\Spontaneous_hist_corrCoef_',num2str(idx_gen),'_noEyes.svg'),'-dsvg','-r0');
    print(Fighandle,strcat('\\shares.rdm.uq.edu.au\scottlab-q0291\Lena\Multi_results\Spontaneous\Spontaneous_hist_corrCoef_',num2str(idx_gen),'_noEyes.eps'),'-depsc','-r0');
    hold off;
end
clearvars rows counter i idx_temp idx_temp2 idx_gen Fighandle yplot

for idx_gen=1:length(Genome)
    idx_temp=find(ismember(FishList,Genome{idx_gen}));    
    h=[];
    for fish_nb=1:length(idx_temp)                
        h(fish_nb,:)=nanmean((1-Corr_matrices_noisy{idx_temp(fish_nb)}));   
    end
    csvwrite(strcat('\\shares.rdm.uq.edu.au\scottlab-q0291\Lena\Multi_results\Spontaneous\MeanCorr_Spontanous_',num2str(idx_gen),'_noEyes.csv'),h');
    hold off;
end
clearvars rows counter i idx_temp idx_temp2 idx_gen Fighandle yplot

rows=4;yplot=1;counter=1;edges=[-1:0.02:1];
Fighandle=figure;
set(Fighandle, 'Position', [200, 200, 1200, 800]);
for idx_gen=1:length(Genome)
    idx_temp=find(ismember(FishList,Genome{idx_gen}));
    %subplot(rows,yplot,idx_gen);
    h=[];
    for fish_nb=1:length(idx_temp)        
        histogram((1-Corr_matrices_denoised{idx_temp(fish_nb)}),edges,'normalization','probability');hold on;
        [h(fish_nb,:),~]=histcounts((1-Corr_matrices_denoised{idx_temp(fish_nb)}),edges);   
    end
    csvwrite(strcat('\\shares.rdm.uq.edu.au\scottlab-q0291\Lena\Multi_results\Spontaneous\BinnedCorr_Spontanous_denoised_',num2str(idx_gen),'_noEyes.csv'),h');
    print(Fighandle,strcat('\\shares.rdm.uq.edu.au\scottlab-q0291\Lena\Multi_results\Spontaneous\Spontaneous_denoised_hist_corrCoef_',num2str(idx_gen),'_noEyes.svg'),'-dsvg','-r0');
    print(Fighandle,strcat('\\shares.rdm.uq.edu.au\scottlab-q0291\Lena\Multi_results\Spontaneous\Spontaneous_denoised_hist_corrCoef_',num2str(idx_gen),'_noEyes.eps'),'-depsc','-r0');
    hold off;
end
clearvars rows counter i idx_temp idx_temp2 idx_gen Fighandle yplot

for idx_gen=1:length(Genome)
    idx_temp=find(ismember(FishList,Genome{idx_gen}));    
    h=[];
    for fish_nb=1:length(idx_temp)                
        h(fish_nb,:)=nanmean((1-Corr_matrices_denoised{idx_temp(fish_nb)}));   
    end
    csvwrite(strcat('\\shares.rdm.uq.edu.au\scottlab-q0291\Lena\Multi_results\Spontaneous\MeanCorr_Spontanous_denoised_',num2str(idx_gen),'_noEyes.csv'),h');
    hold off;
end
clearvars rows counter i idx_temp idx_temp2 idx_gen Fighandle yplot


load('Nodes_N_means_alldatasets.mat', 'Zbrain_brainMask2D');

Graph_info_perGen=struct();
for idx_gen=[3 1 2]
    idx_temp=find(ismember(FishList,Genome{idx_gen}));
    if idx_gen==2;
        idx_gen=4;
        idx_temp2=find(ismember(FishList,Genome{idx_gen}));
        idx_temp=union(idx_temp,idx_temp2);
        idx_gen=2;
    end
    Degrees_fish=cell(length(idx_temp),1);
    Degrees_all=[];
    for fish_nb=1:length(idx_temp)
        if ~isempty(ROI_fish{idx_temp(fish_nb),2})
            Degrees_fish{fish_nb}=degree(Graphs_correlation{idx_temp(fish_nb)});
            Degrees_all=[Degrees_all; Degrees_fish{fish_nb}];
        end
    end
    Graph_info_perGen(idx_gen).degree_ind=Degrees_fish;
    Graph_info_perGen(idx_gen).degree=Degrees_all;
end

Graphs_correlation_WTthr_denoised=cell(length(WT_fishNB),1);
for Fish_nb=1:length(WT_fishNB)
    if ~isempty(ROI_fish{Fish_nb})
        Correlation_matrix=squareform(Corr_matrices_denoised{Fish_nb});
        Correlation_matrix_binary=Correlation_matrix<Threshold_corr_WT_denoised;
        Graphs_correlation_WTthr_denoised{Fish_nb}=graph(Correlation_matrix_binary,'omitselfloops');
    end
end
clearvars Correlation_matrix Correlation_matrix_binary

Graph_info_perGen_denoised=struct();
for idx_gen=[3 1 2]
    idx_temp=find(ismember(FishList,Genome{idx_gen}));
    if idx_gen==2;
        idx_gen=4;
        idx_temp2=find(ismember(FishList,Genome{idx_gen}));
        idx_temp=union(idx_temp,idx_temp2);
        idx_gen=2;
    end
    Degrees_fish=cell(length(idx_temp),1);
    Degrees_all=[];
    for fish_nb=1:length(idx_temp)
        if ~isempty(ROI_fish{idx_temp(fish_nb),2})
            Degrees_fish{fish_nb}=degree(Graphs_correlation_WTthr_denoised{idx_temp(fish_nb)});
            Degrees_all=[Degrees_all; Degrees_fish{fish_nb}];
        end
    end
    Graph_info_perGen_denoised(idx_gen).degree_ind=Degrees_fish;
    Graph_info_perGen_denoised(idx_gen).degree=Degrees_all;
end

edges=[0:1:50];
for idx_gen=[3 1 2]
    h=[];
    for fish_nb=1:length(Graph_info_perGen_denoised(idx_gen).degree_ind)
        [h(fish_nb,:),~]=histcounts(Graph_info_perGen_denoised(idx_gen).degree_ind{fish_nb},edges);%,'normalization','probability');
    end
    csvwrite(strcat('\\shares.rdm.uq.edu.au\scottlab-q0291\Lena\Multi_results\Spontaneous\Degree_GraphsDenoised_WTThreshold_',num2str(idx_gen),'_NoEyes.csv'),h');
end

figure;
edges=[0:1:50];
for idx_gen=[3 1 2]
    histogram(Graph_info_perGen_denoised(idx_gen).degree,edges,'normalization','probability');hold on;
end

edges=[0:1:50];
for idx_gen=[3 1 2]
    h=[];
    for fish_nb=1:length(Graph_info_perGen_denoised(idx_gen).degree_ind)
        h(fish_nb,1)=nanmean(Graph_info_perGen_denoised(idx_gen).degree_ind{fish_nb});
        h(fish_nb,2)=nanstd(Graph_info_perGen_denoised(idx_gen).degree_ind{fish_nb});
    end
    csvwrite(strcat('\\shares.rdm.uq.edu.au\scottlab-q0291\Lena\Multi_results\Spontaneous\Degree_MeanSTD_GraphsDenoised_WTThreshold_',num2str(idx_gen),'_NoEyes.csv'),h');
end

Fighandle=figure;
set(Fighandle, 'Position', [50,50, 1500, 700]);
for idx_gen=[1 2]
    idx_temp=find(ismember(FishList,Genome{idx_gen}));
    if idx_gen==2;
        idx_gen=4;
        idx_temp2=find(ismember(FishList,Genome{idx_gen}));
        idx_temp=union(idx_temp,idx_temp2);
        idx_gen=2;
    end
    for Fish_nb=1:length(idx_temp)
        if ~isempty(ROI_fish{idx_temp(Fish_nb),1})
            Graph_temp=Graphs_correlation_WTthr_denoised{idx_temp(Fish_nb)};
            plot(Zbrain_brainMask2D(:,1),Zbrain_brainMask2D(:,2),'k');hold on;
            ROIs_temp=ROI_fish{idx_temp(Fish_nb),1};
            plot(Graph_temp,'EdgeColor','k','NodeColor',[0.5 0 0],'LineStyle','-','Marker','.','MarkerSize',5,'XData',ROIs_temp(ROI_fish{idx_temp(Fish_nb),2},1),'YData',ROIs_temp(ROI_fish{idx_temp(Fish_nb),2},2));hold off;
            print(Fighandle,strcat('\\shares.rdm.uq.edu.au\scottlab-q0291\Lena\Multi_results\Spontaneous\','Graphs_',num2str(idx_gen),'_Fish_',num2str(Fish_nb),'_NoEyes'),'-dpng','-r0');
        end
    end
end

edges={[0:0.01:2],[0 75 150 50000]};
CorrVsDist_SML=cell(2,length(FishList));
progressbar;
for i=1:length(FishList)   
    progressbar(i/length(FishList));
    if ~isempty(Corr_matrices_denoised{i})
        [CorrVsDist_SML{1,i}, CorrVsDist_SML{2,i}] = hist3([Corr_matrices_denoised{i}' Dist_matrices{i}'],'Edges',edges);    
    end
end

for idx_gen=[3 1 2]
    idx_temp=find(ismember(FishList,Genome{idx_gen}));
    if idx_gen==2;
        idx_gen=4;
        idx_temp2=find(ismember(FishList,Genome{idx_gen}));
        idx_temp=union(idx_temp,idx_temp2);
        idx_gen=2;
    end    
    Histogram_CorrDistBin=[];
    Mean_STD=[];
    for fish_nb=1:length(idx_temp)
        if ~isempty(CorrVsDist_SML{1,idx_temp(fish_nb)})
            temp=CorrVsDist_SML{1,idx_temp(fish_nb)};            
            Histogram_CorrDistBin=[Histogram_CorrDistBin temp(:,1:3)];
            temp=1-Corr_matrices_denoised{idx_temp(fish_nb)};
            idx_short=find(Dist_matrices{idx_temp(fish_nb)}<=75);
            idx_med=find(75<Dist_matrices{idx_temp(fish_nb)}<=150);
            idx_long=find(150<Dist_matrices{idx_temp(fish_nb)});
            Mean_STD=[Mean_STD nanmean(temp(idx_short)) nanstd(temp(idx_short))];
            Mean_STD=[Mean_STD nanmean(temp(idx_med)) nanstd(temp(idx_med))];
            Mean_STD=[Mean_STD nanmean(temp(idx_long)) nanstd(temp(idx_long))];
        end
    end
    csvwrite(strcat('\\shares.rdm.uq.edu.au\scottlab-q0291\Lena\Multi_results\Spontaneous\CorrVsDist_3bins_',num2str(idx_gen),'.csv'),Histogram_CorrDistBin');
    csvwrite(strcat('\\shares.rdm.uq.edu.au\scottlab-q0291\Lena\Multi_results\Spontaneous\CorrVsDist_Mean_STD_',num2str(idx_gen),'.csv'),Mean_STD');
end


for idx_gen=[3 1 2]
    idx_temp=find(ismember(FishList,Genome{idx_gen}));
    if idx_gen==2;
        idx_gen=4;
        idx_temp2=find(ismember(FishList,Genome{idx_gen}));
        idx_temp=union(idx_temp,idx_temp2);
        idx_gen=2;
    end        
    Mean_STD=[];
    for fish_nb=1:length(idx_temp)
        if ~isempty(CorrVsDist_SML{1,idx_temp(fish_nb)})
            temp=CorrVsDist_SML{1,idx_temp(fish_nb)};       
            temp=1-Corr_matrices_denoised{idx_temp(fish_nb)};            
            Mean_STD=[Mean_STD nanmean(temp) nanstd(temp)];            
        end
    end    
    csvwrite(strcat('\\shares.rdm.uq.edu.au\scottlab-q0291\Lena\Multi_results\Spontaneous\CorrVsDist_Mean_STD_',num2str(idx_gen),'_noBin.csv'),Mean_STD');
end

edges={[0:0.01:2],[0 75 150 50000]};
CorrVsDist_SML=cell(2,length(FishList));
progressbar;
for i=1:length(FishList)   
    progressbar(i/length(FishList));
    if ~isempty(Corr_matrices_denoised{i})
        [CorrVsDist_SML{1,i}, CorrVsDist_SML{2,i}] = hist3([Corr_matrices_denoised{i}' Dist_matrices{i}'],'Edges',edges);    
    end
end

for idx_gen=[3 1 2]
    idx_temp=find(ismember(FishList,Genome{idx_gen}));
    if idx_gen==2;
        idx_gen=4;
        idx_temp2=find(ismember(FishList,Genome{idx_gen}));
        idx_temp=union(idx_temp,idx_temp2);
        idx_gen=2;
    end    
    Histogram_CorrDistBin=[];
    Mean_STD=[];
    for fish_nb=1:length(idx_temp)
        if ~isempty(CorrVsDist_SML{1,idx_temp(fish_nb)})
            temp=CorrVsDist_SML{1,idx_temp(fish_nb)};            
            Histogram_CorrDistBin=[Histogram_CorrDistBin temp(:,1:3)];
            temp=1-Corr_matrices_denoised{idx_temp(fish_nb)};
            idx_short=find(Dist_matrices{idx_temp(fish_nb)}<=75);
            idx_med=find(75<Dist_matrices{idx_temp(fish_nb)}<=150);
            idx_long=find(150<Dist_matrices{idx_temp(fish_nb)});
            Mean_STD=[Mean_STD nanmean(temp(idx_short)) nanstd(temp(idx_short))];
            Mean_STD=[Mean_STD nanmean(temp(idx_med)) nanstd(temp(idx_med))];
            Mean_STD=[Mean_STD nanmean(temp(idx_long)) nanstd(temp(idx_long))];
        end
    end
    csvwrite(strcat('\\shares.rdm.uq.edu.au\scottlab-q0291\Lena\Multi_results\Spontaneous\CorrVsDist_3bins_',num2str(idx_gen),'.csv'),Histogram_CorrDistBin');
    csvwrite(strcat('\\shares.rdm.uq.edu.au\scottlab-q0291\Lena\Multi_results\Spontaneous\CorrVsDist_Mean_STD_',num2str(idx_gen),'.csv'),Mean_STD');
end