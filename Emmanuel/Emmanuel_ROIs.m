load('D:\Pictures\processed\Emmanuel\f20_r2050_CL3_idxs_poster.mat')

list_fish=unique(idx_Fish);

CSV_Files=dir('_BetterWarped*.csv');
ROIs=struct();
for i=1:length(CSV_Files);
    temp=csvread(CSV_Files(i).name,1);    
    Fishname=regexp(CSV_Files(i).name,'_BetterWarped(.+).csv','tokens');Fishname=Fishname{1}{1};
    if ismember(str2num(Fishname),list_fish)
        ROIs(str2num(Fishname)).name=Fishname;    
        ROIs(str2num(Fishname)).coord=temp(:,1:3);
    end
end
clearvars i temp CSV_Files Fishname

% MatFiles2=dir('*analysis_matlab.mat');
% MatFiles2_names={MatFiles2.name}';
% ROI_fish_goodIdx={};
% for i=1:length(list_fish)
%     fish_name=strcat('_fish',num2str(list_fish(i)),'_');
%     IndexC=strfind(MatFiles2_names,fish_name);
%     IndexC=find(not(cellfun('isempty', IndexC)));
%     ROI_goodIdx=[];
%     idx_temp=0;
%     for j=1:length(IndexC)
%          name=strcat(MatFiles2(IndexC(j)).name);
%          F=load(name, 'idx_components');
%          F=F.idx_components+1;         
%          b=load(name, 'Baseline');
%          b=b.Baseline;
%          if j==1
%              ROI_goodIdx=F;
%          else             
%              ROI_goodIdx=[ROI_goodIdx F+idx_temp];
%          end         
%          idx_temp=idx_temp+length(b);
%     end
%     ROI_fish_goodIdx{i}=ROI_goodIdx;
% end
% clearvars i j ROI_goodIdx b F name IndexC fish_name idx_temp

ROI_pool=[];
for fish_nb=1:length(list_fish)    
    ROI_temp=ROIs(list_fish(fish_nb)).coord;    
    ROI_pool=[ROI_pool;ROI_temp];
end

Fish_list=list_fish;

%Sort ROIs according to Matfiles order
Sort_ROIs=[];
MatFiles_names={MatFiles.name};
for fish_nb=1:length(Fish_list)
    fish_name=strcat('_fish',num2str(list_fish(fish_nb)),'_');
    IndexC=strfind(MatFiles_names,fish_name);
    MatFiles_fish = find(not(cellfun('isempty', IndexC)));
    for file_nb=1:length(MatFiles_fish)
        if MatFiles_fish(file_nb)==1
            numbersForROIs=[1 MatFiles(1).GoodNumber];
        else
            numbersForROIs=[MatFiles(MatFiles_fish(file_nb)-1).GoodNumber+1 MatFiles(MatFiles_fish(file_nb)).GoodNumber];
        end
        Sort_ROIs=[Sort_ROIs numbersForROIs(1):1:numbersForROIs(2)];        
    end
end
clearvars slice fish_nb roi_nb ROI Centroids IndexC file_nb

ROI_fish(Sort_ROIs,:)=ROI_pool;
ROI_fish=ROI_fish(idx_coef_rsq,:);

% %Hack while the ROIs aren't fixed for the whole strfind issue
% idx_Fish_temp=idx_Fish(idx_coef_rsq);
% idx_Fish_temp=find(idx_Fish_temp<10);
% idxKmeans_final=idxKmeans1_coef_rsq;
% idxKmeans_final(idx_Fish_temp)=0;


for i=1:length(GoodBetas_ZS2)
    idx_temp=find(idxKmeans1_coef_rsq==GoodBetas_ZS2(i));    
    ROI_temp=ROI_fish(idx_temp,:);    
    ROI_temp(:,3)=ROI_temp(:,3)*2;
    ROI_temp(:,4)=1;
    filename=strcat('_S60_',num2str(i),'.csv');
    csvwrite(filename,ROI_temp);
end

%All brain regions
load('D:\Pictures\processed\Itia\_BrainRegAndROIs.mat', 'Zbrain_Masks')
Zbrain_AllMask=vertcat(Zbrain_Masks{[1:1:77 79:1:294],3});
Zbrain_AllMask=unique(Zbrain_AllMask,'rows');
%Removing the eyes




% % Create rotation matrix
% theta = 180; % to rotate 90 counterclockwise
% R = [cosd(theta) -sind(theta); sind(theta) cosd(theta)];
% 
% ROI_temp=round(ROI_fish);
% ROI_rot=R*ROI_temp(:,1:2)';
% ROI_rot=ROI_rot';
% ROI_rot(:,1)=ROI_rot(:,1)+621;
% ROI_rot(:,2)=ROI_rot(:,2)+1406;
% ROI_rot(:,3)=ROI_temp(:,3);
% ROI_rot_correct=[];
% ROI_rot_correct(:,2)=ROI_rot(:,1);
% ROI_rot_correct(:,1)=ROI_rot(:,2);
% ROI_rot_correct(:,3)=ROI_rot(:,3)-10;
ROI_temp=round(ROI_fish);
ROI_rot_correct=ROI_temp(:,[2 1 3]);


idx_brain=ismember(ROI_rot_correct,Zbrain_AllMask,'rows');

idxKmeans_final=idxKmeans1_coef_rsq;
idxKmeans_final(~idx_brain)=0;

GoodBetas_ZS2=goodmaps([3 2 1]);
for i=1:length(GoodBetas_ZS2)
    idx_temp=find(idxKmeans_final==GoodBetas_ZS2(i));    
    ROI_temp=ROI_fish(idx_temp,:);    
    ROI_temp(:,3)=ROI_temp(:,3)*2;
    ROI_temp(:,4)=1;
    filename=strcat('_F20_better_',num2str(i),'.csv');
    csvwrite(filename,ROI_temp);
end

% progressbar();
% for roi_nb=1:size(ROI_fish,1)
%     ROI_temp=round(ROI_fish(roi_nb,:));    
%     if ~ismember(ROI_temp,Zbrain_AllMask,'rows');
%         idxKmeans_final(roi_nb)=0;    
%     end
%    progressbar(roi_nb/size(ROI_fish,1));
% end
% clearvars GoodROIs Mask counter idx_temp Coords IsInMask IsInMask_temp