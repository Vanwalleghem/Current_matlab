File_list=dir('suite2p*');
MatFiles_names={File_list.name};

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

[~,idx_sort]=sortrows(names',[2 3]);
File_list=File_list(idx_sort'); %This orders the Matlab files per fish, it helps with the indexing of the ROIs after ANTs
File_list = rmfield(File_list,'isdir');File_list = rmfield(File_list,'folder');File_list = rmfield(File_list,'bytes');File_list = rmfield(File_list,'date');
File_list = rmfield(File_list,'datenum');

Suite2p_traces=[];IsCell=[];Suite2p_centroids=[];Numbers=[];
idx_DPF=[];idx_Fish=[];idx_Treatment=[];
for i=1:length(File_list)    
    load(strcat(File_list(i).name,'\plane0\Fall.mat'),'F');    
    load(strcat(File_list(i).name,'\plane0\Fall.mat'),'stat');        
    load(strcat(File_list(i).name,'\plane0\Fall.mat'),'iscell');        
    Centroids=zeros(size(F,1),2);
    for ij=1:length(stat)
        Centroids(ij,:)=stat{ij}.med;        
    end    
    if i==1
        Suite2p_traces=F;
        IsCell=iscell(:,2);
        Suite2p_centroids=Centroids;
        Numbers=size(F,1);   
        idx_DPF(1:Numbers)=names(2,i);
        idx_Fish(1:Numbers)=names(1,i);
        idx_Treatment(1:Numbers)=names(3,i);
    else        
        if size(F,2)<1200
            F=padarray(F,[0 1200-size(F,2)],'replicate','post');            
        end
        Suite2p_traces=vertcat(F,Suite2p_traces);
        IsCell=vertcat(iscell(:,2),IsCell);
        Suite2p_centroids=vertcat(Centroids,Suite2p_centroids);
        Numbers=vertcat(size(F,1),Numbers);
        idx_DPF(Numbers(i-1)+1:Numbers(i))=names(2,i);
        idx_Fish(Numbers(i-1)+1:Numbers(i))=names(1,i);
        idx_Treatment(Numbers(i-1)+1:Numbers(i))=names(3,i);
    end        
end
close all;
clearvars i ij F Centroids ans ij iscell stat

Numbers=[0; cumsum(Numbers)];

ZS=zscore(Suite2p_traces,1,2);

save('RAW_suite2p','-v7.3');

figure;imagesc(ZS(IsCell>0.05,:),[-0.5 3]);colormap hot;

fid = fopen('Gut limits');
Limits=inf(2,length(File_list));
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


idx_DPF=[];idx_Fish=[];idx_Treatment=[];
for i=1:length(File_list)   
    idx_DPF(Numbers(i)+1:Numbers(i+1))=names(2,i);
    idx_Fish(Numbers(i)+1:Numbers(i+1))=names(1,i);
    idx_Treatment(Numbers(i)+1:Numbers(i+1))=names(3,i);            
end
close all;
clearvars i ij F Centroids ans ij iscell stat

%You can also build an index of the Fish and the plane of each ROI
idx_Gut=nan(size(ZS,1),1);
for i=1:length(File_list)
    ROI_temp=Suite2p_centroids(Numbers(i)+1:Numbers(i+1),:);
    idx_ROI=(ROI_temp(:,1)<Limits(1,i) & ROI_temp(:,1)>20) & (ROI_temp(:,2)<Limits(2,i) & ROI_temp(:,2)>20);
    
    idx_Gut(Numbers(i)+1:Numbers(i+1))=idx_ROI;
    
end
clearvars i Fish Plane name counter
idx_Gut=boolean(idx_Gut);
idx_Cell=IsCell>0.05;

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
        idx_temp=idx_temp(idx_Gut(idx_temp) & idx_Cell(idx_temp));
        imagesc(ZS(idx_temp(randperm(length(idx_temp))),:),[-0.5 4]);colormap(hot);
        counter=counter+1;
        set(gca,'XTickLabel',[]);set(gca,'YTickLabel',[]);
    end    
end

ZS_guts={};
for treat=1:3
    for dpf=3:7  
        idx_temp=find(idx_DPF==dpf & idx_Treatment==treat);
        idx_temp=idx_temp(idx_Gut(idx_temp) & idx_Cell(idx_temp));
        ZS_guts{treat,dpf-2,1}=detrend(ZS(idx_temp,:)')'; 
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
            ZS_temp=detrend(ZS(idx_temp(idx_Fish(idx_temp)==fish_nb),:)')';
            if ~isempty(ZS_temp)
                ops.iPC=1:min([size(ZS_temp,1) 1200]);
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

results_ENS=struct();
results_ENS_raw=struct();
for treat=1:3
    for dpf=3:7
        idx_temp=ZS_guts{treat,dpf-2,2};
        for fish_nb=1:max(unique(idx_Fish(idx_temp)))
            ZS_temp=detrend(ZS(idx_temp(idx_Fish(idx_temp)==fish_nb),:)')';
            if ~isempty(ZS_temp)
                Sm=smooth_ZS_fish{treat,dpf-2,fish_nb,3};
                pks_smooth=cell(1,size(Sm,1));
                locs_smooth=cell(1,size(Sm,1));                
                Periodicity_smooth=nan(size(Sm,1),1);
                pks_raw=cell(1,size(Sm,1));
                locs_raw=cell(1,size(Sm,1));                
                Periodicity_raw=nan(size(Sm,1),1);
                for roi_nb=1:size(ZS_temp,1)
                    [pks_smooth{roi_nb},locs_smooth{roi_nb}] = findpeaks(Sm(roi_nb,:),'MinPeakProminence',1,'MinPeakDistance',5);
                    [pks_raw{roi_nb},locs_raw{roi_nb}] = findpeaks(ZS_temp(roi_nb,:),'MinPeakProminence',1,'MinPeakDistance',5);
                    if length(locs_smooth{roi_nb})>5
                        locs_temp=locs_smooth{roi_nb}(2:end)-locs_smooth{roi_nb}(1:end-1);
                        locs_temp=locs_temp/2;
                        Periodicity_smooth(roi_nb)=mean(locs_temp);
                    end
                    if length(locs_raw{roi_nb})>5
                        locs_temp=locs_raw{roi_nb}(2:end)-locs_raw{roi_nb}(1:end-1);
                        locs_temp=locs_temp/2;
                        Periodicity_raw(roi_nb)=mean(locs_temp);
                    end
                end
                results_ENS(treat,dpf-2,fish_nb).peaks=pks_smooth;
                results_ENS(treat,dpf-2,fish_nb).locs=locs_smooth;
                results_ENS(treat,dpf-2,fish_nb).Periodicity=Periodicity_smooth;
                results_ENS_raw(treat,dpf-2,fish_nb).peaks=pks_raw;
                results_ENS_raw(treat,dpf-2,fish_nb).locs=locs_raw;
                results_ENS_raw(treat,dpf-2,fish_nb).Periodicity=Periodicity_raw;
            end
        end
    end
end

edge=[0:2:40];
treat=3;
h=nan(length(edge)-1,5*8);
for dpf=3:7
    for fish_nb=1:8
        h(:,fish_nb+(dpf-3)*8)=histcounts(results_ENS_raw(treat,dpf-2,fish_nb).Periodicity,edge,'Normalization','probability')';
    end
end


PrismTemp=nan(3,5*8);
for treat=1:3
    for dpf=3:7
        for fish_nb=1:8
            PrismTemp(treat,fish_nb+(dpf-3)*8)=nanmean(results_ENS_raw(treat,dpf-2,fish_nb).Periodicity);
        end
    end
end
PrismTemp(PrismTemp==0)=nan;

peak_nb=nan(3,5,8);
for treat=1:3
    for dpf=3:7
        idx_temp=ZS_guts{treat,dpf-2,2};
        for fish_nb=1:max(unique(idx_Fish(idx_temp)))
            ZS_temp=detrend(ZS(idx_temp(idx_Fish(idx_temp)==fish_nb),:)')';  
            if ~isempty(ZS_temp)
                peak_temp_nb=0;
                peak_temp=results_ENS_raw(treat,dpf-2,fish_nb).peaks;
                for i=1:length(peak_temp)
                    peak_temp_nb=peak_temp_nb+length(peak_temp{i});
                end
                peak_nb(treat,dpf-2,fish_nb)=peak_temp_nb/length(peak_temp);
            end            
        end
    end
end

PrismTemp=nan(3,5*8);
for treat=1:3
    for dpf=3:7
        for fish_nb=1:8
            PrismTemp(treat,fish_nb+(dpf-3)*8)=peak_nb(treat,dpf-2,fish_nb);
        end
    end
end
PrismTemp=PrismTemp/10;

edge=[0:4:60];
h=nan(length(edge)-1,5*8);
treat=1;
for dpf=3:7
    idx_temp=ZS_guts{treat,dpf-2,2};
    for fish_nb=1:max(unique(idx_Fish(idx_temp)))
        ZS_temp=detrend(ZS(idx_temp(idx_Fish(idx_temp)==fish_nb),:)')';
        if ~isempty(ZS_temp)
            peak_temp_nb=[];
            peak_temp=results_ENS_raw(treat,dpf-2,fish_nb).peaks;
            for i=1:length(peak_temp)
                peak_temp_nb(i)=length(peak_temp{i});
            end            
            h(:,fish_nb+(dpf-3)*8)=histcounts(peak_temp_nb,edge,'Normalization','probability')';
        end        
    end
end