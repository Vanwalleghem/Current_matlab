%Exclude blurry fish
% Genome{1}=[201810222 201810223 201810224 201810232 201810234 201810292 201811054 201811062 201811127];%mutant
Genome{1}=[201810223, 201810224, 201810234, 201811054, 201811056, 201811062, 201811127];%mutant
Genome{2}=[201810225 201810291 201810294 201810301 201811053 201811224];
% Genome{3}=[201810233 201810297 201810302 201811055 201811064 201811124];
Genome{3}=[ 201810233 201810297 201811055 201811064 201811124];
Genome{4}=[201810226 201810293 201810295  201811052 201811063 201811125 201811191 201811196];%Hets

Corr_matrices_noisy={};
progressbar;
for idx_gen=[3 1]
    idx_Fish_temp=find(ismember(20180000+FishList,Genome{idx_gen}));
    for fish_nb=1:length(idx_Fish_temp)
        if ~isempty(ROI_fish{idx_Fish_temp(fish_nb),2})    
            idx_temp=find(idx_Fish==FishList(idx_Fish_temp(fish_nb)));
            if length(idx_temp)==length(ROI_fish{idx_Fish_temp(fish_nb),2})
                isInBrain=find(ROI_fish{idx_Fish_temp(fish_nb),2});
                Corr_matrices_noisy{idx_gen,fish_nb}=pdist(ZS_CN(idx_temp(isInBrain),:),'correlation');
            else
                idx_Fish_temp(fish_nb)
                break
            end
        end
        progressbar(fish_nb/length(idx_Fish_temp));
    end
end

Dist_matrices={};
progressbar;
for idx_gen=[3 1]
    idx_Fish_temp=find(ismember(FishList,Genome{idx_gen}));
    for fish_nb=1:length(idx_Fish_temp)
        ROI_temp=ROI_fish{idx_Fish_temp(fish_nb),1};
        isInBrain=find(ROI_fish{idx_Fish_temp(fish_nb),2});
        Dist_matrices{idx_gen,fish_nb}=pdist(ROI_temp(isInBrain,:));
        progressbar(fish_nb/length(idx_Fish_temp));
    end
end

Max_distance_ROIs=0;
for idx_gen=[3 1]
    idx_Fish_temp=find(ismember(FishList,Genome{idx_gen}));
    for fish_nb=1:length(idx_Fish_temp)    
        Max_temp=max(Dist_matrices{idx_gen,fish_nb});
        if Max_distance_ROIs<Max_temp
            Max_distance_ROIs=Max_temp;
        end
    end
end

edges={[0:0.005:2],[0:0.005:2]};
CorrVsDist={};
for idx_gen=[3 1]
    idx_Fish_temp=find(ismember(FishList,Genome{idx_gen}));
    for fish_nb=1:length(idx_Fish_temp)
        progressbar(fish_nb/length(idx_Fish_temp));
        CorrVsDist{idx_gen,fish_nb} = hist3([Corr_matrices_noisy{idx_gen,fish_nb}' Dist_matrices{idx_gen,fish_nb}'/(Max_distance_ROIs/2)],'Edges',edges);
    end
end

edges={[0:0.005:2],[0:0.005:2]};
CorrVsDist={};
for idx_gen=[3 1]
    idx_Fish_temp=find(ismember(FishList,Genome{idx_gen}));
    for fish_nb=1:length(idx_Fish_temp)
        progressbar(fish_nb/length(idx_Fish_temp));
        Corr_temp=Corr_matrices_noisy{idx_gen,fish_nb}';
        Dist_temp=Dist_matrices{idx_gen,fish_nb}'/(Max_distance_ROIs/2);
        idx_relevant=union(find(Corr_temp<(mean(Corr_temp)-2*std(Corr_temp))),find(Corr_temp>(mean(Corr_temp)+2*std(Corr_temp))));
        CorrVsDist{idx_gen,fish_nb} = hist3([Corr_temp(idx_relevant) Dist_temp(idx_relevant)],'Edges',edges);        
    end
end

edges=[0:5:500];
CorrVsDist_values={};
for idx_gen=[3 1]
    idx_Fish_temp=find(ismember(FishList,Genome{idx_gen}));
    values_temp=nan(length(idx_Fish_temp),length(edges)-1);
    for fish_nb=1:length(idx_Fish_temp)
        progressbar(fish_nb/length(idx_Fish_temp));        
        Corr_temp=Corr_matrices_noisy{idx_gen,fish_nb}';
        Dist_temp=Dist_matrices{idx_gen,fish_nb}';
        idx_relevant=union(find(Corr_temp<(mean(Corr_temp)-2*std(Corr_temp))),find(Corr_temp>(mean(Corr_temp)+2*std(Corr_temp))));
        for bin_nb=2:length(edges)
            idx_dist=edges(bin_nb-1)<Dist_temp(idx_relevant) & Dist_temp(idx_relevant)<edges(bin_nb);
            values_temp(fish_nb,bin_nb-1) = mean(Corr_temp(idx_relevant(idx_dist)));        
        end
        CorrVsDist_values{idx_gen}=values_temp;
    end
end

rows=3;yplot=3;counter=1;
for idx_gen=[3 1]
    Fighandle=figure;
    set(Fighandle, 'Position', [200, 200, 1500, 1500]);
    idx_temp=find(ismember(FishList,Genome{idx_gen}));
    ha = tight_subplot(rows,yplot,[.02 .02],[.02 .02],[.02 .02]);
    h=[];
    for fish_nb=1:length(idx_temp)
        axes(ha(fish_nb));
        imagesc(CorrVsDist{idx_gen,fish_nb})
    end
    print(Fighandle,strcat('\\shares.rdm.uq.edu.au\scottlab-q0291\Lena\Multi_results\Spontaneous\','WTvsFMR_CorrelationVsDistance_',num2str(idx_gen),'_NoEyes'),'-dpng','-r0');
end
clearvars rows counter i idx_temp idx_temp2 idx_gen Fighandle yplot
idx_gen=3;
idx_temp=find(ismember(FishList,Genome{idx_gen}));
CorrVDist_temp=zeros(size(CorrVsDist{idx_gen,1},1),size(CorrVsDist{idx_gen,1},1),length(idx_temp));
for fish_nb=1:length(idx_temp)
    temp=CorrVsDist{idx_gen,(fish_nb)};
    if ~isempty(temp)
        CorrVDist_temp(:,:,fish_nb)=temp/sum(temp(:));
    end
end
AVG_WT=squeeze(mean(CorrVDist_temp,3));
idx_gen=1;
idx_temp=find(ismember(FishList,Genome{idx_gen}));
CorrVDist_temp=zeros(size(CorrVsDist{1,1},1),size(CorrVsDist{1,1},1),length(idx_temp));
for fish_nb=1:length(idx_temp)
    temp=CorrVsDist{idx_gen,(fish_nb)};
    if ~isempty(temp)
        CorrVDist_temp(:,:,fish_nb)=temp/sum(temp(:));
    end
end
AVG_FMR=squeeze(mean(CorrVDist_temp,3));
Fighandle=figure;
set(Fighandle, 'Position', [200, 200, 3200, 800]);
ha = tight_subplot(1,2,[.02 .02],[.02 .02],[.02 .02]);
axes(ha(1));imagesc(AVG_WT,[0 8e-5]);
axes(ha(2));imagesc(AVG_FMR,[0 8e-5]);

Fighandle=figure;
set(Fighandle, 'Position', [200, 200, 1000, 1000]);
plot(sum(AVG_WT,1));hold on;plot(sum(AVG_FMR,1));

Fighandle=figure;
AVG_WT=CorrVsDist_values{3};
x=linspace(0,Max_distance_ROIs/2,size(AVG_WT,2));
set(Fighandle, 'Position', [200, 200, 1000, 1000]);
meanToPlot=mean(1-AVG_WT,1);
stdToPlot=std(1-AVG_WT,1,1);
H=shadedErrorBar(x, meanToPlot, stdToPlot);
H.mainLine.LineWidth=3;
H.mainLine.Color=[0 0 1];
H.patch.FaceColor=[0 0 1];
H.edge(1).Color=[0 0 1];
H.edge(2).Color=[0 0 1];
hold on;
AVG_FMR=CorrVsDist_values{1};
meanToPlot=mean(1-AVG_FMR,1);
stdToPlot=std(1-AVG_FMR,1,1);
H=shadedErrorBar(x, meanToPlot, stdToPlot);
H.mainLine.LineWidth=3;
H.mainLine.Color=[1 0.5 0];
H.patch.FaceColor=[1 0.5 0];
H.edge(1).Color=[1 0.5 0];
H.edge(2).Color=[1 0.5 0];

edges=[0:5:500];
CorrVsDist_values_All={};
for idx_gen=[3 1]
    idx_Fish_temp=find(ismember(FishList,Genome{idx_gen}));
    values_temp=nan(length(idx_Fish_temp),length(edges)-1);
    for fish_nb=1:length(idx_Fish_temp)
        progressbar(fish_nb/length(idx_Fish_temp));        
        Corr_temp=Corr_matrices_noisy{idx_gen,fish_nb}';
        Dist_temp=Dist_matrices{idx_gen,fish_nb}';
        
        for bin_nb=2:length(edges)
            idx_dist=edges(bin_nb-1)<Dist_temp() & Dist_temp()<edges(bin_nb);
            values_temp(fish_nb,bin_nb-1) = mean(Corr_temp((idx_dist)));        
        end
        CorrVsDist_values_All{idx_gen}=values_temp;
    end
end

Fighandle=figure;
AVG_WT=CorrVsDist_values_All{3};
x=linspace(0,Max_distance_ROIs/2,size(AVG_WT,2));
set(Fighandle, 'Position', [200, 200, 1000, 1000]);
meanToPlot=mean(1-AVG_WT,1);
stdToPlot=std(1-AVG_WT,1,1);
H=shadedErrorBar(x, meanToPlot, stdToPlot);
H.mainLine.LineWidth=3;
H.mainLine.Color=[0 0 1];
H.patch.FaceColor=[0 0 1];
H.edge(1).Color=[0 0 1];
H.edge(2).Color=[0 0 1];
hold on;
AVG_FMR=CorrVsDist_values_All{1};
meanToPlot=mean(1-AVG_FMR,1);
stdToPlot=std(1-AVG_FMR,1,1);
H=shadedErrorBar(x, meanToPlot, stdToPlot);
H.mainLine.LineWidth=3;
H.mainLine.Color=[1 0.5 0];
H.patch.FaceColor=[1 0.5 0];
H.edge(1).Color=[1 0.5 0];
H.edge(2).Color=[1 0.5 0];

FiringRate_noisy={};
progressbar;
for idx_gen=[3 1]
    idx_Fish_temp=find(ismember(FishList,Genome{idx_gen}));
    for fish_nb=1:length(idx_Fish_temp)
        if ~isempty(ROI_fish{idx_Fish_temp(fish_nb),2})
            idx_temp=find(idx_Fish==FishList(idx_Fish_temp(fish_nb)));
            if length(idx_temp)==length(ROI_fish{idx_Fish_temp(fish_nb),2})
                isInBrain=find(ROI_fish{idx_Fish_temp(fish_nb),2});
                ZS_temp=ZS_CN(idx_temp(isInBrain),:);
                temp=nan(size(ZS_temp,1),2);
                for roi_nb=1:length(temp)
                    [~,loc]=findpeaks(ZS_temp(roi_nb,:),'MinPeakDistance',5,'MinPeakProminence',1.5);
                    temp(roi_nb,1)=length(loc);
                    [~,loc]=findpeaks(-ZS_temp(roi_nb,:),'MinPeakDistance',5,'MinPeakHeight',1.5);
                    temp(roi_nb,2)=length(loc);
                    progressbar(roi_nb/length(temp),[]);
                end
                FiringRate_noisy{idx_gen,fish_nb}=temp;
            end
            progressbar(fish_nb/length(idx_Fish_temp));
        end
    end
end
