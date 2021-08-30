results_ENS=struct();
results_ENS_raw=struct();
for treat=1:3
    for dpf=3:7
        idx_temp=ZS_guts{treat,dpf-2,2};
        for fish_nb=1:max(unique(idx_Fish(idx_temp)))
            ZS_temp=ZS(idx_temp(idx_Fish(idx_temp)==fish_nb),:);
            if ~isempty(ZS_temp)
                Sm=smooth_ZS_fish{treat,dpf-2,fish_nb,3};
                pks_smooth=cell(1,size(Sm,1));
                locs_smooth=cell(1,size(Sm,1));                
                Periodicity_smooth=nan(size(Sm,1),1);
                pks_raw=cell(1,size(Sm,1));
                locs_raw=wcell(1,size(Sm,1));                
                Periodicity_raw=nan(size(Sm,1),1);
                for roi_nb=1:size(ZS_temp,1)
                    [pks_smooth{roi_nb},locs_smooth{roi_nb}] = findpeaks(Sm(roi_nb,:),'MinPeakProminence',1.5,'MinPeakDistance',5);
                    [pks_raw{roi_nb},locs_raw{roi_nb}] = findpeaks(ZS_temp(roi_nb,:),'MinPeakProminence',1.5,'MinPeakDistance',5);
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


edge=[0:2:60];
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

treat=3;
h=nan(length(edge)-1,5*8);
for dpf=3:7
    for fish_nb=1:8
        h(:,fish_nb+(dpf-3)*8)=histcounts(results_ENS_raw(treat,dpf-2,fish_nb).Periodicity,edge,'Normalization','probability')';
    end
end

peak_nb=nan(3,5,8);
for treat=1:3
    for dpf=3:7
        idx_temp=ZS_guts{treat,dpf-2,2};
        for fish_nb=1:max(unique(idx_Fish(idx_temp)))
            ZS_temp=DF(idx_temp(idx_Fish(idx_temp)==fish_nb),:);            
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
treat=3;
for dpf=3:7
    idx_temp=ZS_guts{treat,dpf-2,2};
    for fish_nb=1:max(unique(idx_Fish(idx_temp)))
        ZS_temp=DF(idx_temp(idx_Fish(idx_temp)==fish_nb),:);
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