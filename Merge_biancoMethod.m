Merged_DF=ZS2_rsq;
counter=1;
progressbar % can be skipped if you haven't progressbar installed    
while counter<=size(Merged_DF,1) % can also be for counter=1:size(Dataset,1)    
    if isnan(Merged_DF(counter,1)) % if the timeserie has been merged, it's set to nan, so this prevents spending time on it
        max_corr=0;
        MergedIDX{counter}=[];
    else
        max_corr=1;    % hack to get the equivalent of a do... while loop (at least one iteration)
        MergedIDX{counter}=counter;
    end    
    while max_corr>0.75 % standard threshold of Bianco et al
        corr_temp=zeros(size(Merged_DF,1),1);        
        parfor j=1:size(Merged_DF,1)    % parallelized, can replace by for loop if not needed
        	if j>counter && ~isnan(Merged_DF(j,1))
                temp=corrcoef(Merged_DF(counter,:), Merged_DF(j,:));
                corr_temp(j)=temp(1,2);                                
            end
        end        
        [max_corr idx_max]=nanmax(corr_temp);
        if max_corr>0.75
            if max_corr>0.85 % shortcut to merge in one go all time series with >0.85 correlation, can be skipped or changed
                idx_max=find(corr_temp>0.85);
                MergedIDX{counter}=[MergedIDX{counter} idx_max'];
            else
                MergedIDX{counter}=[MergedIDX{counter} idx_max]; % merge things one by one, so it's long
            end
            Merged_DF(counter,:)=nanmean(ZS2_rsq(MergedIDX{counter},:),1); %average from dataset, easier than adding one at a time with a factor                
            Merged_DF(idx_max,:)=nan;
            corr_temp(idx_max)=nan;            
        end
    end
    counter=counter+1;
    progressbar(counter/size(Merged_DF,1)); % can be skipped if you haven't progressbar installed    
end
clearvars i j idx_max counter max_corr temp

number_merge=[];
for i=1:length(MergedIDX)
    number_merge(i)=length(MergedIDX{i});    
end

clusters=Merged_DF(find(number_merge>100),:);

figure;
for i=1:size(clusters,1)
    plot(clusters(i,:)),pause
end

figure;imagesc(ZS2_rsq(MergedIDX{1},:));