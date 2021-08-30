CorrelationMatrices_minusFish=nan(size(Mean_allNodes_perFish,1),size(Mean_allNodes_perFish,1),size(Mean_allNodes_perFish,2));
for fish_nb=1:size(Mean_allNodes_perFish,2)
    fish_idx=[1:1:size(Mean_allNodes_perFish,2)];
    fish_idx(fish_nb)=[];
    CorrelationMatrices_minusFish(:,:,fish_nb)=squareform(pdist(squeeze(nanmean(Mean_allNodes_perFish(:,fish_idx,:),2)),'correlation'));
end

Strength_perClusterAndFish=zeros(8,size(Mean_allNodes_perFish,2));
for fish_nb=1:size(Mean_allNodes_perFish,2)
    Corr_temp=squeeze(CorrelationMatrices_minusFish(:,:,fish_nb));
    Corr_temp=abs(1-Corr_temp);
    Corr_temp=weight_conversion(Corr_temp,'autofix');
    Corr_temp=weight_conversion(Corr_temp,'normalize');
    graphTemp=threshold_proportional(Corr_temp,0.25);
    DensTemp=strengths_und(graphTemp);
    for i=1:8
        Strength_perClusterAndFish(i,fish_nb)=mean(DensTemp(Node_ID==i));
    end
end


Degree_perClusterAndFish=zeros(8,size(Mean_allNodes_perFish,2));
for fish_nb=1:size(Mean_allNodes_perFish,2)
    Corr_temp=squeeze(CorrelationMatrices_minusFish(:,:,fish_nb));
    Corr_temp=abs(1-Corr_temp);
    Corr_temp=weight_conversion(Corr_temp,'autofix');
    Corr_temp=weight_conversion(Corr_temp,'normalize');
    graphTemp=threshold_proportional(Corr_temp,0.25);
    DensTemp=degrees_und(graphTemp);
    for i=1:8
        Degree_perClusterAndFish(i,fish_nb)=mean(DensTemp(Node_ID==i));
    end
end

Participation_perClusterAndFish=zeros(8,size(Mean_allNodes_perFish,2));
for fish_nb=1:size(Mean_allNodes_perFish,2)
    Corr_temp=squeeze(CorrelationMatrices_minusFish(:,:,fish_nb));
    Corr_temp=abs(1-Corr_temp);
    Corr_temp=weight_conversion(Corr_temp,'autofix');
    Corr_temp=weight_conversion(Corr_temp,'normalize');
    graphTemp=threshold_proportional(Corr_temp,0.25);
    DensTemp=participation_coef(graphTemp,Node_ID);    
    for i=1:8
        Participation_perClusterAndFish(i,fish_nb)=mean(DensTemp(Node_ID==i));
    end
end