LM_FMR=[];
parfor i=1:length(ZS)
    %mdl=stepwiselm(Regressor',DataSet(i,:),'linear','Criterion','adjrsquared','Intercept',false,'Upper','interactions','Verbose',0);
    mdl=fitlm(Looms(1,:)',ZS(i,:),'interactions');
    LM_FMR(i).coef=mdl.Coefficients;
    %ModelResultsSeg_ZS(i).MSE=mdl.MSE;
    %ModelResultsSeg_ZS(i).Fitted=mdl.Fitted;
    LM_FMR(i).rsquared=mdl.Rsquared.Adjusted;
end
idx_rsq=find([LM_FMR.rsquared]>0.05);
ZS_rsq=ZS(idx_rsq,:);
figure;imagesc(ZS_rsq, [-0.5 4]);colormap hot;

options = statset('UseParallel',1); [idxKmeans_ZS_rsq Cmap_ZS_rsq]=kmeans(ZS_rsq,10,'Options',options,'Distance','cityblock','Replicates',3,'MaxIter',1000,'Display','final');