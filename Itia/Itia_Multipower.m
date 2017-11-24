Stimuli=zeros(6,size(ZS,2));
GCaMP6=[0,1.69644104899772,5.13796058542217,8.27886020152244,10.3756715204800,11.8173714529814,12.2425184714093,10.8571417354877,8.80831829681196,6.91339112244670,5.46959264663869,4.30868766622567,3.42533619066766,2.75378443486879,2.18017250852183,1.72816235135824,1.32732537295463,1.00684435500268,0.730210038304555,0.530242444093118,0.362253250339685,0.227668255288566,0.0869242416152502,0.000718266708050853,-0.0828334873368325]';
idxStart=58;
counter=0;
for i=1:18
    if mod(i,3)==1
        counter=counter+1;
    end
    Stimuli(counter,(idxStart+(i-1)*40):(idxStart+(i-1)*40)+size(GCaMP6,1)-1)=GCaMP6;
end

ModelMultipower=[];
parfor i=1:length(ZS)
    %mdl=stepwiselm(Regressor',DataSet(i,:),'linear','Criterion','adjrsquared','Intercept',false,'Upper','interactions','Verbose',0);
    mdl=fitlm(Stimuli',ZS(i,:),'interactions');
    ModelMultipower(i).coef=mdl.Coefficients;
    %ModelResultsSeg_ZS(i).MSE=mdl.MSE;
    %ModelResultsSeg_ZS(i).Fitted=mdl.Fitted;
    ModelMultipower(i).rsquared=mdl.Rsquared.Adjusted;
end

coefficients={};
for idx=1:length(ModelMultipower)
    coef=[ModelMultipower(idx).coef];
    temp=coef.Properties.RowNames;temp=regexp(temp,'x(\d+)','tokens');
    if ~isempty(temp)
        %temp=[temp{:}];temp=[temp{:}];temp=[temp{:}];%temp=str2num(temp);
        for coef_idx=2:height(coef)
            if coef.pValue(coef_idx)<0.05
                coefficients{idx,str2num(temp{coef_idx}{1}{1})}=coef.Estimate(coef_idx);
            end
        end
    end
end
idxempty=cellfun('isempty',coefficients);
coefficients(idxempty)={0};
clearvars idxempty idx coef_idx coef temp
coefficients=cell2mat(coefficients);

ZS_AVG=zeros(size(ZS,1),246);
ZS_STD=zeros(size(ZS,1),246);
parfor idx_ZS=1:size(ZS,1)
    start=58;
    AVG=[];
    for i=1:3
        AVG(i,:)=ZS(idx_ZS,start:start+40);
        start=start+40;
    end
    STD=std(AVG,1,1);
    AVG=mean(AVG,1);
    AVG=AVG-min(AVG);
    j=1;
    for j=2:6
        for i=1:3
            temp(i,:)=ZS(idx_ZS,start:start+40);
            start=start+40;
        end
        temp_STD=std(temp,1,1);
        temp=mean(temp,1);
        temp=temp-min(temp);
        STD=[STD temp_STD];
        AVG=[AVG temp];
    end
    ZS_AVG(idx_ZS,:)=AVG;
    ZS_STD(idx_ZS,:)=STD;
end

AVG_Stimuli=zeros(6,size(ZS_AVG,2));
GCaMP6=[0,1.69644104899772,5.13796058542217,8.27886020152244,10.3756715204800,11.8173714529814,12.2425184714093,10.8571417354877,8.80831829681196,6.91339112244670,5.46959264663869,4.30868766622567,3.42533619066766,2.75378443486879,2.18017250852183,1.72816235135824,1.32732537295463,1.00684435500268,0.730210038304555,0.530242444093118,0.362253250339685,0.227668255288566,0.0869242416152502,0.000718266708050853,-0.0828334873368325]';
idxStart=1;
for i=1:6    
    AVG_Stimuli(i,(idxStart+(i-1)*40):(idxStart+(i-1)*40)+size(GCaMP6,1)-1)=GCaMP6;
end

ModelMultipower_AVG=[];
parfor i=1:length(ZS_AVG)
    %mdl=stepwiselm(Regressor',DataSet(i,:),'linear','Criterion','adjrsquared','Intercept',false,'Upper','interactions','Verbose',0);
    mdl=fitlm(AVG_Stimuli',ZS_AVG(i,:),'interactions');
    ModelMultipower_AVG(i).coef=mdl.Coefficients;
    %ModelResultsSeg_ZS(i).MSE=mdl.MSE;
    %ModelResultsSeg_ZS(i).Fitted=mdl.Fitted;
    ModelMultipower_AVG(i).rsquared=mdl.Rsquared.Adjusted;
end

options = statset('UseParallel',1); [idxKmeans_ZS_AVG_rsq Cmap_ZS_AVG_rsq]=kmeans(ZS_AVG(idx_rsq_AVG,:),10,'Options',options,'Distance','cityblock','Replicates',3,'MaxIter',1000,'Display','final');
options = statset('UseParallel',1); [idxKmeans_ZS_AVG_rsq2 Cmap_ZS_AVG_rsq2]=kmeans(ZS_AVG(idx_rsq_AVG,:),20,'Options',options,'Distance','cityblock','Replicates',3,'MaxIter',1000,'Display','final');
options = statset('UseParallel',1); [idxKmeans_ZS_AVG_rsq3 Cmap_ZS_AVG_rsq3]=kmeans(ZS_AVG(idx_rsq_AVG,:),5,'Options',options,'Distance','cityblock','Replicates',3,'MaxIter',1000,'Display','final');
[Model_ZS,GoodBetas_ZS]=Test_Regress(Cmap_ZS_AVG_rsq2,AVG_Stimuli,idxKmeans_ZS_AVG_rsq2,0.2);
[Model_ZS2,GoodBetas_ZS2]=Test_Regress(Cmap_ZS_AVG_rsq,AVG_Stimuli,idxKmeans_ZS_AVG_rsq,0.2);
[Model_ZS3,GoodBetas_ZS3]=Test_Regress(Cmap_ZS_AVG_rsq3,AVG_Stimuli,idxKmeans_ZS_AVG_rsq3,0.2);

GoodBetasSelect=GoodBetas_ZS([1 2 3 6 7 8 10 12 13 14 18 20]);

figure;
for i=1:20
    idx_temp=find(idxKmeans_ZS_AVG_rsq2==i);
    plot(mean(ZS_AVG(idx_rsq_AVG(idx_temp),:),1));
    pause
end

