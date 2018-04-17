ZS=zscore(GoodCalcium+GoodNoise,1,2);
Stimuli=zeros(3,size(ZS,2));
start=43;
spike=[0,1.69644104899772,5.13796058542217,8.27886020152244,10.3756715204800,12.2425184714093,10.8571417354877,8.80831829681196,6.91339112244670,5.46959264663869,4.30868766622567,3.42533619066766,2.75378443486879,2.18017250852183,1.72816235135824,1.32732537295463,1.00684435500268,0.730210038304555,0.530242444093118,0.362253250339685,0.0869242416152502,0.000718266708050853]';
for i = 1:3
    for j = 0:2
        Stimuli(i,start+j*120:start+j*120+length(spike)-1)=spike';
    end
    start=start+20;
end
clearvars GoodCalcium GoodNoise

ModelBothSides=[];
parfor i=1:length(ZS)    
    mdl=fitlm(Stimuli',ZS(i,:));%,'interactions');
    ModelBothSides(i).coef=mdl.Coefficients;        
    ModelBothSides(i).rsquared=mdl.Rsquared.Adjusted;
end
idx_rsq=find([ModelBothSides.rsquared]>0.1);
figure;histogram([ModelBothSides.rsquared]);
ZS_rsq=ZS(idx_rsq,:);

options = statset('UseParallel',1); [idxKmeans_ZS_rsq Cmap_ZS_rsq]=kmeans(ZS_rsq,20,'Options',options,'Distance','cityblock','Replicates',3,'MaxIter',1000,'Display','final');
[ModelCmap,GoodBetas]=Test_Regress(Cmap_ZS_rsq,Stimuli,idxKmeans_ZS_rsq,0.3);

%Building the AVG across 3 presentations
Firststart=30;interstimulus=120;Nb_stimuli=3;    
ZS_AVG=zeros(size(ZS_rsq,1),Nb_stimuli*41);
parfor idx_ZS=1:size(ZS_rsq,1)
    start=Firststart;
    AVG=[];
    for i=1:3
        AVG(i,:)=ZS_rsq(idx_ZS,start:start+40);
        start=start+interstimulus;
    end
    AVG=mean(AVG,1);
    AVG=AVG-mean(AVG(1:5));%-min(AVG);
    temp=[];
    for j=2:3
        start=Firststart+40*(j-1);
        for i=1:3
            temp(i,:)=ZS_rsq(idx_ZS,start:start+40);
            start=start+interstimulus;
        end
        temp=mean(temp,1);
        temp=temp-mean(temp(1:5));%-min(temp);
        AVG=[AVG temp];
    end
    ZS_AVG(idx_ZS,:)=AVG;
end


Stimuli_AVG=zeros(3,123);
GCaMP6=[0,1.69644104899772,5.13796058542217,8.27886020152244,10.3756715204800,11.8173714529814,12.2425184714093,10.8571417354877,8.80831829681196,6.91339112244670,5.46959264663869,4.30868766622567,3.42533619066766,2.75378443486879,2.18017250852183,1.72816235135824,1.32732537295463,1.00684435500268,0.730210038304555,0.530242444093118,0.362253250339685,0.227668255288566,0.0869242416152502,0.000718266708050853,-0.0828334873368325]';
idxStart=10;
for i=1:3    
    Stimuli_AVG(i,(idxStart+(i-1)*41):(idxStart+(i-1)*41)+size(GCaMP6,1)-1)=GCaMP6;
end

options = statset('UseParallel',1); [idxKmeans_ZS_AVG Cmap_ZS_AVG]=kmeans(ZS_AVG,20,'Options',options,'Distance','cityblock','Replicates',3,'MaxIter',1000,'Display','final');
[ModelCmap_AVG,GoodBetas_AVG]=Test_Regress(Cmap_ZS_AVG,Stimuli_AVG,idxKmeans_ZS_AVG,0.3);

figure;
counter=1;xplot=floor(sqrt(length(GoodBetas_AVG)));yplot=ceil(length(GoodBetas_AVG)/xplot);
for i=GoodBetas_AVG
    idx_temp=find(idxKmeans_ZS_AVG==i);
    subplot(xplot,yplot,counter);plot(mean(ZS_rsq(idx_temp,:),1));ylim([-3 3]);
    counter=counter+1;
end