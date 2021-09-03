MatFiles=dir('*Vestib*analysis_matlab.mat');
name=strcat(MatFiles(1).name);
Calcium=load(name, 'DenoisedTraces');
Calcium=Calcium.DenoisedTraces;
MatFiles(1).number=size(Calcium,1);
Noise=load(name, 'Noise');
Noise=Noise.Noise;
Fitness=load(name, 'idx_components');
Fitness=Fitness.idx_components+1;
GoodCalcium=Calcium(Fitness,:);
GoodNoise=Noise(Fitness,:);
MatFiles(1).GoodNumber=length(Fitness);
for i = 2:length(MatFiles)
    name=strcat(MatFiles(i).name);
    C=load(name, 'DenoisedTraces');
    C=C.DenoisedTraces;
    N=load(name, 'Noise');
    N=N.Noise;
    F=load(name, 'idx_components');
    F=F.idx_components+1;
    GC=C(F,:);
    GN=N(F,:);
    Noise=vertcat(Noise,N);
    Calcium=vertcat(Calcium,C);
    Fitness=horzcat(Fitness,F);
    GoodCalcium=vertcat(GoodCalcium,GC);
    GoodNoise=vertcat(GoodNoise,GN);
    MatFiles(i).number=size(Calcium,1);
    MatFiles(i).GoodNumber=MatFiles(i-1).GoodNumber+length(F);
end
clearvars GC C S F N name i GS Calcium Noise Fitness
ZS=zscore(GoodCalcium+GoodNoise,1,2);
options = statset('UseParallel',1);
[idxKmeans Cmap]=kmeans(ZS,10,'Options',options,'Replicates',3,'MaxIter',1000,'Display','final');

%Vestibular: Stimuli onset at timepoints (32 stimuli *3 repetitions)
[38:20:658,708:1328,1378:20:1998]

%Audio: stimuli onset at timepoints (16 stimuli *3 repetitions) ; [48:25:1223]
GCaMP=[0,1.69644104899772,5.13796058542217,8.27886020152244,10.3756715204800,11.8173714529814,12.2425184714093,10.8571417354877,8.80831829681196,6.91339112244670,5.46959264663869,4.30868766622567,3.42533619066766,2.75378443486879,2.18017250852183,1.72816235135824,1.32732537295463,1.00684435500268,0.730210038304555,0.530242444093118,0.362253250339685,0.227668255288566,0.0869242416152502,0.000718266708050853,-0.0828334873368325]';
Stimuli_Audio=zeros(32,size(ZS,2));
counter=1;start=37;
for i=1:32    
    Stimuli_Audio(counter,start+(i-1)*20:start+(i-1)*20+length(GCaMP)-1)=GCaMP;
    counter=counter+1;
end

counter=1;start=657;
for i=1:32    
    Stimuli_Audio(counter,start+(i-1)*20:start+(i-1)*20+length(GCaMP)-1)=GCaMP;
    counter=counter+1;
end

counter=1;start=1377;
for i=1:32    
    Stimuli_Audio(counter,start+(i-1)*20:start+(i-1)*20+length(GCaMP)-1)=GCaMP;
    counter=counter+1;
end

Stimuli_Audio=Stimuli_Audio(:,1:size(ZS,2));

ModelLinReg=[];
parfor i=1:length(ZS)    
    mdl=fitlm(Stimuli_Audio',ZS(i,:));%,'interactions');
    ModelLinReg(i).coef=mdl.Coefficients;    
    ModelLinReg(i).rsquared=mdl.Rsquared.Adjusted;
end

figure;histogram([ModelLinReg.rsquared]);
idx_rsq=find([ModelLinReg.rsquared]>0.1);
[idxKmeans Cmap]=kmeans(ZS(idx_rsq,:),10,'Options',options,'Replicates',3,'MaxIter',1000,'Display','final');
[idxKmeans_20 Cmap_20]=kmeans(ZS(idx_rsq,:),20,'Options',options,'Replicates',3,'MaxIter',1000,'Display','final');

coefficients={};
for idx=1:length(ModelLinReg)
    coef=[ModelLinReg(idx).coef];
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

[idxKmeans_coef Cmap_coef]=kmeans(coefficients(idx_rsq,:),10,'Options',options,'Replicates',3,'MaxIter',1000,'Display','final');
[idxKmeans_coef_20 Cmap_coef_20]=kmeans(coefficients(idx_rsq,:),20,'Options',options,'Replicates',3,'MaxIter',1000,'Display','final');

Cluster_mean=zeros(max(idxKmeans_coef(:)), size(ZS,2));
for i=1:max(idxKmeans_coef)
    Cluster_mean(i,:)=mean(ZS(idx_rsq(find(idxKmeans_coef==i)),:),1);
end

Cluster_mean_20=zeros(max(idxKmeans_coef_20(:)), size(ZS,2));
for i=1:max(idxKmeans_coef_20)
    Cluster_mean_20(i,:)=mean(ZS(idx_rsq(find(idxKmeans_coef_20==i)),:),1);
end

figure;
for i=1:max(idxKmeans_coef_20)
    plot(mean(ZS(idx_rsq(find(idxKmeans_coef_20==i)),:),1));title(num2str(i));
    hold on;plot((Stimuli_Audio'/20)-1);hold off;
    pause
end
