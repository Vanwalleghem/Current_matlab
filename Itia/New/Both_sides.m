MatFiles=dir('*analysis_matlab.mat');
name=strcat(MatFiles(1).name);
Calcium=load(name, 'DenoisedTraces');
Calcium=Calcium.DenoisedTraces;
MatFiles(1).number=size(Calcium,1);
%Spikes=load(name, 'Spikes');
%Spikes=Spikes.Spikes;
Noise=load(name, 'Noise');
Noise=Noise.Noise;
Fitness=load(name, 'idx_components');
Fitness=Fitness.idx_components+1;
GoodCalcium=Calcium(Fitness,:);
%GoodSpikes=Spikes(Fitness,:);
GoodNoise=Noise(Fitness,:);
MatFiles(1).GoodNumber=length(Fitness);
%MatFiles(1).GC=GoodCalcium;
for i = 2:length(MatFiles)
    name=strcat(MatFiles(i).name);
    C=load(name, 'DenoisedTraces');
    C=C.DenoisedTraces;
%     if i==3
%         C=[C(:,1) C(:,1) C(:,1:58)];
%     end
    %S=load(name, 'Spikes');
    %S=S.Spikes;
    N=load(name, 'Noise');
    N=N.Noise;
    F=load(name, 'idx_components');
    F=F.idx_components+1;
    GC=C(F,:);
    GN=N(F,:);
    %GS=S(F,:);
    Noise=vertcat(Noise,N);
    Calcium=vertcat(Calcium,C);
    %Spikes=vertcat(Spikes,S);
    Fitness=horzcat(Fitness,F);
    GoodCalcium=vertcat(GoodCalcium,GC);
    GoodNoise=vertcat(GoodNoise,GN);
    %GoodSpikes=vertcat(GoodSpikes,GS);
    MatFiles(i).number=size(Calcium,1);
    MatFiles(i).GoodNumber=MatFiles(i-1).GoodNumber+length(F);
    %MatFiles(i).GC=GC;
end
clearvars GC C S F N name i GS Calcium Noise Fitness

MatFiles=dir('*analysis_matlab.mat');
name=strcat(MatFiles(1).name);
Baseline=load(name, 'Baseline');
Baseline=Baseline.Baseline;
Baseline=cell2mat(Baseline);
Fitness=load(name, 'idx_components');
Fitness=Fitness.idx_components+1;
Baseline=Baseline(Fitness);
for i = 2:length(MatFiles)
    name=strcat(MatFiles(i).name);
    B=load(name, 'Baseline');
    B=B.Baseline;
    B=cell2mat(B);
    Fitness=load(name, 'idx_components');
    Fitness=Fitness.idx_components+1;
    B=B(Fitness);
    Baseline=horzcat(Baseline,B);    
    %MatFiles(i).Baseline=Baseline;    
end
clearvars GC C S F N name i GS Calcium Noise Fitness

ZS=zscore(GoodCalcium,1,2);
ZS2=zscore(bsxfun(@plus,GoodCalcium+GoodNoise,Baseline'),1,2);
x = linspace(0.2,size(ZS,2)/5,size(ZS,2));y = linspace(1,size(ZS,1),size(ZS,1));
Fighandle=figure;
set(Fighandle, 'Position', [100, 100, 1200, 900]);
imagesc(x,y,ZS(randperm(size(ZS,1)),:),[-0.5 5]);colormap hot;set(gca,'YTickLabel',[]);

Numbers=[1 [MatFiles.GoodNumber]];
counter=1;
idx_Plane=nan(length(ZS),1);
idx_Position=nan(length(ZS),1);
idx_Fish=nan(length(ZS),1);
name=strcat(MatFiles(1).name);
for i=1:length(MatFiles)	
    name=strcat(MatFiles(i).name);
    [Plane,~]=regexp(name,'_(\d+)_','tokens','match');Plane=str2num(Plane{1}{1});
    if strfind(name,'ILO')
        Fish=2;
    elseif strfind(name,'IRO')
        Fish=3;
    elseif strfind(name,'ERO')
        Fish=1;
    end
    idx_Plane(Numbers(i):Numbers(i+1))=Plane;
    idx_Position(Numbers(i):Numbers(i+1))=Fish;
    [Fish,~]=regexp(name,'Fish2017(\d+)_','tokens','match');Fish=str2num(Fish{1}{1});
    idx_Fish(Numbers(i):Numbers(i+1))=Fish;
end
clearvars i Fish Plane name counter

options = statset('UseParallel',1); [idxKmeans_ZS Cmap_ZS]=kmeans(ZS,20,'Options',options,'MaxIter',1000,'Display','final');
figure;imagesc(Cmap_ZS,[-0.5 3]);colormap hot


Fighandle=figure;
set(Fighandle, 'Position', [100, 100, 1200, 900]);
imagesc(x,y,ZS(randperm(size(ZS_ILO,1)),:),[-0.5 5]);colormap hot;set(gca,'YTickLabel',[]);
figure;histogram(idx_Fish(idx_ILO));

Fighandle=figure;
set(Fighandle, 'Position', [100, 100, 1200, 900]);counter=1;
for i=1:3
    idx_temp=find(idx_Position==i);
    subplot(3,2,counter);imagesc(ZS(idx_temp,:),[-0.5 4]);colormap hot;
    subplot(3,2,counter+1);plot(mean(ZS(idx_temp,:),1));
    counter=counter+2;
end    

Clusters={};
options = statset('UseParallel',1);
for i=1:3
    idx_temp=find(idx_Position==i);
    [Clusters{1,i}, Clusters{2,i}]=kmeans(ZS(idx_temp,:),20,'Options',options,'MaxIter',1000,'Display','final');
end    
Fighandle=figure;
set(Fighandle, 'Position', [100, 100, 1200, 900]);counter=1;
for i=1:3
    subplot(3,1,i);imagesc(Clusters{2,i},[-0.5 4]);colormap hot;
end
GoodBetas_pos={};
GoodBetas_pos{1}=[1 5 7 20];
GoodBetas_pos{2}=[13 15];
GoodBetas_pos{3}=[2 7 10 17];
Fighandle=figure;
set(Fighandle, 'Position', [100, 100, 1200, 900]);counter=1;
for i=1:3
    Cmap=Clusters{2,i};
    subplot(3,1,i);plot(Cmap(GoodBetas_pos{i},:)');axis([0 400 -0.5 3]);
end

Fighandle=figure;
set(Fighandle, 'Position', [100, 100, 1200, 900]);counter=1;
for i=1:3
    Cmap=Clusters{2,i};
    plot(mean(Cmap(GoodBetas_pos{i},:),1));axis([0 400 -0.5 3]);hold  on;
end
legend('ERO-ELO','ILO-ERO','IRO-ELO');


counter=1;
x = linspace(0.5,size(ZS,2)/2,size(ZS,2));
rows=length(GoodBetas);
counter=1;
Fighandle=figure;
set(Fighandle, 'Position', [100, 100, 1300, 900]);
for i=GoodBetas
    idx_temp=find(idxKmeans_final==i);
    subplot(rows,5,counter);plot(mean(ZS(idx_temp,:),1));
    subplot(rows,5,counter+1);imagesc(ZS(idx_temp,:),[-0.5 4]);
    subplot(rows,5,counter+2);histogram(idx_Plane(idx_temp));
    subplot(rows,5,counter+3);histogram(idx_Position(idx_temp),[0.5:1:3.5]);h = gca;h.XTick=[1 2 3];h.XTickLabel={'ERO-ELO','ERO-ILO','IRO-ELO'};
    subplot(rows,5,counter+4);histogram(idx_Fish(idx_temp));%h = gca;h.XTickLabel={'1','2'};
    counter=counter+5;
end


Stimuli=zeros(6,size(ZS,2));
start=40;
spike=[0,1.69644104899772,5.13796058542217,8.27886020152244,10.3756715204800,12.2425184714093,10.8571417354877,8.80831829681196,6.91339112244670,5.46959264663869,4.30868766622567,3.42533619066766,2.75378443486879,2.18017250852183,1.72816235135824,1.32732537295463,1.00684435500268,0.730210038304555,0.530242444093118,0.362253250339685,0.0869242416152502,0.000718266708050853]';
for i = 1:3
    for j = 0:2
        Stimuli(i,start+j*120:start+j*120+length(spike)-1)=spike';
    end
    start=start+40;
end
start=44;
for i = 1:3
    for j = 0:2
        Stimuli(i+3,start+j*120:start+j*120+length(spike)-1)=spike';
    end
    start=start+40;
end

ModelBothSides=[];
parfor i=1:length(ZS)
    %mdl=stepwiselm(Regressor',DataSet(i,:),'linear','Criterion','adjrsquared','Intercept',false,'Upper','interactions','Verbose',0);
    mdl=fitlm(Stimuli',ZS(i,:),'interactions');
    ModelBothSides(i).coef=mdl.Coefficients;
    %ModelResultsSeg_ZS(i).MSE=mdl.MSE;
    %ModelResultsSeg_ZS(i).Fitted=mdl.Fitted;
    ModelBothSides(i).rsquared=mdl.Rsquared.Adjusted;
end
idx_rsq=find([ModelBothSides.rsquared]>0.15);
ZS_rsq=ZS(idx_rsq,:);
Fighandle=figure;
set(Fighandle, 'Position', [100, 100, 1200, 900]);
imagesc(x,y,ZS_rsq(randperm(size(ZS_rsq,1)),:),[-0.5 5]);colormap hot;set(gca,'YTickLabel',[]);

[idxKmeans_ZS_rsq Cmap_ZS_rsq]=kmeans(ZS_rsq,20,'Options',options,'MaxIter',1000,'Display','final');
Fighandle=figure;
set(Fighandle, 'Position', [100, 100, 1200, 900]);
imagesc(Cmap_ZS_rsq,[-0.5 4]);colormap hot;

figure;
for i=1:20
    plot(Cmap_ZS_rsq(i,:));pause;
end

coefficients={};
for idx=1:length(ModelBothSides)
    coef=[ModelBothSides(idx).coef];
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

coef_stats=zeros(2,size(coefficients,2))
coef_stats(1,:)=mean(coefficients,1);
coef_stats(2,:)=std(coefficients,1,1);

idx_left=intersect(find(max(ZS(:,40:50),[],2)>0.5),find(max(ZS(:,80:90),[],2)<-0.2));
figure;imagesc(ZS(idx_left,:),[-2 3]);colormap hot;
figure;plot(mean(ZS(idx_left,:),1));

ZS_AVG=zeros(size(ZS,1),123);
ZS_STD=zeros(size(ZS,1),123);
parfor idx_ZS=1:size(ZS,1)
    start=30;
    AVG=[];
    for i=1:3
        AVG(i,:)=ZS(idx_ZS,start+(i-1)*120:start+40+(i-1)*120);
        %start=start+40;
    end
    STD=std(AVG,1,1);
    AVG=mean(AVG,1);
    AVG=AVG-min(AVG);
    j=1;
    for j=2:3
        start=start+40;
        for i=1:3
            temp(i,:)=ZS(idx_ZS,start+(i-1)*120:start+40+(i-1)*120);            
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

[idxKmeans_ZS_AVG Cmap_ZS_AVG]=kmeans(ZS_AVG,20,'Options',options,'MaxIter',1000,'Display','final');
Fighandle=figure;
set(Fighandle, 'Position', [100, 100, 1200, 900]);
imagesc(Cmap_ZS_AVG,[0.5 3]);colormap gray;

counter=1;
x = linspace(0.5,size(ZS,2)/2,size(ZS,2));
rows=length(GoodBetas_ZS);
counter=1;
Fighandle=figure;
set(Fighandle, 'Position', [100, 100, 300, 1300]);
colors = distinguishable_colors(length(GoodBetas_ZS),[1 1 1; 0 0 0]);
colors = colors*256;
for i=GoodBetas_ZS
    idx=find(idxKmeans_ZS_rsq==i);
    subplot(rows,1,counter);plot(mean(ZS_AVG(idx_rsq(idx),:),1),'color',colors(counter,:)/256);xlim([0 120]);
    counter=counter+1;
end    

idxKmeans_final=size(idxKmeans_ZS);
idxKmeans_final(idx_rsq)=idxKmeans_ZS_rsq;

[idxKmeans_ZS_rsq Cmap_ZS_rsq]=kmeans(ZS_rsq,10,'Options',options,'MaxIter',1000,'Display','final');
Fighandle=figure;
set(Fighandle, 'Position', [100, 100, 1200, 900]);
imagesc(Cmap_ZS_rsq,[-0.5 4]);colormap hot;

[idxKmeans_ZS_rsq2 Cmap_ZS_rsq2]=kmeans(ZS2(idx_rsq,:),10,'Options',options,'MaxIter',1000,'Display','final');
Fighandle=figure;
set(Fighandle, 'Position', [100, 100, 1200, 900]);
imagesc(Cmap_ZS_rsq2,[-0.5 4]);colormap hot;

figure;
for i=1:size(Cmap_ZS_rsq,1)
    plot(Cmap_ZS_rsq(i,:));pause;
end