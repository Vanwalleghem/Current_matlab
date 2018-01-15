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
clearvars GC C S F N name i GS Calcium Noise Fitness GN
ZS=zscore(GoodCalcium,1,2);

options = statset('UseParallel',1); [idxKmeans_ZS Cmap_ZS]=kmeans(ZS,30,'Options',options,'Replicates',3,'MaxIter',1000,'Display','final');
figure;imagesc(Cmap_ZS);

[coeff,score,~,~,explained,~] = pca(ZS);

spike=[0,1.69644104899772,10.3756715204800,12.2425184714093,10.8571417354877,8.80831829681196,6.91339112244670,5.46959264663869,4.30868766622567,3.42533619066766,2.75378443486879,2.18017250852183,1.72816235135824,1.32732537295463,1.00684435500268,0.730210038304555,0.530242444093118,0.362253250339685,0.0869242416152502,0.000718266708050853]';
framerate=2;
Lena_regressor=zeros(4,size(ZS,2));
for i = [68, 82, 115,150,157,189,200,218,232]
    idx=i*framerate;
    Lena_regressor(1,idx-1:idx-1+length(spike)-1)=spike';
end
for i = [75, 89, 108,164,171,182,207,225,239]
    idx=i*framerate;
    Lena_regressor(2,idx-1:idx-1+length(spike)-1)=spike';
end
for i = [64, 100, 122,157,174,178,189,196,210,214,228,232]
    idx=i*framerate;
    Lena_regressor(3,idx-1:idx-1+length(spike)-1)=spike';
end
for i = [60, 96, 104,150,167,178,185,196,200,214,218,242]
    idx=i*framerate;
    Lena_regressor(4,idx-1:idx-1+length(spike)-1)=spike';
end
x = linspace(1/framerate,size(ZS,2)/framerate,size(ZS,2));y = linspace(1,size(ZS,1),size(ZS,1));

LM_FMR=[];
parfor i=1:length(ZS)
    %mdl=stepwiselm(Regressor',DataSet(i,:),'linear','Criterion','adjrsquared','Intercept',false,'Upper','interactions','Verbose',0);
    mdl=fitlm(Lena_regressor',ZS(i,:),'interactions');
    LM_FMR(i).coef=mdl.Coefficients;
    %ModelResultsSeg_ZS(i).MSE=mdl.MSE;
    %ModelResultsSeg_ZS(i).Fitted=mdl.Fitted;
    LM_FMR(i).rsquared=mdl.Rsquared.Adjusted;
end
idx_rsq=find([LM_FMR.rsquared]>0.15);
ZS_rsq=ZS(idx_rsq,:);

Fighandle=figure;
set(Fighandle, 'Position', [100, 100, 1200, 900]);
imagesc(x,y,ZS_rsq(randperm(size(ZS_rsq,1)),:),[-0.5 5]);colormap hot;set(gca,'YTickLabel',[]);

[Model_ZS,GoodBetas_ZS]=Test_Regress(Cmap_ZS,Lena_regressor,idxKmeans_ZS,0.2);

figure;plot(Cmap_ZS(11,:));hold on;plot(Cmap_ZS(20,:));
for i = [68, 82, 115,150,157,189,200,218,232]
        i=i*framerate;
        coloring='r';duration=10;
        rectangle('EdgeColor',coloring,'FaceColor',coloring,'Position',[i -1 duration 0.25]);
end
for i = [75, 89, 108,164,171,182,207,225,239]
        i=i*framerate;
        coloring='b';duration=5;
        rectangle('EdgeColor',coloring,'FaceColor',coloring,'Position',[i -1 duration 0.25]);
end
for i = [64, 100, 122,157,174,178,189,196,210,214,228,232]
        i=i*framerate;
        coloring='g';duration=1;
        rectangle('EdgeColor',coloring,'FaceColor',coloring,'Position',[i -0.8 duration 0.25]);
end
for i = [60, 96, 104,150,167,178,185,196,200,214,218,242]
        i=i*framerate;
        coloring='m';duration=1;
        rectangle('EdgeColor',coloring,'FaceColor',coloring,'Position',[i -0.6 duration 0.25]);
end

options = statset('UseParallel',1); [idxKmeans_ZS_rsq Cmap_ZS_rsq]=kmeans(ZS_rsq,20,'Options',options,'Distance','cityblock','Replicates',3,'MaxIter',1000,'Display','final');
figure;imagesc(Cmap_ZS_rsq);
[Model_ZS_rsq,GoodBetas_ZS_rsq]=Test_Regress(Cmap_ZS_rsq,Lena_regressor,idxKmeans_ZS_rsq,0.3);

Tail_fish2=zscore(Tail_fish2);
Tail_fish2b=Tail_fish2;Tail_fish2b(Tail_fish2>-2.5)=0;
 [pks,locs] = findpeaks(Tail_fish2b,'MinPeakDistance',10);
 Tail_fish2=round((locs-30)/5);
 
 Lena_regressor(5,:)=Lena_regressor(4,:)*0;
 for idx = Tail_fish2'
    Lena_regressor(5,idx-1:idx-1+length(spike)-1)=spike';
 end
figure; plot(Lena_regressor(end,:)/10); hold on;plot(Cmap_ZS(20,:));

Tail_fish4=zscore(Tail_fish4);
Tail_fish4b=Tail_fish4-Background;
Tail_fish4b(Tail_fish4b>-0.5)=0;
 [pks,locs] = findpeaks(Tail_fish4b,'MinPeakDistance',10);
 locs=round((locs-30)/5);
 
 Lena_regressor(6,:)=Lena_regressor(4,:)*0;
 for idx = locs'
    Lena_regressor(6,idx-1:idx-1+length(spike)-1)=spike';
 end
figure; plot(Lena_regressor(end,:)/10); hold on;plot(Cmap_ZS(20,:));

Tail_fish3=zscore(Tail_fish3);
Tail_fish3b=Tail_fish3-zscore(Background2);
Tail_fish3b(Tail_fish3b>-1)=0;
 [pks,locs] = findpeaks(Tail_fish3b,'MinPeakDistance',10);
 locs=round((locs-30)/5);
 
 Lena_regressor(7,:)=Lena_regressor(4,:)*0;
 for idx = locs'
    Lena_regressor(7,idx-1:idx-1+length(spike)-1)=spike';
 end
figure; plot(Lena_regressor(end,:)/10); hold on;plot(Cmap_ZS(20,:));
 