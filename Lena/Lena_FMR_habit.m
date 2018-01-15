MatFiles=dir('*analysis_matlab.mat');
name=strcat(MatFiles(1).name);
Calcium=load(name, 'DenoisedTraces');
Calcium=Calcium.DenoisedTraces;
MatFiles(1).number=size(Calcium,1);
%Spikes=load(name, 'Spikes');
%Spikes=Spikes.Spikes;
%Noise=load(name, 'Noise');
%Noise=Noise.Noise;
Fitness=load(name, 'idx_components');
Fitness=Fitness.idx_components+1;
GoodCalcium=Calcium(Fitness,:);
%GoodSpikes=Spikes(Fitness,:);
%GoodNoise=Noise(Fitness,:);
MatFiles(1).GoodNumber=length(Fitness);
%MatFiles(1).GC=GoodCalcium;
progressbar;
for i = 2:length(MatFiles)
    name=strcat(MatFiles(i).name);
    C=load(name, 'DenoisedTraces');
    C=C.DenoisedTraces;
%     if i==3
%         C=[C(:,1) C(:,1) C(:,1:58)];
%     end
    %S=load(name, 'Spikes');
    %S=S.Spikes;
%    N=load(name, 'Noise');
   % N=N.Noise;
    F=load(name, 'idx_components');
    F=F.idx_components+1;
    GC=C(F,:);
    %GN=N(F,:);
    %GS=S(F,:);
    %Noise=vertcat(Noise,N);
    %Calcium=vertcat(Calcium,C);
    %Spikes=vertcat(Spikes,S);
    Fitness=horzcat(Fitness,F);
    GoodCalcium=vertcat(GoodCalcium,GC);
    %GoodNoise=vertcat(GoodNoise,GN);
    %GoodSpikes=vertcat(GoodSpikes,GS);
    MatFiles(i).number=MatFiles(i-1).number+size(C,1);
    MatFiles(i).GoodNumber=MatFiles(i-1).GoodNumber+length(F);
    %MatFiles(i).GC=GC;
    progressbar(i/length(MatFiles));
end
clearvars GC C S F N name i GS Calcium Noise Fitness GN
save('FMR_GoodCalcium.mat','GoodCalcium','-v7.3');
save('FMR_GoodNoise.mat','GoodNoise','-v7.3');
clearvars GoodCalcium GoodNoise
ZS=zscore(GoodCalcium,1,2);

Numbers=[1 [MatFiles.GoodNumber]];
counter=1;
idx_Plane=nan(length(ZS),1);
idx_Fish=nan(length(ZS),1);
name=strcat(MatFiles(1).name);
for i=1:length(MatFiles)	
    name=strcat(MatFiles(i).name);
    if strfind(name,'clutch2')
        [Fish,~]=regexp(name,'clutch2_fish(\d)_6','tokens','match');Fish=str2num(Fish{1}{1});          
        Fish=Fish+4;
    else
         [Fish,~]=regexp(name,'s_(\d)_6','tokens','match');Fish=str2num(Fish{1}{1});   
    end
    idx_Fish(Numbers(i):Numbers(i+1))=Fish;
    [Plane,~]=regexp(name,'Slice(\d+)_','tokens','match');Plane=str2num(Plane{1}{1});
    idx_Plane(Numbers(i):Numbers(i+1))=Plane;
end
clearvars i Fish Plane name counter


options = statset('UseParallel',0); [idxKmeans_ZS Cmap_ZS]=kmeans(ZS,20,'Options',options,'Replicates',3,'MaxIter',1000,'Display','final');
figure;imagesc(Cmap_ZS);

[coeff,score,~,~,explained,~] = pca(ZS);

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


framerate=2;
Tail_fish1=zscore(Tail_fish1);
Tail_fish1b=Tail_fish1-zscore(Background1);
Tail_fish1b(Tail_fish1b<0.5)=0;
 [pks,locs] = findpeaks(Tail_fish1b,'MinPeakDistance',10);
 locs=round(locs/(10/framerate));

 [~,looms]=findpeaks(-zscore(Background1),'MinPeakHeight',1,'MinPeakDistance',120);
looms(looms<1000)=[];
looms=round(looms/(10/framerate));
locs=locs-looms(1)+315;
locs(locs<1)=[];
looms=looms-looms(1)+315;
 
spike=[0,1.69644104899772,10.3756715204800,12.2425184714093,10.8571417354877,8.80831829681196,6.91339112244670,5.46959264663869,4.30868766622567,3.42533619066766,2.75378443486879,2.18017250852183,1.72816235135824,1.32732537295463,1.00684435500268,0.730210038304555,0.530242444093118,0.362253250339685,0.0869242416152502,0.000718266708050853]';
framerate=2;
Lena_regressor=zeros(2,size(ZS,2));
for idx = looms'    
    Lena_regressor(1,idx-1:idx-1+length(spike)-1)=spike';
end 
 
 for idx = locs'
    Lena_regressor(2,idx-1:idx-1+length(spike)-1)=spike';
 end

 LM_FMR=[];
parfor i=1:length(ZS)
    %mdl=stepwiselm(Regressor',DataSet(i,:),'linear','Criterion','adjrsquared','Intercept',false,'Upper','interactions','Verbose',0);
    mdl=fitlm(Lena_regressor(1,:)',ZS(i,:),'interactions');
    LM_FMR(i).coef=mdl.Coefficients;
    %ModelResultsSeg_ZS(i).MSE=mdl.MSE;
    %ModelResultsSeg_ZS(i).Fitted=mdl.Fitted;
    LM_FMR(i).rsquared=mdl.Rsquared.Adjusted;
end
idx_rsq=find([LM_FMR.rsquared]>0.05);
ZS_rsq=ZS(idx_rsq,:);

Fighandle=figure;
set(Fighandle, 'Position', [100, 100, 1200, 900]);
x = linspace(1/framerate,size(ZS,2)/framerate,size(ZS,2));y = linspace(1,size(ZS,1),size(ZS,1));
imagesc(x,y,ZS_rsq(randperm(size(ZS_rsq,1)),:),[-0.5 5]);colormap hot;set(gca,'YTickLabel',[]);

options = statset('UseParallel',1); [idxKmeans_ZS_rsq Cmap_ZS_rsq]=kmeans(ZS_rsq,15,'Options',options,'Distance','cityblock','Replicates',3,'MaxIter',1000,'Display','final');

[Model_ZS,GoodBetas_ZS]=Test_Regress(Cmap_ZS_rsq,Lena_regressor,idxKmeans_ZS_rsq,0.15);
GoodBetas=GoodBetas_ZS([1 2 3]);

idxKmeans_final=zeros(size(ZS,1),1);
idxKmeans_final(idx_rsq)=idxKmeans_ZS_rsq;

Fighandle=figure;
set(Fighandle, 'Position', [100, 100, 1200, 900]);counter=1;yplot=length(GoodBetas);
for i=GoodBetas
    idx_temp=find(idxKmeans_final==i);
    subplot(yplot,3,counter);plot(mean(ZS(idx_temp,:),1));ylim([-1 7]);title(num2str(length(idx_temp)))
    subplot(yplot,3,counter+1);histogram(idx_Fish(idx_temp));
    subplot(yplot,3,counter+2);histogram(idx_Plane(idx_temp));
    counter=counter+yplot;
end


