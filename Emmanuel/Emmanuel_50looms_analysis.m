MatFiles=dir('*50looms*matlab*.mat');
name=strcat(MatFiles(1).name);
Calcium=load(name, 'DenoisedTraces');
Calcium=Calcium.DenoisedTraces;
MatFiles(1).number=size(Calcium,1);
Spikes=load(name, 'Spikes');
Spikes=Spikes.Spikes;
Noise=load(name, 'Noise');
Noise=Noise.Noise;
%DF=load(name, 'dFonF');
%DF=DF.dFonF;
Fitness=load(name, 'idx_components');
Fitness=Fitness.idx_components+1;
GoodCalcium=Calcium(Fitness,:);
GoodSpikes=Spikes(Fitness,:);
GoodNoise=Noise(Fitness,:);
%GoodDF=DF(Fitness,:);
MatFiles(1).GoodNumber=length(Fitness);
for i = 2:length(MatFiles)
name=strcat(MatFiles(i).name);
C=load(name, 'DenoisedTraces');
C=C.DenoisedTraces;
%     if i==3
%         C=[C(:,1) C(:,1) C(:,1:58)];
%     end
S=load(name, 'Spikes');
S=S.Spikes;
N=load(name, 'Noise');
N=N.Noise;
F=load(name, 'idx_components');
F=F.idx_components+1;
%D=load(name, 'dFonF');
%D=D.dFonF;
GC=C(F,:);
GS=S(F,:);
%GD=D(F,:);
Noise=vertcat(Noise,N);
GN=N(F,:);
Calcium=vertcat(Calcium,C);
%DF=vertcat(DF,D);
Spikes=vertcat(Spikes,S);
Fitness=horzcat(Fitness,F);
GoodCalcium=vertcat(GoodCalcium,GC);
GoodNoise=vertcat(GoodNoise,GN);
%GoodDF=vertcat(GoodDF,GD);
GoodSpikes=vertcat(GoodSpikes,GS);
MatFiles(i).number=size(Calcium,1);
MatFiles(i).GoodNumber=MatFiles(i-1).GoodNumber+length(F);
end
clearvars GC C S F N name i GS GN N;
ZS=zscore(GoodCalcium,1,2);
options = statset('UseParallel',1); [idxKmeans_ZS Cmap_ZS]=kmeans(ZS,50,'Options',options,'Distance','cityblock','Replicates',5,'MaxIter',1000,'Display','final');

Stimuli=zeros(1,size(GoodCalcium,2));
GCaMP6=[-0.104392135015146,1.69644104899772,5.13796058542217,8.27886020152244,10.3756715204800,11.8173714529814,12.2425184714093,10.8571417354877,8.80831829681196,6.91339112244670,5.46959264663869,4.30868766622567,3.42533619066766,2.75378443486879,2.18017250852183,1.72816235135824,1.32732537295463,1.00684435500268,0.730210038304555,0.530242444093118,0.362253250339685,0.227668255288566,0.0869242416152502,0.000718266708050853,-0.0828334873368325]';
GCaMP6=interp(GCaMP6,2);
idx=105;
for i=1:20
    Stimuli(1,idx:idx+size(GCaMP6,1)-1)=GCaMP6';
    idx=idx+200;
end

Stimuli2=zeros(20,size(GoodCalcium,2));
GCaMP6=[-0.104392135015146,1.69644104899772,5.13796058542217,8.27886020152244,10.3756715204800,11.8173714529814,12.2425184714093,10.8571417354877,8.80831829681196,6.91339112244670,5.46959264663869,4.30868766622567,3.42533619066766,2.75378443486879,2.18017250852183,1.72816235135824,1.32732537295463,1.00684435500268,0.730210038304555,0.530242444093118,0.362253250339685,0.227668255288566,0.0869242416152502,0.000718266708050853,-0.0828334873368325]';
GCaMP6=interp(GCaMP6,2);
idx=105;
for i=1:20
    Stimuli2(i,idx:idx+size(GCaMP6,1)-1)=GCaMP6';
    idx=idx+200;
end