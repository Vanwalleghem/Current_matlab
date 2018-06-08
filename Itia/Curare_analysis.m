MatFiles=dir('*analysis_matlab.mat');
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
    %Noise=vertcat(Noise,N);
    %Calcium=vertcat(Calcium,C);
    Fitness=horzcat(Fitness,F);
    GoodCalcium=vertcat(GoodCalcium,GC);
    GoodNoise=vertcat(GoodNoise,GN);    
    %MatFiles(i).number=size(Calcium,1);
    MatFiles(i).GoodNumber=MatFiles(i-1).GoodNumber+length(F);    
end
clearvars GC C S F N name i GS Calcium Noise Fitness

ZS2=zscore(GoodCalcium+GoodNoise,1,2);

Numbers=[1 [MatFiles.GoodNumber]];
counter=1;
idx_Plane=nan(length(ZS2),1);
idx_Position=nan(length(ZS2),1);
idx_Fish=nan(length(ZS2),1);
name=strcat(MatFiles(1).name);
for i=1:length(MatFiles)	
    name=strcat(MatFiles(i).name);
    [Plane,~]=regexp(name,'Slice(\d+)_','tokens','match');Plane=str2num(Plane{1}{1});
    if strfind(name,'control')
        Fish=1;
    elseif strfind(name,'curare')
        Fish=2;  
    end
    idx_Plane(Numbers(i):Numbers(i+1))=Plane;
    idx_Position(Numbers(i):Numbers(i+1))=Fish;
    [Fish,~]=regexp(name,'Fish2018(\d+)_fish(\d)','tokens');Fish=str2num(strcat(Fish{1}{1},Fish{1}{2}));
    idx_Fish(Numbers(i):Numbers(i+1))=Fish;
end
clearvars i Fish Plane name counter

options = statset('UseParallel',1); [idxKmeans_ZS Cmap_ZS]=kmeans(ZS2,30,'Options',options,'Distance','cityblock','Replicates',5,'MaxIter',1000,'Display','final');

idxStart=120;
interval=40;
repetitions=10;
Stimuli=zeros(2,size(ZS2,2));
GCaMP6=[0,1.69644104899772,5.13796058542217,8.27886020152244,10.3756715204800,11.8173714529814,12.2425184714093,10.8571417354877,8.80831829681196,6.91339112244670,5.46959264663869,4.30868766622567,3.42533619066766,2.75378443486879,2.18017250852183,1.72816235135824,1.32732537295463,1.00684435500268,0.730210038304555,0.530242444093118,0.362253250339685,0.227668255288566,0.0869242416152502,0.000718266708050853,-0.0828334873368325]';
for i=1:repetitions    
    Stimuli(1,(idxStart+(i-1)*interval):(idxStart+(i-1)*interval)+size(GCaMP6,1)-1)=GCaMP6;
    Stimuli(2,(idxStart+4+(i-1)*interval):(idxStart+4+(i-1)*interval)+size(GCaMP6,1)-1)=GCaMP6;
end

idxStart=840;
interval=40;
repetitions=10;
for i=1:repetitions    
    Stimuli(3,(idxStart+(i-1)*interval):(idxStart+(i-1)*interval)+size(GCaMP6,1)-1)=GCaMP6;
    Stimuli(4,(idxStart+4+(i-1)*interval):(idxStart+4+(i-1)*interval)+size(GCaMP6,1)-1)=GCaMP6;
end



ModelMultipower2=[];
parfor i=1:length(ZS2)    
    mdl=fitlm(Stimuli',ZS2(i,:));%,'interactions');
    ModelMultipower2(i).coef=mdl.Coefficients;    
    ModelMultipower2(i).rsquared=mdl.Rsquared.Adjusted;
end

idx_rsq=find([ModelMultipower2.rsquared]>0.1);
ZS2_rsq=ZS2(idx_rsq,:);

options = statset('UseParallel',1); [idxKmeans_ZS_rsq Cmap_ZS_rsq]=kmeans(ZS2_rsq,10,'Options',options,'Distance','cityblock','Replicates',5,'MaxIter',1000,'Display','final');
[Model_ZS,GoodBetas_ZS]=Test_Regress(Cmap_ZS_rsq,Stimuli,idxKmeans_ZS_rsq,0.2);
[Model,GoodBetas]=Test_Regress(Cmap_ZS,Stimuli,idxKmeans_ZS,0.4);

figure;
ha = tight_subplot(length(GoodBetas_ZS),1,[.01 .01],[.01 .01],[.01 .01]);
for i=1:length(GoodBetas_ZS)
    axes(ha(i));
    plot(mean(ZS2_rsq(find(idxKmeans_ZS_rsq==GoodBetas_ZS(i)),:),1));ylim([-3 5]);
    rectangle('EdgeColor','none','FaceColor','k','Position',[120 -3 0.5 7]);
end

FishList=unique(idx_Fish);
idx_Fish_unique=idx_Fish;
for fish_nb=1:length(FishList)
    idx_Fish_unique(idx_Fish==FishList(fish_nb))=fish_nb;
end
idx_Fish_unique_rsq=idx_Fish_unique(idx_rsq);

FigHandle=figure;
set(FigHandle, 'Position', [100, 100, 1800, 600]);
ha = tight_subplot(2,3,[.01 .01],[.01 .01],[.01 .01]);framerate=4;
x = linspace(1/framerate,size(ZS2,2)/framerate,size(ZS2,2));
colors=[0.7 0.4 1;0.14 1 0.14];
for i=1:2
    axes(ha(-2+i*3));
    idx_temp=find(idxKmeans_ZS_rsq==GoodBetas_ZS(i));
    idx_ctrl=find(idx_Position(idx_rsq(idx_temp))==1);
    idx_curare=find(idx_Position(idx_rsq(idx_temp))==2);
    meanToPlot=mean(ZS2_rsq(idx_temp(idx_ctrl),:),1);%,'color','k','LineWidth',3);ylim([-2 4]);%title(num2str(length(idx_ctrl)));
    std_95=std(ZS2_rsq(idx_temp(idx_ctrl),:),1,1);ylim([-3 5]);
    H=shadedErrorBar(x,meanToPlot,std_95);
    H.mainLine.Color=colors(i,:);
    H.mainLine.LineWidth=3;
    set(gca,'XTickLabel',[]);set(gca,'YTickLabel',[]);
    axes(ha(-1+i*3));
    meanToPlot=mean(ZS2_rsq(idx_temp(idx_curare),:),1);%,'color','k','LineWidth',3);ylim([-2 4]);%title(num2str(length(idx_curare)));
    std_95=std(ZS2_rsq(idx_temp(idx_curare),:),1,1);
    H=shadedErrorBar(x,meanToPlot,std_95);ylim([-3 5]);
    H.mainLine.Color=colors(i,:);
    H.mainLine.LineWidth=3;
    set(gca,'XTickLabel',[]);set(gca,'YTickLabel',[]);
    axes(ha(i*3));
    histogram(idx_Fish_unique_rsq(idx_temp(idx_ctrl)));
    hold on;
    histogram(idx_Fish_unique_rsq(idx_temp(idx_curare)));
    set(gca,'XTickLabel',[]);set(gca,'YTickLabel',[]);
%     mean_distrib(1)=mean(idx_Fish_unique_rsq(idx_temp(idx_ctrl)));
%     mean_distrib(2)=mean(idx_Fish_unique_rsq(idx_temp(idx_curare)));
%     std_distrib(1)=std(idx_Fish_unique_rsq(idx_temp(idx_ctrl)));
%     std_distrib(2)=std(idx_Fish_unique_rsq(idx_temp(idx_curare)));
%     barwitherr(std_distrib,mean_distrib);
end


meanToPlot=mean(LinReg.(regionName).ZS_rsq(idx_temp2,:),1);
        MeanClusPerBrain{counter3,j}=meanToPlot;
        std_95=std(LinReg.(regionName).ZS_rsq(idx_temp2,:),1,1);%confidence_intervals(LinReg.(regionName).ZS_rsq(idx_temp2,:),95);
        %std(LinReg.(regionName).ZS_rsq(idx_temp2,:),1,1);
        STDClusPerBrain{counter3,j}=std_95;
        H=shadedErrorBar(x(:),meanToPlot-mean(meanToPlot(1:4)), std_95);axis([0 190 -5 5]);
        H.mainLine.Color=colors{i}(j,:);
        H.patch.FaceColor=colors{i}(j,:)/2;
        H.edge(1).Color=colors{i}(j,:)/2;
        H.edge(2).Color=colors{i}(j,:)/2;
        set(gca,'XTickLabel',[]);set(gca,'YTickLabel',[]);