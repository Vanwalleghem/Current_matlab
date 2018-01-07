idx_Fish=zeros(1,length(MatFiles));
idx_Types=zeros(1,length(MatFiles));;
for i=1:length(MatFiles)    
    name=strcat(MatFiles(i).name);
    [Fish,~]=regexp(name,'Fish2017(\d+)_','tokens','match');Fish=str2num(Fish{1}{1});
    idx_Fish(i)=Fish;
     if strfind(name,'ELO_IRO')
        Fish=1;
    elseif strfind(name,'ILO_ERO')
        Fish=2;
     else
         Fish=0;
    end
    idx_Types(i)=Fish;
end
clearvars i Fish Plane name counter

figure;histogram(idx_Types)

%MatFiles=dir('*analysis_matlab.mat');
MatFiles=dir('*ELO_IRO*analysis_matlab.mat');
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

ZS=zscore(GoodCalcium+GoodNoise,1,2);

Numbers=[1 [MatFiles.GoodNumber]];
counter=1;
idx_Plane=nan(length(ZS),1);
idx_Fish=nan(length(ZS),1);
name=strcat(MatFiles(1).name);
for i=1:length(MatFiles)	
    name=strcat(MatFiles(i).name);
    [Plane,~]=regexp(name,'_(\d+)_','tokens','match');Plane=str2num(Plane{1}{1});    
    idx_Plane(Numbers(i):Numbers(i+1))=Plane;    
    [Fish,~]=regexp(name,'Fish2017(\d+)_','tokens','match');Fish=str2num(Fish{1}{1});
    idx_Fish(Numbers(i):Numbers(i+1))=Fish;
end
clearvars i Fish Plane name counter

[~,idx_fish_both,~]=intersect(idx_Fish,list_fish_both);


Stimuli=zeros(3,size(ZS,2));
start=42;
spike=[-0.104392135015146,1.69644104899772,5.13796058542217,8.27886020152244,10.3756715204800,11.8173714529814,12.2425184714093,10.8571417354877,8.80831829681196,6.91339112244670,5.46959264663869,4.30868766622567,3.42533619066766,2.75378443486879,2.18017250852183,1.72816235135824,1.32732537295463,1.00684435500268,0.730210038304555,0.530242444093118,0.362253250339685,0.227668255288566,0.0869242416152502,0.000718266708050853,-0.0828334873368325]';
for i = 1:3
    for j = 0:2
        Stimuli(i,start+j*120:start+j*120+length(spike)-1)=spike';
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
idx_rsq=find([ModelBothSides.rsquared]>0.1);
figure;histogram([ModelBothSides.rsquared]);
ZS_rsq=ZS(idx_rsq,:);

x = linspace(0.2,size(ZS,2)/5,size(ZS,2));y = linspace(1,size(ZS,1),size(ZS,1));
Fighandle=figure;
set(Fighandle, 'Position', [100, 100, 1200, 900]);
imagesc(x,y,ZS(randperm(size(ZS,1)),:),[-0.5 5]);colormap hot;set(gca,'YTickLabel',[]);
Fighandle=figure;
set(Fighandle, 'Position', [100, 100, 1200, 900]);
imagesc(x,y,ZS_rsq(randperm(size(ZS_rsq,1)),:),[-0.5 5]);colormap hot;set(gca,'YTickLabel',[]);

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

options = statset('UseParallel',1); [idxKmeans_ZS_rsq Cmap_ZS_rsq]=kmeans(ZS_rsq,20,'Options',options,'Distance','cityblock','Replicates',3,'MaxIter',1000,'Display','final');
[ModelCmap,GoodBetas]=Test_Regress(Cmap_ZS_rsq,Stimuli,idxKmeans_ZS_rsq,0.3);

options = statset('UseParallel',1); [idxKmeans_ZS_rsq_10 Cmap_ZS_rsq_10]=kmeans(ZS_rsq,10,'Options',options,'Distance','cityblock','Replicates',3,'MaxIter',1000,'Display','final');
[ModelCmap_10,GoodBetas_10]=Test_Regress(Cmap_ZS_rsq_10,Stimuli,idxKmeans_ZS_rsq_10,0.3);


options = statset('UseParallel',1); [idxKmeans_ZS_rsq_15 Cmap_ZS_rsq_15]=kmeans(ZS_rsq,15,'Options',options,'Distance','cityblock','Replicates',3,'MaxIter',1500,'Display','final');
[ModelCmap_15,GoodBetas_15]=Test_Regress(Cmap_ZS_rsq_15,Stimuli,idxKmeans_ZS_rsq_15,0.3);

options = statset('UseParallel',1); [idxKmeans_ZS_rsq_15b Cmap_ZS_rsq_15b]=kmeans(ZS_rsq,15,'Options',options,'Distance','correlation','Replicates',3,'MaxIter',1500,'Display','final');
[ModelCmap_15b,GoodBetas_15b]=Test_Regress(Cmap_ZS_rsq_15b,Stimuli,idxKmeans_ZS_rsq_15b,0.3);


idx_weird=find(coefficients(:,1)>0.05 & coefficients(:,2)<0);
figure;plot(mean(ZS(idx_weird,:),1));
idx_weird2=find(coefficients(:,1)<0 & coefficients(:,2)>0.05);
figure;plot(mean(ZS(idx_weird,:),1));hold on;plot(mean(ZS(idx_weird2,:),1));
figure;
subplot(1,2,1);imagesc(ZS(idx_weird,:),[-1 3]);colormap hot
subplot(1,2,2);imagesc(ZS(idx_weird2,:),[-1 3]);colormap hot

idxKmeans_final=zeros(size(ZS,1),1);
idxKmeans_final(idx_rsq)=idxKmeans_ZS_rsq;

temp=unique(idx_Fish);
for i=1:length(temp);
    idx_Fish(find(idx_Fish==temp(i)))=i;
end

GoodBetas=GoodBetas([1 5 11 17]);

counter=1;
x = linspace(0.5,size(ZS,2)/2,size(ZS,2));
rows=length(GoodBetas);
counter=1;
Fighandle=figure;
set(Fighandle, 'Position', [100, 100, 1300, 900]);xplot=4;
for i=GoodBetas
    idx_temp=find(idxKmeans_final==i);
    subplot(rows,xplot,counter);plot(mean(ZS(idx_temp,:),1));
    subplot(rows,xplot,counter+1);imagesc(ZS(idx_temp,:),[-0.5 4]);
    subplot(rows,xplot,counter+2);histogram(idx_Plane(idx_temp));    
    subplot(rows,xplot,counter+3);histogram(idx_Fish(idx_temp));%h = gca;h.XTickLabel={'1','2'};
    counter=counter+4;
end

options = statset('UseParallel',1); [idxCoef_10 Cmap_coef_10]=kmeans(coefficients(idx_rsq,:),10,'Replicates',3,'MaxIter',1500,'Display','final');

figure;
for i=1:size(Cmap_coef_10,1)
    subplot(5,2,i);plot(mean(ZS_rsq(find(idxCoef_10==i),:),1));
end

ZS_both=vertcat(ZS_rsq,ZS_rsq_center);
options = statset('UseParallel',1); [idxKmeans_ZS_rsq_15b Cmap_ZS_rsq_15b]=kmeans(ZS_both,15,'Options',options,'Distance','cityblock','Replicates',3,'MaxIter',1500,'Display','final');
[ModelCmap_15b,GoodBetas_15b]=Test_Regress(Cmap_ZS_rsq_15b,Stimuli,idxKmeans_ZS_rsq_15b,0.3);

idx_Position(1:length(ZS_rsq))=1;
idx_Position(length(ZS_rsq)+1:length(ZS_rsq)+length(ZS_rsq_center))=2;

GoodBetas_select=GoodBetas_15b([9 1 3 4 5 12 7 10 13 14]);
counter=1;
x = linspace(0.5,size(ZS,2)/2,size(ZS,2));
rows=length(GoodBetas_select);
counter=1;
Fighandle=figure;
set(Fighandle, 'Position', [100, 100, 1300, 900]);xplot=3;
for i=GoodBetas_select
    idx_temp=find(idxKmeans_ZS_rsq_15b==i);
    subplot(rows,xplot,counter);plot(mean(ZS_both(idx_temp,:),1));ylim([-2 4]);
    subplot(rows,xplot,counter+1);imagesc(ZS_both(idx_temp,:),[-2 4]);colormap hot;    
    subplot(rows,xplot,counter+2);histogram(idx_Position(idx_temp));%h = gca;h.XTickLabel={'1','2'};
    counter=counter+xplot;
end

Fighandle=figure;
set(Fighandle, 'Position', [100, 100, 1300, 900]);xplot=3;
counter=1;xplot=floor(sqrt(length(GoodBetas_select)));yplot=ceil(length(GoodBetas_select)/xplot);
for i=GoodBetas_select
    idx_temp=find(idxKmeans_ZS_rsq_15b==i);
    idx_exterior=find(idx_Position(idx_temp)==1);
    idx_center=find(idx_Position(idx_temp)==2);
    subplot(xplot,yplot,counter);plot(mean(ZS_both(idx_temp(idx_exterior),:),1));hold on;plot(mean(ZS_both(idx_temp(idx_center),:),1));hold off;
    counter=counter+1;
end