%load('Complex_Noise.mat');
MatFiles=dir('*_analysis_matlab.mat');
MatFiles = rmfield(MatFiles,'isdir');MatFiles = rmfield(MatFiles,'bytes');MatFiles = rmfield(MatFiles,'date');
MatFiles = rmfield(MatFiles,'datenum');
name=strcat(MatFiles(1).name);
Calcium=load(name, 'DenoisedTraces');
Calcium=Calcium.DenoisedTraces(:,1:2500);
MatFiles(1).number=size(Calcium,1);
Fitness=load(name, 'idx_components');
Fitness=Fitness.idx_components+1;
GoodCalcium=Calcium(Fitness,:);%should be (Fitness,:)
MatFiles(1).GoodNumber=length(Fitness);
for i = 2:length(MatFiles)
name=strcat(MatFiles(i).name);
C=load(name, 'DenoisedTraces');
C=C.DenoisedTraces(:,1:2500);
F=load(name, 'idx_components');
F=F.idx_components+1;
GC=C(F,:);
Fitness=horzcat(Fitness,F);
GoodCalcium=vertcat(GoodCalcium,GC);
MatFiles(i).number=size(Calcium,1);
MatFiles(i).GoodNumber=MatFiles(i-1).GoodNumber+length(F);
end
clearvars GC C S F N name i GS DF Noise Calcium Spikes;

ZS2=zscore(GoodCalcium,1,2);

options = statset('UseParallel',1); [idxKmeans_ZS Cmap_ZS]=kmeans(ZS2,50,'Options',options,'Distance','cityblock','Replicates',3,'MaxIter',1000,'Display','final');

NewFlow=zeros(6,size(ZS2,2));
GCaMP6=[-0.104392135015146,1.69644104899772,5.13796058542217,8.27886020152244,10.3756715204800,11.8173714529814,12.2425184714093,10.8571417354877,8.80831829681196,6.91339112244670,5.46959264663869,4.30868766622567,3.42533619066766,2.75378443486879,2.18017250852183,1.72816235135824,1.32732537295463,1.00684435500268,0.730210038304555,0.530242444093118,0.362253250339685,0.227668255288566,0.0869242416152502,0.000718266708050853,-0.0828334873368325]';
GCaMP6s=interp(GCaMP6,2);GCaMP6s=GCaMP6s/max(GCaMP6s);
GCaMP6s(GCaMP6s<0)=0;
back=    [56 256 557 1006 1106 1466 1827 2086]; %Infuse
back_off=[106 356 706 1056 1206 1516 1926 2176];
fwd=    [156 407 757 856 1256 1316 1526 1626 1986 2236]; %Withdraw
fwd_off=[206 506 806 956 1306 1366 1576 1776 2036 2286];
for i=1:length(back)
NewFlow(1,back(i):back(i)+size(GCaMP6s,1)-1)=GCaMP6s';
NewFlow(2,back_off(i):back_off(i)+size(GCaMP6s,1)-1)=GCaMP6s';
NewFlow(5,back(i):back_off(i))=1;
end
for i=1:length(fwd)
NewFlow(3,fwd(i):fwd(i)+size(GCaMP6s,1)-1)=GCaMP6s';
NewFlow(4,fwd_off(i):fwd_off(i)+size(GCaMP6s,1)-1)=GCaMP6s';
NewFlow(6,fwd(i):fwd_off(i))=1;
end
clearvars GCaMP6 back back_off fwd fwd_off;
NewFlow(7,:)=[0:size(ZS2,2):1];

NewFlow_basic=zeros(8,size(ZS2,2));
load('D:\Pictures\processed\Flow\BasicClusters.mat');
GCaMP6=[-0.104392135015146,1.69644104899772,5.13796058542217,8.27886020152244,10.3756715204800,11.8173714529814,12.2425184714093,10.8571417354877,8.80831829681196,6.91339112244670,5.46959264663869,4.30868766622567,3.42533619066766,2.75378443486879,2.18017250852183,1.72816235135824,1.32732537295463,1.00684435500268,0.730210038304555,0.530242444093118,0.362253250339685,0.227668255288566,0.0869242416152502,0.000718266708050853,-0.0828334873368325]';
GCaMP6s=interp(GCaMP6,2);GCaMP6s=GCaMP6s/max(GCaMP6s);
back=    [56 256 557 1006 1106 1466 1827 2086]; %Infuse
back_off=[106 356 706 1056 1206 1516 1926 2176];
fwd=    [156 407 757 856 1256 1316 1526 1626 1986 2236]; %Withdraw
fwd_off=[206 506 806 956 1306 1366 1576 1776 2036 2286];
for i=1:length(back)
    NewFlow_basic(1,back(i):back(i)+100)=Basic_Clusters(1,56:156)+abs(min(Basic_Clusters(1,56:156)));
    NewFlow_basic(2,back(i):back(i)+100)=Basic_Clusters(2,56:156)+abs(min(Basic_Clusters(2,56:156)));
    NewFlow_basic(3,back(i):back(i)+100)=Basic_Clusters(3,56:156)+abs(min(Basic_Clusters(3,56:156)));
    NewFlow_basic(7,back(i):back(i)+100)=Basic_Clusters(7,56:156)+abs(min(Basic_Clusters(7,56:156)));
    NewFlow_basic(4,back(i):back(i)+100)=Basic_Clusters(4,56:156)+abs(min(Basic_Clusters(4,56:156)));
    NewFlow_basic(5,back(i):back(i)+100)=Basic_Clusters(5,56:156)+abs(min(Basic_Clusters(5,56:156)));
    NewFlow_basic(6,back(i):back(i)+100)=Basic_Clusters(6,56:156)+abs(min(Basic_Clusters(6,56:156)));
    NewFlow_basic(8,back(i):back(i)+100)=Basic_Clusters(8,56:156)+abs(min(Basic_Clusters(8,56:156)));
end
for i=1:length(fwd)
    NewFlow_basic(1,fwd(i):fwd(i)+100)=Basic_Clusters(1,156:256)+abs(min(Basic_Clusters(1,156:256)));
    NewFlow_basic(2,fwd(i):fwd(i)+100)=Basic_Clusters(2,156:256)+abs(min(Basic_Clusters(2,156:256)));
    NewFlow_basic(3,fwd(i):fwd(i)+100)=Basic_Clusters(3,156:256)+abs(min(Basic_Clusters(3,156:256)));
    NewFlow_basic(7,fwd(i):fwd(i)+100)=Basic_Clusters(7,156:256)+abs(min(Basic_Clusters(7,156:256)));
    NewFlow_basic(4,fwd(i):fwd(i)+100)=Basic_Clusters(4,156:256)+abs(min(Basic_Clusters(4,156:256)));
    NewFlow_basic(5,fwd(i):fwd(i)+100)=Basic_Clusters(5,156:256)+abs(min(Basic_Clusters(5,156:256)));
    NewFlow_basic(6,fwd(i):fwd(i)+100)=Basic_Clusters(6,156:256)+abs(min(Basic_Clusters(6,156:256)));
    NewFlow_basic(8,fwd(i):fwd(i)+100)=Basic_Clusters(8,156:256)+abs(min(Basic_Clusters(8,156:256)));
end
clearvars GCaMP6 GCaMP6s back back_off fwd fwd_off fwd_basic fwd_long fwd_three NewFlow2
NewFlow_basic(9,:)=[0:1/size(ZS2,2):1-1/size(ZS2,2)];

PredictFlow=zeros(8,size(ZS2,2));
GCaMP6=[-0.104392135015146,1.69644104899772,5.13796058542217,8.27886020152244,10.3756715204800,11.8173714529814,12.2425184714093,10.8571417354877,8.80831829681196,6.91339112244670,5.46959264663869,4.30868766622567,3.42533619066766,2.75378443486879,2.18017250852183,1.72816235135824,1.32732537295463,1.00684435500268,0.730210038304555,0.530242444093118,0.362253250339685,0.227668255288566,0.0869242416152502,0.000718266708050853,-0.0828334873368325]';
GCaMP6s=interp(GCaMP6,2);GCaMP6s=GCaMP6s/max(GCaMP6s);
back=      [56 1006 1466]; %Infuse
back_long= [256 1106 1827 2086];
back_three=557;
back_off=[106 356 706 1056 1206 1516 1926 2176];
fwd=    [156 757 1256 1316 1526 1986 2236]; %Withdraw
fwd_long= [407 856];
fwd_three=1626;
fwd_off=[206 506 806 956 1306 1366 1576 1776 2036 2286];
for clust_nb = 1:8
    for i=1:length(back)
        PredictFlow(clust_nb,back(i):back(i)+100)=Basic_Clusters(clust_nb,56:156)+abs(min(Basic_Clusters(clust_nb,56:156)));
    end
    for i=1:length(back_long)
        PredictFlow(clust_nb,back_long(i):back_long(i)+201)=interp(Basic_Clusters(clust_nb,56:156)+abs(min(Basic_Clusters(clust_nb,56:156))),2);
    end
    PredictFlow(clust_nb,back_three:back_three+302)=interp(Basic_Clusters(clust_nb,56:156)+abs(min(Basic_Clusters(clust_nb,56:156))),3);
    for i=1:length(fwd)
        PredictFlow(clust_nb,fwd(i):fwd(i)+100)=Basic_Clusters(clust_nb,156:256)+abs(min(Basic_Clusters(clust_nb,156:256)));
    end
    for i=1:length(fwd_long)
        PredictFlow(clust_nb,fwd_long(i):fwd_long(i)+201)=interp(Basic_Clusters(clust_nb,156:256)+abs(min(Basic_Clusters(clust_nb,156:256))),2);
    end
    PredictFlow(clust_nb,fwd_three:fwd_three+302)=interp(Basic_Clusters(clust_nb,156:256)+abs(min(Basic_Clusters(clust_nb,156:256))),3);
end
clearvars GCaMP6 GCaMP6s back back_off fwd fwd_off back_long back_three fwd_long;
PredictFlow(9,:)=[0:1/size(ZS2,2):1-1/size(ZS2,2)];

parfor i=1:size(ZS2,1)
    mdl=fitlm(NewFlow_basic',ZS2(i,:));    
    model_basic(i).coef=mdl.Coefficients;        
    model_basic(i).rsquared=mdl.Rsquared.Adjusted;
end
clearvars slow_fast slow_fast_fwd i

parfor i=1:size(ZS2,1)
    mdl=fitlm(PredictFlow',ZS2(i,:));    
    model_predict(i).coef=mdl.Coefficients;        
    model_predict(i).rsquared=mdl.Rsquared.Adjusted;
end
clearvars slow_fast slow_fast_fwd i

figure;
subplot(3,1,1);histogram([model_NewFlow.rsquared]);
subplot(3,1,2);histogram([model_predict.rsquared]);
subplot(3,1,3);histogram([model_basic.rsquared]);

idx_rsq_NewFlow=find([model_NewFlow.rsquared]>0.1);
idx_rsq_predict=find([model_predict.rsquared]>0.1);
idx_rsq_basic=find([model_basic.rsquared]>0.1);

[idxKmeans_NewFlow Cmap_NewFlow]=kmeans(ZS2(idx_rsq_NewFlow,:),40,'Options',options,'Distance','cityblock','Replicates',3,'MaxIter',1000,'Display','final');
[idxKmeans_predict Cmap_predict]=kmeans(ZS2(idx_rsq_predict,:),40,'Options',options,'Distance','cityblock','Replicates',3,'MaxIter',1000,'Display','final');
[idxKmeans_basic Cmap_basic]=kmeans(ZS2(idx_rsq_basic,:),40,'Options',options,'Distance','cityblock','Replicates',3,'MaxIter',1000,'Display','final');

figure;
subplot(3,1,1);imagesc(Cmap_NewFlow);
subplot(3,1,2);imagesc(Cmap_basic);
subplot(3,1,3);imagesc(Cmap_predict);

Threshold=0.3;
[Model_NewFlow,GoodBetas_NewFlow]=Test_Regress(Cmap_NewFlow,NewFlow,idxKmeans_NewFlow,Threshold);
[Model_basic,GoodBetas_basic]=Test_Regress(Cmap_basic,NewFlow_basic,idxKmeans_basic,Threshold);
[Model_predict,GoodBetas_predict]=Test_Regress(Cmap_predict,PredictFlow,idxKmeans_predict,Threshold);


idx_rsq_all=union(union(idx_rsq_NewFlow,idx_rsq_predict),idx_rsq_basic);
ZS_rsq=ZS2(idx_rsq_all,:);

Corr_BasicCLust=zeros(size(Basic_Clusters,1),size(ZS2,1));
start_complex=[155 1005];start_basic=[155 258];duration=90;

parfor i=1:size(ZS2,1)
    corr_temp=[];
    ZS_rsq_temp=detrend([ZS2(i,start_complex(1):start_complex(1)+duration) ZS2(i,start_complex(2):start_complex(2)+duration)]);
    for j=1:size(Basic_Clusters,1)
        Basic_temp=[Basic_Clusters(j,start_basic(1):start_basic(1)+duration) Basic_Clusters(j,start_basic(2):start_basic(2)+duration)];
        temp=corrcoef(ZS_rsq_temp, Basic_temp);
        %temp=pdist([ZS_rsq_temp; Basic_temp],'cityblock');
        corr_temp(j)=temp(1,2);        
        %corr_temp(j)=temp;
    end
    Corr_BasicCLust(:,i)=corr_temp;
end

[MaxCorr_BasiClust MaxCorr_BasiClust_ind]=max(abs(Corr_BasicCLust),[],1);
idx_BasicClust_corr=find(MaxCorr_BasiClust>0.4);
BasicClustData=[];
for basic=1:max(MaxCorr_BasiClust_ind)
    idx_BasicClust_temp=find(MaxCorr_BasiClust_ind(idx_BasicClust_corr)==basic);
    BasicClustData(basic).ZS=ZS2(idx_BasicClust_corr(idx_BasicClust_temp),:);
    BasicClustData(basic).mean=mean(BasicClustData(basic).ZS,1);      
    ZS_temp=BasicClustData(basic).ZS;
    [BasicClustData(basic).idxKmeans BasicClustData(basic).Cmap]=kmeans(ZS_temp,10,'Options',options,'Distance','cityblock','Replicates',3,'MaxIter',1000,'Display','final');
    [BasicClustData(basic).Model,BasicClustData(basic).GoodBetas]=Test_Regress(BasicClustData(basic).Cmap,NewFlow_basic,BasicClustData(basic).idxKmeans,Threshold);
end
clearvars basic ZS_temp idx_BasicClust_temp

Fighandle=figure;
set(Fighandle, 'Position', [100, 100, 1300, 900]);
counter=1;xplot=floor(sqrt(max(MaxCorr_BasiClust_ind)));yplot=ceil(max(MaxCorr_BasiClust_ind)/xplot);
for i=1:max(MaxCorr_BasiClust_ind)    
    subplot(xplot,yplot,counter);plot(x,BasicClustData(i).mean);xlim([0 150]); %hold on;plot(x,(double(Flow_profile2)/10)-1)
    counter=counter+1;
end

Fighandle=figure;
set(Fighandle, 'Position', [100, 100, 1300, 900]);
counter=1;xplot=floor(sqrt(max(MaxCorr_BasiClust_ind)));yplot=ceil(max(MaxCorr_BasiClust_ind)/xplot);
for i=1:max(MaxCorr_BasiClust_ind)    
    subplot(xplot,yplot,counter);imagesc(BasicClustData(i).ZS);xlim([0 500]); %hold on;plot(x,(double(Flow_profile2)/10)-1)
    counter=counter+1;
end

GoodBetas_NewFlow_select=GoodBetas_NewFlow([1 2 3 4 5 6 9 12 16 17 18 19 22 23 25 26]);
GoodBetas_NewFlow_select=GoodBetas_NewFlow_select([1 4 5 6 7 8 9 10 12 13 14 15 16]);
Fighandle=figure;
set(Fighandle, 'Position', [100, 100, 1300, 900]);
counter=1;xplot=4;yplot=ceil(length(GoodBetas_NewFlow_select)/2);
ha=tight_subplot(yplot,xplot);
for i=1:length(GoodBetas_NewFlow_select)
    idx_temp=find(idxKmeans_NewFlow==GoodBetas_NewFlow_select(i));
    axes(ha(counter));
    temp=mean(ZS2(idx_rsq_NewFlow(idx_temp),:),1);
    std_temp=std(ZS2(idx_rsq_NewFlow(idx_temp),:),1,1);
    H=shadedErrorBar(x, temp, std_temp);axis([0 450 -2 5]);title(num2str((i)));hold on;
    for jj=1:length(back)
        rectangle('EdgeColor','k','FaceColor','g','Position',[back(jj)/5 -1 abs(back(jj)-back_off(jj))/5 0.25]);
    end
    for jj=1:length(fwd)
        rectangle('EdgeColor','k','FaceColor','m','Position',[fwd(jj)/5 -1 abs(fwd(jj)-fwd_off(jj))/5 0.25]);
    end    
    hold off;
    axes(ha(counter+1));
    imagesc(ZS2(idx_rsq_NewFlow(idx_temp),:),[-2 5]);colormap hot;
    counter=counter+2;    
end

GoodBetas_basic_select=GoodBetas_basic([1 2 4 9 11 12 15 16 17 19 21 22 23 25 26 27 29 30]);
GoodBetas_basic_select=GoodBetas_basic_select([1 2 3 4 5 6 7 8 11 12 14 15 16 17 18]);
Fighandle=figure;
set(Fighandle, 'Position', [100, 100, 1300, 900]);
counter=1;xplot=4;yplot=ceil(length(GoodBetas_basic_select)/2);
ha=tight_subplot(yplot,xplot);
for i=1:length(GoodBetas_basic_select)
    idx_temp=find(idxKmeans_basic==GoodBetas_basic_select(i));
    axes(ha(counter));
    temp=mean(ZS2(idx_rsq_basic(idx_temp),:),1);
    std_temp=std(ZS2(idx_rsq_basic(idx_temp),:),1,1);
    H=shadedErrorBar(x, temp, std_temp);axis([0 300 -2 5]);title(num2str((i)));hold on;
    for jj=1:length(back)
        rectangle('EdgeColor','k','FaceColor','g','Position',[back(jj)/5 -1 abs(back(jj)-back_off(jj))/5 0.25]);
    end
    for jj=1:length(fwd)
        rectangle('EdgeColor','k','FaceColor','m','Position',[fwd(jj)/5 -1 abs(fwd(jj)-fwd_off(jj))/5 0.25]);
    end    
    hold off;
    axes(ha(counter+1));
    imagesc(ZS2(idx_rsq_basic(idx_temp),1:150*5),[-2 5]);colormap hot;
    counter=counter+2;
end

GoodBetas_predict_select=GoodBetas_predict([1 2 6 8 10 12 13 15 17 18 19 20 24]);
%GoodBetas_predict_select=GoodBetas_predict_select([2 3 4 5 7 10]);
Fighandle=figure;
set(Fighandle, 'Position', [100, 100, 1300, 900]);
counter=1;xplot=4;yplot=ceil(length(GoodBetas_predict_select)/2);
ha=tight_subplot(yplot,xplot);
for i=1:length(GoodBetas_predict_select)
    idx_temp=find(idxKmeans_predict==GoodBetas_predict_select(i));
    axes(ha(counter));
    temp=mean(ZS2(idx_rsq_predict(idx_temp),:),1);
    std_temp=std(ZS2(idx_rsq_predict(idx_temp),:),1,1);
    H=shadedErrorBar(x, temp, std_temp);axis([0 150 -2 5]);title(num2str((i)));hold on;
    for jj=1:length(back)
        rectangle('EdgeColor','k','FaceColor','g','Position',[back(jj)/5 -1 abs(back(jj)-back_off(jj))/5 0.25]);
    end
    for jj=1:length(fwd)
        rectangle('EdgeColor','k','FaceColor','m','Position',[fwd(jj)/5 -1 abs(fwd(jj)-fwd_off(jj))/5 0.25]);
    end    
    hold off;
    axes(ha(counter+1));
    imagesc(ZS2(idx_rsq_predict(idx_temp),1:150*5),[-2 5]);colormap hot;
    counter=counter+2;
end


%Test Speed encoding
Speed_flow=zeros(2,size(ZS2,2));
back=    [56 256 557 1006 1106 1466 1827 2086]-5; %Infuse
back_off=[106 356 706 1056 1206 1516 1926 2176];
fwd=    [156 407 757 856 1256 1316 1526 1626 1986 2236]-5; %Withdraw
fwd_off=[206 506 806 956 1306 1366 1576 1776 2036 2286];
%back=back/5;back_off=back_off/5;
%fwd=fwd/5;fwd_off=fwd_off/5;
slow_fast=[1 2 3 2 1 2 2 1];
slow_fast_fwd=[2 1 1 2 2 2 2 3 1 1];
for i=1:length(back)
    Speed_flow(1,back(i):back_off(i))=slow_fast(i);
end
for i=1:length(fwd)    
    Speed_flow(2,fwd(i):fwd_off(i))=slow_fast_fwd(i);
end


idx_temp=find(Speed_flow(1,:)==3);
Speed_flow(1,idx_temp(50):idx_temp(100))=2;
idx_temp=find(Speed_flow(2,:)==3);
Speed_flow(2,idx_temp(50):idx_temp(100))=2;
Speed_flow(Speed_flow==3)=1;

Max_back=zeros(length(back),size(ZS2,1));
for i=1:length(back)
    Max_back(i,:)=max(ZS2(:,back(i):back_off(i)),[],2);    
end
Max_fwd=zeros(length(fwd),size(ZS2,1));
for i=1:length(fwd)    
    Max_fwd(i,:)=max(ZS2(:,fwd(i):fwd_off(i)),[],2);    
end

Fighandle=figure;
set(Fighandle, 'Position', [100, 100, 1300, 900]);
counter=1;xplot=4;yplot=ceil(length(GoodBetas_NewFlow_select)/2);
ha=tight_subplot(yplot,xplot);
for i=1:length(GoodBetas_NewFlow_select)
    idx_temp=find(idxKmeans_NewFlow==GoodBetas_NewFlow_select(i));    
    temp=mean(ZS2(idx_rsq_NewFlow(idx_temp),1:end),1);
    std_temp=std(ZS2(idx_rsq_NewFlow(idx_temp),1:end),1,1);
    idx_temp=find(idxKmeans_NewFlow==GoodBetas_NewFlow_select(i));
    idx_speed=find(sum(Speed_flow,1)>0);
    axes(ha(counter));
    H=shadedErrorBar(x(1:length(idx_speed)), temp(idx_speed), std_temp(idx_speed));axis([0 300 -2 5]);title(num2str((i)));hold on;
    plot(x(1:length(idx_speed)),Speed_flow(1,idx_speed),'Color','m');hold on;
    plot(x(1:length(idx_speed)),Speed_flow(2,idx_speed),'Color','g');hold off;
    axes(ha(counter+1));
    imagesc(ZS2(idx_rsq_NewFlow(idx_temp),idx_speed),[-2 5]);colormap hot;
    counter=counter+2;
end

Fighandle=figure;
set(Fighandle, 'Position', [100, 100, 1300, 900]);
counter=1;xplot=4;yplot=ceil(length(GoodBetas_basic_select)/2);
ha=tight_subplot(yplot,xplot);
for i=1:length(GoodBetas_basic_select)
    idx_temp=find(idxKmeans_basic==GoodBetas_basic_select(i));
    axes(ha(counter));
    temp=mean(ZS2(idx_rsq_basic(idx_temp),:),1);
    std_temp=std(ZS2(idx_rsq_basic(idx_temp),:),1,1);    
    idx_speed=find(sum(Speed_flow,1)>0);
    axes(ha(counter));
    H=shadedErrorBar(x(1:length(idx_speed)), temp(idx_speed), std_temp(idx_speed));axis([0 300 -2 5]);title(num2str((i)));hold on;
    plot(x(1:length(idx_speed)),Speed_flow(1,idx_speed),'Color','m');hold on;
    plot(x(1:length(idx_speed)),Speed_flow(2,idx_speed),'Color','g');hold off;
    axes(ha(counter+1));
    imagesc(ZS2(idx_rsq_basic(idx_temp),idx_speed),[-2 5]);colormap hot;
    counter=counter+2;
end

for i=1:numel(GoodBetas_NewFlow_select)
    idx_temp=find(idxKmeans_NewFlow==GoodBetas_NewFlow_select(i));    
    corr_temp=zeros(length(idx_temp),1);
    mean_temp=mean(ZS2(idx_rsq_NewFlow(idx_temp),1:end),1);
    for j=1:length(idx_temp)
        temp=corrcoef(mean_temp, ZS2(idx_rsq_NewFlow(idx_temp(j)),:));
        corr_temp(j)=temp(1,2);
    end
    NewFlowData(i).CorrCoef=corr_temp;
end

GoodClusters_goodmembers_NewFlow=[];Threshold=0.4;
idxKmeans_ZS_goodmembers=zeros(1,size(ZS2,1));
for i=1:length(GoodBetas_NewFlow_select)    
    GoodClusters_goodmembers_NewFlow(i).ZS=ZS2(idx_rsq_NewFlow(find(NewFlowData(i).CorrCoef>=Threshold)),:);
    temp=find(idxKmeans_NewFlow==GoodBetas_NewFlow_select(i));
    GoodClusters_goodmembers_NewFlow(i).idx=temp(find(NewFlowData(i).CorrCoef>=Threshold));
    GoodClusters_goodmembers_NewFlow(i).mean=mean(GoodClusters_goodmembers_NewFlow(i).ZS,1);
    GoodClusters_goodmembers_NewFlow(i).STD=std(GoodClusters_goodmembers_NewFlow(i).ZS,1,1);
    idx=find(idxKmeans_NewFlow==GoodBetas_NewFlow_select(i));
    idx=idx(find(NewFlowData(i).CorrCoef>=Threshold));
    idxKmeans_ZS_goodmembers(idx_rsq_NewFlow(idx))=GoodBetas_NewFlow_select(i);    
end



Fighandle=figure;
set(Fighandle, 'Position', [100, 100, 1300, 900]);
counter=1;xplot=2;yplot=2;
ha=tight_subplot(yplot,xplot);
idx_temp=find(idxKmeans_NewFlow==GoodBetas_NewFlow_select(4) |  idxKmeans_NewFlow==GoodBetas_NewFlow_select(6));    
temp=mean(ZS2(idx_rsq_NewFlow(idx_temp),1:end),1);
std_temp=std(ZS2(idx_rsq_NewFlow(idx_temp),1:end),1,1);
idx_speed=find(sum(Speed_flow,1)>0);
axes(ha(counter));
H=shadedErrorBar(x(1:length(idx_speed)), temp(idx_speed), std_temp(idx_speed));axis([0 300 -2 5]);title(length(idx_temp));hold on;
plot(x(1:length(idx_speed)),Speed_flow(1,idx_speed),'Color','m');hold on;
plot(x(1:length(idx_speed)),Speed_flow(2,idx_speed),'Color','g');hold off;
axes(ha(counter+1));
imagesc(ZS2(idx_rsq_NewFlow(idx_temp),idx_speed),[-2 5]);colormap hot;
counter=counter+2;
idx_temp=find(idxKmeans_NewFlow==GoodBetas_NewFlow_select(5));    
temp=mean(ZS2(idx_rsq_NewFlow(idx_temp),1:end),1);
std_temp=std(ZS2(idx_rsq_NewFlow(idx_temp),1:end),1,1);
idx_speed=find(sum(Speed_flow,1)>0);
axes(ha(counter));
H=shadedErrorBar(x(1:length(idx_speed)), temp(idx_speed), std_temp(idx_speed));axis([0 300 -2 5]);title(length(idx_temp));hold on;
plot(x(1:length(idx_speed)),Speed_flow(1,idx_speed),'Color','m');hold on;
plot(x(1:length(idx_speed)),Speed_flow(2,idx_speed),'Color','g');hold off;
axes(ha(counter+1));
imagesc(ZS2(idx_rsq_NewFlow(idx_temp),idx_speed),[-2 5]);colormap hot;

%Kmeans on the max response to fwd flow (4 slow, 5 fast, 1 mix)
[idxKmeans_MaxFwd Cmap_MaxFwd]=kmeans(Max_fwd([2 3 9 10 1 4 5 6 7 8],idx_rsq_NewFlow)',10,'Options',options,'Replicates',3,'MaxIter',1000,'Display','final');
%Kmeans on the max response to bwd flow (3 slow, 4 fast, 1 mix)
[idxKmeans_MaxBwd Cmap_MaxBwd]=kmeans(Max_back([1 5 8 2 4 6 7 3],idx_rsq_NewFlow)',10,'Options',options,'Replicates',3,'MaxIter',1000,'Display','final');

Fighandle=figure;
set(Fighandle, 'Position', [100, 100, 1300, 900]);
idx_speed=idx_rsq_NewFlow(find(idxKmeans_MaxFwd==6));
temp=mean(ZS2((idx_speed),1:end),1);
std_temp=std(ZS2((idx_speed),1:end),1,1);
H=shadedErrorBar(x, temp, std_temp);axis([0 450 -2 5]);title(length(idx_temp));hold on;
plot(x,Speed_flow(1,:),'Color','m');hold on;
plot(x,Speed_flow(2,:),'Color','g');hold off;

Fighandle=figure;
set(Fighandle, 'Position', [100, 100, 1300, 900]);
idx_speed=idx_rsq_NewFlow(find(idxKmeans_MaxBwd==5));
temp=mean(ZS2((idx_speed),1:end),1);
std_temp=std(ZS2((idx_speed),1:end),1,1);
H=shadedErrorBar(x, temp, std_temp);axis([0 450 -2 5]);title(length(idx_temp));hold on;
plot(x,Speed_flow(1,:),'Color','m');hold on;
plot(x,Speed_flow(2,:),'Color','g');hold off;

idx_temp=find(idxKmeans_NewFlow==GoodBetas_NewFlow_select(4) |  idxKmeans_NewFlow==GoodBetas_NewFlow_select(6));    
Prism_temp=Max_back(:,idx_rsq_NewFlow(idx_temp))';
csvwrite('__Max_back_4-6_NewFlow.csv',Prism_temp);
Prism_temp=Max_fwd(:,idx_rsq_NewFlow(idx_temp))';
csvwrite('__Max_fwd_4-6_NewFlow.csv',Prism_temp);
Prism_temp=[];
idx_temp=find(idxKmeans_NewFlow==GoodBetas_NewFlow_select(5));    
Prism_temp=Max_back(:,idx_rsq_NewFlow(idx_temp))';
csvwrite('__Max_back_5_NewFlow.csv',Prism_temp);
Prism_temp=Max_fwd(:,idx_rsq_NewFlow(idx_temp))';
csvwrite('__Max_fwd_5_NewFlow.csv',Prism_temp);

Fighandle=figure;
set(Fighandle, 'Position', [100, 100, 1300, 900]);
counter=1;xplot=2;yplot=2;
ha=tight_subplot(yplot,xplot);
idx_temp=find(idxKmeans_ZS_goodmembers==GoodBetas_NewFlow_select(4) |  idxKmeans_ZS_goodmembers==GoodBetas_NewFlow_select(6));    
temp=mean(ZS2((idx_temp),1:end),1);
std_temp=std(ZS2((idx_temp),1:end),1,1);
idx_speed=find(sum(Speed_flow,1)>0);
axes(ha(counter));
H=shadedErrorBar(x(1:length(idx_speed)), temp(idx_speed), std_temp(idx_speed));axis([0 300 -2 5]);title(length(idx_temp));hold on;
plot(x(1:length(idx_speed)),Speed_flow(1,idx_speed),'Color','m');hold on;
plot(x(1:length(idx_speed)),Speed_flow(2,idx_speed),'Color','g');hold off;
axes(ha(counter+1));
imagesc(ZS2((idx_temp),idx_speed),[-2 5]);colormap hot;
counter=counter+2;
idx_temp=find(idxKmeans_ZS_goodmembers==GoodBetas_NewFlow_select(5));    
temp=mean(ZS2((idx_temp),1:end),1);
std_temp=std(ZS2((idx_temp),1:end),1,1);
idx_speed=find(sum(Speed_flow,1)>0);
axes(ha(counter));
H=shadedErrorBar(x(1:length(idx_speed)), temp(idx_speed), std_temp(idx_speed));axis([0 300 -2 5]);title(length(idx_temp));hold on;
plot(x(1:length(idx_speed)),Speed_flow(1,idx_speed),'Color','m');hold on;
plot(x(1:length(idx_speed)),Speed_flow(2,idx_speed),'Color','g');hold off;
axes(ha(counter+1));
imagesc(ZS2((idx_temp),idx_speed),[-2 5]);colormap hot;

idx_temp=find(idxKmeans_ZS_goodmembers==GoodBetas_NewFlow_select(4) |  idxKmeans_ZS_goodmembers==GoodBetas_NewFlow_select(6));    
Prism_temp=Max_back(:,(idx_temp))';
csvwrite('__Max_back_4-6_NewFlow.csv',Prism_temp);
Prism_temp=Max_fwd(:,(idx_temp))';
csvwrite('__Max_fwd_4-6_NewFlow.csv',Prism_temp);
Prism_temp=[];
idx_temp=find(idxKmeans_ZS_goodmembers==GoodBetas_NewFlow_select(5));    
Prism_temp=Max_back(:,(idx_temp))';
csvwrite('__Max_back_5_NewFlow.csv',Prism_temp);
Prism_temp=Max_fwd(:,(idx_temp))';
csvwrite('__Max_fwd_5_NewFlow.csv',Prism_temp);

Fighandle=figure;
set(Fighandle, 'Position', [100, 100, 1300, 900]);
counter=1;xplot=2;yplot=2;
ha=tight_subplot(yplot,xplot);
axes(ha(counter));
idx_speed=idx_rsq_NewFlow(find(idxKmeans_MaxFwd==6));
temp=mean(ZS2((idx_speed),1:end),1);
std_temp=std(ZS2((idx_speed),1:end),1,1);
H=shadedErrorBar(x, temp, std_temp);axis([0 450 -2 5]);title(length(idx_temp));hold on;
plot(x,Speed_flow(1,:),'Color','m');hold on;
plot(x,Speed_flow(2,:),'Color','g');hold off;
axes(ha(counter+1));
imagesc(ZS2((idx_speed),1:end),[-2 5]);colormap hot;
counter=counter+2;
axes(ha(counter));
idx_speed=idx_rsq_NewFlow(find(idxKmeans_MaxBwd==5));
temp=mean(ZS2((idx_speed),1:end),1);
std_temp=std(ZS2((idx_speed),1:end),1,1);
H=shadedErrorBar(x, temp, std_temp);axis([0 450 -2 5]);title(length(idx_temp));hold on;
plot(x,Speed_flow(1,:),'Color','m');hold on;
plot(x,Speed_flow(2,:),'Color','g');hold off;
axes(ha(counter+1));
imagesc(ZS2((idx_speed),1:end),[-2 5]);colormap hot;

%Alternative look for speed encoding
slow_fast=[1 2 3 2 1 2 2 1];
slow_fast_fwd=[2 1 1 2 2 2 2 3 1 1];
Fast_speed_idx=slow_fast==2;
Slow_speed_idx=slow_fast==1;

Speed_encode_idx=find(mean(Max_back(Fast_speed_idx,:),1)>(mean(Max_back(Slow_speed_idx,:),1)+2*std(Max_back(Slow_speed_idx,:),1,1)));
figure;
imagesc(ZS2(Speed_encode_idx,:),[-1 4]);

Slow_Speed_encode_idx=find(mean(Max_back(Slow_speed_idx,:),1)>(mean(Max_back(Fast_speed_idx,:),1)+2*std(Max_back(Fast_speed_idx,:),1,1)));
figure;
imagesc(ZS2(Slow_Speed_encode_idx,:),[-1 4]);

[idxKmeans_FastBack Cmap_FastBack]=kmeans(ZS2(Speed_encode_idx,:),10,'Options',options,'Distance','cityblock','Replicates',3,'MaxIter',1000,'Display','final');
Threshold=0.2;
[Model_FastBack,GoodBetas_FastBack]=Test_Regress(Cmap_FastBack,NewFlow,idxKmeans_FastBack,Threshold);

[idxKmeans_SlowBack Cmap_SlowBack]=kmeans(ZS2(Slow_Speed_encode_idx,:),10,'Options',options,'Distance','cityblock','Replicates',3,'MaxIter',1000,'Display','final');
Threshold=0.2;
[Model_SlowBack,GoodBetas_SlowBack]=Test_Regress(Cmap_SlowBack,NewFlow,idxKmeans_SlowBack,Threshold);


Fighandle=figure;
set(Fighandle, 'Position', [100, 100, 1300, 900]);
counter=1;xplot=2;yplot=2;
ha=tight_subplot(yplot,xplot);
axes(ha(counter));
idx_speed=Speed_encode_idx(find(idxKmeans_FastBack==GoodBetas_FastBack));
temp=mean(ZS2((idx_speed),1:end),1);
std_temp=std(ZS2((idx_speed),1:end),1,1);
H=shadedErrorBar(x, temp, std_temp);axis([0 450 -2 5]);title(length(idx_temp));hold on;
plot(x,Speed_flow(1,:),'Color','m');hold on;
plot(x,Speed_flow(2,:),'Color','g');hold off;
axes(ha(counter+1));
imagesc(ZS2((idx_speed),1:end),[-2 5]);colormap hot;
counter=counter+2;
axes(ha(counter));
idx_speed=Slow_Speed_encode_idx(find(idxKmeans_SlowBack==GoodBetas_SlowBack));
temp=mean(ZS2((idx_speed),1:end),1);
std_temp=std(ZS2((idx_speed),1:end),1,1);
H=shadedErrorBar(x, temp, std_temp);axis([0 450 -2 5]);title(length(idx_temp));hold on;
plot(x,Speed_flow(1,:),'Color','m');hold on;
plot(x,Speed_flow(2,:),'Color','g');hold off;
axes(ha(counter+1));
imagesc(ZS2((idx_speed),1:end),[-2 5]);colormap hot;

Speed_encode_idx_fwd=find(mean(Max_fwd_clear(Fast_speed_idx,:),1)>(mean(Max_fwd_clear(Slow_speed_idx,:),1)+2*std(Max_fwd_clear(Slow_speed_idx,:),1,1)));
figure;
imagesc(ZS2(Speed_encode_idx_fwd,:),[-1 4]);

Slow_Speed_encode_idx_fwd=find(mean(Max_fwd_clear(Slow_speed_idx,:),1)>(mean(Max_fwd_clear(Fast_speed_idx,:),1)+2*std(Max_fwd_clear(Fast_speed_idx,:),1,1)));
figure;
imagesc(ZS2(Slow_Speed_encode_idx_fwd,:),[-1 4]);


idx_temp=find(Speed_flow(1,:)==3);
Speed_flow(1,idx_temp(50):idx_temp(100))=2;
idx_temp=find(Speed_flow(2,:)==3);
Speed_flow(2,idx_temp(50):idx_temp(100))=2;
Speed_flow(Speed_flow==3)=1;

counter=1;
Max_flow=zeros(length(back)+length(fwd),size(ZS2,1));
for i=1:length(back)
    Max_flow(counter,:)=max(ZS2(:,back(i):back_off(i)),[],2);    
    counter=counter+1;
end

for i=1:length(fwd)    
    Max_flow(counter,:)=max(ZS2(:,fwd(i):fwd_off(i)),[],2);    
    counter=counter+1;
end

%Kmeans on the max response to fwd flow (4 slow, 5 fast, 1 mix)
[idxKmeans_MaxFlow Cmap_MaxFlow]=kmeans(Max_flow([2 3 9 10 11 15 18 1 4 5 6 7 12 14 16 17 8 13],idx_rsq_NewFlow)',20,'Options',options,'Replicates',3,'MaxIter',1000,'Display','final');
%Kmeans on the max response to bwd flow (3 slow, 4 fast, 1 mix)


%Test Speed encoding
Speed_flow=zeros(2,size(ZS2,2));
back=    [56 256 557 1006 1106 1466 1827 2086]-5; %Infuse
back_off=[106 356 706 1056 1206 1516 1926 2176];
fwd=    [156 407 757 856 1256 1316 1526 1626 1986 2236]-5; %Withdraw
fwd_off=[206 506 806 956 1306 1366 1576 1776 2036 2286];
%back=back/5;back_off=back_off/5;
%fwd=fwd/5;fwd_off=fwd_off/5;
slow_fast=[1 2 3 2 1 2 2 1];
slow_fast_fwd=[2 1 1 2 2 2 2 3 1 1];
for i=1:length(back)
    Speed_flow(1,back(i):back_off(i))=slow_fast(i);
end
for i=1:length(fwd)    
    Speed_flow(2,fwd(i):fwd_off(i))=slow_fast_fwd(i);
end


Speed_encode_idx=find(mean(Max_flow(Fast_speed_idx,:),1)>(mean(Max_flow(Slow_speed_idx,:),1)+2*std(Max_flow(Slow_speed_idx,:),1,1)));
figure;
imagesc(ZS2(Speed_encode_idx,:),[-1 4]);

Slow_Speed_encode_idx=find(mean(Max_flow(Slow_speed_idx,:),1)>(mean(Max_flow(Fast_speed_idx,:),1)+2*std(Max_flow(Fast_speed_idx,:),1,1)));
figure;
imagesc(ZS2(Slow_Speed_encode_idx,:),[-1 4]);

[idxKmeans_FastBack Cmap_FastBack]=kmeans(ZS2(Speed_encode_idx,:),20,'Options',options,'Distance','cityblock','Replicates',3,'MaxIter',1000,'Display','final');
Threshold=0.4;
[Model_FastBack,GoodBetas_FastBack]=Test_Regress(Cmap_FastBack,NewFlow,idxKmeans_FastBack,Threshold);

[idxKmeans_SlowBack Cmap_SlowBack]=kmeans(ZS2(Slow_Speed_encode_idx,:),20,'Options',options,'Distance','cityblock','Replicates',3,'MaxIter',1000,'Display','final');
Threshold=0.4;
[Model_SlowBack,GoodBetas_SlowBack]=Test_Regress(Cmap_SlowBack,NewFlow,idxKmeans_SlowBack,Threshold);

Fighandle=figure;
set(Fighandle, 'Position', [100, 100, 1300, 900]);
counter=1;xplot=2;yplot=2;
ha=tight_subplot(yplot,xplot);
axes(ha(counter));
idx_speed=Speed_encode_idx(find(idxKmeans_FastBack==GoodBetas_FastBack));
temp=mean(ZS2((idx_speed),1:end),1);
std_temp=std(ZS2((idx_speed),1:end),1,1);
H=shadedErrorBar(x, temp, std_temp);axis([0 450 -2 5]);title(length(idx_speed));hold on;
plot(x,Speed_flow(1,:),'Color','m');hold on;
plot(x,Speed_flow(2,:),'Color','g');hold off;
axes(ha(counter+1));
imagesc(ZS2(idx_speed(randperm(length(idx_speed))),1:end),[-2 5]);colormap hot;
counter=counter+2;
axes(ha(counter));
idx_speed=Slow_Speed_encode_idx(find(idxKmeans_SlowBack==GoodBetas_SlowBack));
temp=mean(ZS2((idx_speed),1:end),1);
std_temp=std(ZS2((idx_speed),1:end),1,1);
H=shadedErrorBar(x, temp, std_temp);axis([0 450 -2 5]);title(length(idx_speed));hold on;
plot(x,Speed_flow(1,:),'Color','m');hold on;
plot(x,Speed_flow(2,:),'Color','g');hold off;
axes(ha(counter+1));
imagesc(ZS2(idx_speed(randperm(length(idx_speed))),1:end),[-2 5]);colormap hot;

Speed_encoders=[];Threshold=0.4;
idx_speed=Speed_encode_idx(find(idxKmeans_FastBack==GoodBetas_FastBack));
mean_temp=mean(ZS2((idx_speed),1:end),1);
ZS_temp=ZS2((idx_speed),1:end);
corr_temp=zeros(size(ZS_temp,1),1);
parfor i=1:length(idx_speed)
	temp=corrcoef(mean_temp, ZS_temp(i,:));
	corr_temp(i)=temp(1,2);
end
Speed_encoders(1).idx=idx_speed(corr_temp>=Threshold);
Speed_encoders(1).mean=mean(ZS2(Speed_encoders.idx,:),1);
Speed_encoders(1).std=std(ZS2(Speed_encoders.idx,:),1,1);

idx_speed=Slow_Speed_encode_idx(find(idxKmeans_SlowBack==GoodBetas_SlowBack));
mean_temp=mean(ZS2((idx_speed),1:end),1);
ZS_temp=ZS2((idx_speed),1:end);
corr_temp=zeros(size(ZS_temp,1),1);
parfor i=1:length(idx_speed)
	temp=corrcoef(mean_temp, ZS_temp(i,:));
	corr_temp(i)=temp(1,2);
end
Speed_encoders(2).idx=idx_speed(corr_temp>=Threshold);
Speed_encoders(2).mean=mean(ZS2(Speed_encoders.idx,:),1);
Speed_encoders(2).std=std(ZS2(Speed_encoders.idx,:),1,1);

Fighandle=figure;
set(Fighandle, 'Position', [100, 100, 800, 800]);
idx_speed=Speed_encoders(1).idx;
imagesc(ZS2(idx_speed(randperm(length(idx_speed))),1:end),[-2 5]);colormap hot;

Fighandle=figure;
set(Fighandle, 'Position', [100, 100, 1500, 800]);
temp=mean(ZS2((idx_speed),1:end),1);
std_temp=std(ZS2((idx_speed),1:end),1,1);
H=shadedErrorBar(x, temp, std_temp);axis([0 475 -2 6]);title(length(idx_speed));hold on;
plot(x,(Speed_flow(1,:)/2)-2.01,'Color','m');hold on;
plot(x,(Speed_flow(2,:)/2)-2.01,'Color','g');hold off;

Fighandle=figure;
set(Fighandle, 'Position', [100, 100, 800, 800]);
idx_speed=Speed_encoders(2).idx;
imagesc(ZS2(idx_speed(randperm(length(idx_speed))),1:end),[-2 5]);colormap hot;

Fighandle=figure;
set(Fighandle, 'Position', [100, 100, 1500, 800]);
temp=mean(ZS2((idx_speed),1:end),1);
std_temp=std(ZS2((idx_speed),1:end),1,1);
H=shadedErrorBar(x, temp, std_temp);axis([0 475 -2 6]);title(length(idx_speed));hold on;
plot(x,(Speed_flow(1,:)/2)-2.01,'Color','m');hold on;
plot(x,(Speed_flow(2,:)/2)-2.01,'Color','g');hold off;


idx_speed=Speed_encoders(2).idx;
counter=1;
Max_flow_speed=zeros(length(back)+length(fwd),length(idx_speed));
for i=1:length(back)
    Max_flow_speed(counter,:)=max(ZS2(idx_speed,back(i):back_off(i)),[],2);    
    counter=counter+1;
end
for i=1:length(fwd)    
    Max_flow_speed(counter,:)=max(ZS2(idx_speed,fwd(i):fwd_off(i)),[],2);    
    counter=counter+1;
end
Max_flow_speed=Max_flow_speed';

%idx_Fish and ROIs
%idx_Fish and ROIs
Numbers=[0 [MatFiles.GoodNumber]];
counter=1;
idx_Plane=nan(length(ZS2),1);
idx_Fish=nan(length(ZS2),1);
for i=1:length(MatFiles)
    name=strcat(MatFiles(i).name);
    [Plane,~]=regexp(name,'_(\d+)um_','tokens','match');Plane=str2num(Plane{1}{1});
    if ~isempty(regexp(name,'2p_complex_f(\d)proper_','tokens'))
        [Fish,~]=regexp(name,'2p_complex_f(\d)proper_.+(\d)_output','tokens','match');
        if isempty(Fish)
            Fish=regexp(name,'2p_complex_f(\d)proper_','tokens');
            Fish=strcat(Fish{1}{1},'1234','0');
            Fish=str2num(Fish);
        else
            Fish=strcat(Fish{1}{1},'1234',Fish{1}{2});
            Fish=str2num(Fish);
        end
    elseif ~isempty(regexp(name,'Fish(\d)_','tokens'))
        [Fish,~]=regexp(name,'Fish(\d)_.+_(\d)_output','tokens','match');
        if isempty(Fish)
            Fish=regexp(name,'Fish(\d)_','tokens');
            Fish=strcat(Fish{1}{1},'1234','0');
            Fish=str2num(Fish);
        else
            Fish=strcat(Fish{1}{1},'1234',Fish{1}{2});
            Fish=str2num(Fish);
        end
    else
        [Fish,~]=regexp(name,'F(\d+)_.+_(\d)_output','tokens','match');
        if isempty(Fish)
            Fish=regexp(name,'F(\d+)_','tokens');
            Fish=strcat(Fish{1}{1},'1234','0');
            Fish=str2num(Fish);
        else
            Fish=strcat(Fish{1}{1},'1234',Fish{1}{2});
            Fish=str2num(Fish);
        end
        
    end
    idx_Plane(Numbers(i)+1:Numbers(i+1))=Plane;
    idx_Fish(Numbers(i)+1:Numbers(i+1))=Fish;
end
clearvars i Fish Plane name counter

All_ROIs=[];
ROIs_idx=[];
for i = 1:length(MatFiles)
    name=strcat(MatFiles(i).name);
    Rs=load(name, 'ROIs');
    Rs=Rs.ROIs;
    F=load(name, 'idx_components');
    F=F.idx_components+1;
    Rs=Rs(:,F);
    All_ROIs{i}=Rs;
    if i==1
        ROIs_idx(i)=length(F);
    else
        ROIs_idx(i)=ROIs_idx(i-1)+length(F);
    end
end
clearvars GC C S F N name i;


Fish_list=unique(idx_Fish);
Errored_ROI={};
progressbar(0,0,0);
for fish_nb=1:length(Fish_list)
    IndexC=find(idx_Fish==Fish_list(fish_nb));
    IndexC=ismember(ROIs_idx,IndexC);
    MatFiles_fish = find(IndexC>0);
    if MatFiles_fish(1)==1
        numbersForROIs=[1 [MatFiles(MatFiles_fish).GoodNumber]];
    else
        numbersForROIs=[MatFiles(MatFiles_fish(1)-1).GoodNumber [MatFiles(MatFiles_fish).GoodNumber]];
    end
    Centroids=[];
    counter=1;
    for plane=1:length(MatFiles_fish)
        filename=MatFiles(MatFiles_fish(plane)).name;
        [slice,~]=regexp(filename,'_(\d+)um_','tokens','match');slice=str2num(slice{1}{1});
        idx_name=strcat(num2str(Fish_list(fish_nb)),'_',num2str(slice));
         if plane==1
            temp_roi=0;
        else
            temp_roi=temp_roi+size(ROI,3);
        end
        ROI=All_ROIs{MatFiles_fish(plane)};
        imagename=regexp(filename,'_output_analysis','split');
        imagename=strcat(imagename{1},'_mean.tiff');
        image=double(imread(imagename));image=image/max(max(image));image=image*128;
        ROI=reshape(full(ROI),size(image,1),size(image,2),size(ROI,2));                
        for roi_nb=1:size(ROI,3)
            progressbar([],[],roi_nb/size(ROI,3));
            temp=regionprops(uint16(squeeze(ROI(:,:,roi_nb)))==max(max(uint16(squeeze(ROI(:,:,roi_nb))))),'Centroid');
            Centroids(roi_nb+temp_roi,5)=roi_nb+temp_roi;
            temp=temp.Centroid;
            Centroids(roi_nb+temp_roi,1:2)=temp;
            Centroids(roi_nb+temp_roi,3)=(slice/20);
        end     
        progressbar([],fish_nb/length(MatFiles_fish));
    end
    if iscell(Fish_list)
        image_name=strcat('_ROIsFish_',Fish_list{fish_nb},'.csv');
    else
        image_name=strcat('_ROIsFish_',num2str(Fish_list(fish_nb)),'.csv');
    end
    csvwrite(image_name,Centroids);
    progressbar(fish_nb/length(Fish_list));
end
clearvars i j k ROI Centroids M I I_row I_col test temp image_name filename fish fish_nb

idx_speed=Speed_encoders(1).idx;
ROI_speed=ROI_rotated(idx_speed,:);
ROI_speed(:,4)=1;
csvwrite('Speed_encoding_ROIs.csv',ROI_speed);

idx_speed=Speed_encoders(2).idx;
ROI_speed=ROI_rotated(idx_speed,:);
ROI_speed(:,4)=1;
csvwrite('Slow_encoding_ROIs.csv',ROI_speed);

%Onset encodes change of speed
[idxKmeans_Speed Cmap_Speed]=kmeans(ZS2(:,sum(Speed_flow,1)>0),40,'Options',options,'Distance','cityblock','Replicates',3,'MaxIter',1000,'Display','final');
[Model_OnsetSpeed,GoodBetas_OnsetSpeed]=Test_Regress(Cmap_Speed,NewFlow(:,sum(Speed_flow,1)>0),idxKmeans_Speed,Threshold);



%Test Acceleration encoding
Accel_flow=zeros(2,size(ZS2,2));
back=    [56 256 557 1006 1106 1466 1827 2086]-5; %Infuse
back_off=[106 356 706 1056 1206 1516 1926 2176];
fwd=    [156 407 757 856 1256 1316 1526 1626 1986 2236]-5; %Withdraw
fwd_off=[206 506 806 956 1306 1366 1576 1776 2036 2286];
%back=back/5;back_off=back_off/5;
%fwd=fwd/5;fwd_off=fwd_off/5;
Decel_fast=[1 2 3 2 1 2 2 1];
Decel_fast_fwd=[2 1 1 2 2 2 2 3 1 1];
for i=1:length(back)
    Accel_flow(1,back(i):back_off(i))=Decel_fast(i);
end
for i=1:length(fwd)    
    Accel_flow(2,fwd(i):fwd_off(i))=Decel_fast_fwd(i);
end

Accel_times=[];
temp=find(Accel_flow(1,:)==3);
Accel_times(1,1)=temp(find(Speed_flow(1,temp)==1,1));
Accel_times(1,2)=temp(find(Speed_flow(1,temp)==2,1));
temp=find(Accel_flow(2,:)==3);
Accel_times(2,1)=temp(find(Speed_flow(2,temp)==1,1));
Accel_times(2,2)=temp(find(Speed_flow(2,temp)==2,1));

Max_back_accel=zeros(2,size(ZS2,1));
Max_fwd_accel=zeros(2,size(ZS2,1));
for i=1:2
    Max_back_accel(i,:)=max(ZS2(:,Accel_times(1,i):Accel_times(1,i)+40),[],2);    
    Max_fwd_accel(i,:)=max(ZS2(:,Accel_times(2,i):Accel_times(2,i)+40),[],2);    
end

Threshold=0.4;
Accel_encode_idx=find(Max_back_accel(1,:)<Max_back_accel(2,:));
Decel_encode_idx=find(Max_back_accel(1,:)>Max_back_accel(2,:));
[idxKmeans_AccelBack Cmap_AccelBack]=kmeans(ZS2(Accel_encode_idx,1:1000),10,'Options',options,'Distance','cityblock','Replicates',3,'MaxIter',1000,'Display','final');Threshold=0.2;
[Model_AccelBack,GoodBetas_AccelBack]=Test_Regress(Cmap_AccelBack,NewFlow(:,1:1000),idxKmeans_AccelBack,Threshold);
[idxKmeans_DecelBack Cmap_DecelBack]=kmeans(ZS2(Decel_encode_idx,1:1000),10,'Options',options,'Distance','cityblock','Replicates',3,'MaxIter',1000,'Display','final');
[Model_DecelBack,GoodBetas_DecelBack]=Test_Regress(Cmap_DecelBack,NewFlow(:,1:1000),idxKmeans_DecelBack,Threshold);

Accel_FWD_encode_idx=find(Max_fwd_accel(1,:)<Max_fwd_accel(2,:));
Decel_FWD_encode_idx=find(Max_fwd_accel(1,:)>Max_fwd_accel(2,:));

[idxKmeans_AccelFwd Cmap_AccelFwd]=kmeans(ZS2(Accel_FWD_encode_idx,1200:end),10,'Options',options,'Distance','cityblock','Replicates',3,'MaxIter',1000,'Display','final');
[Model_AccelFwd,GoodBetas_AccelFwd]=Test_Regress(Cmap_AccelFwd,NewFlow(:,1200:end),idxKmeans_AccelFwd,Threshold);
[idxKmeans_DecelFwd Cmap_DecelFwd]=kmeans(ZS2(Decel_FWD_encode_idx,1200:end),10,'Options',options,'Distance','cityblock','Replicates',3,'MaxIter',1000,'Display','final');
[Model_DecelFwd,GoodBetas_DecelFwd]=Test_Regress(Cmap_DecelFwd,NewFlow(:,1200:end),idxKmeans_DecelFwd,Threshold);

Fighandle=figure;
set(Fighandle, 'Position', [100, 100, 800, 800]);
ha=tight_subplot(1,length(GoodBetas_AccelBack));
for i=1:length(GoodBetas_AccelBack)
    axes(ha(i));
    plot(Cmap_AccelBack(GoodBetas_AccelBack(i),:));hold on;plot(sum(Speed_flow(:,1:1000),1)/5);
end

Fighandle=figure;
set(Fighandle, 'Position', [100, 100, 1600, 800]);
ha=tight_subplot(1,length(GoodBetas_DecelBack));
for i=1:length(GoodBetas_DecelBack)
    axes(ha(i));
    plot(Cmap_DecelBack(GoodBetas_DecelBack(i),:));hold on;plot(sum(Speed_flow(:,1:1000),1)/5);
end

Fighandle=figure;
set(Fighandle, 'Position', [100, 100, 1600, 800]);
ha=tight_subplot(1,length(GoodBetas_AccelFwd));
for i=1:length(GoodBetas_AccelFwd)
    axes(ha(i));
    plot(Cmap_AccelFwd(GoodBetas_AccelFwd(i),:));hold on;plot(sum(Speed_flow(:,1200:end),1)/5);
end

Fighandle=figure;
set(Fighandle, 'Position', [100, 100, 1600, 800]);
ha=tight_subplot(1,length(GoodBetas_DecelFwd));
for i=1:length(GoodBetas_DecelFwd)
    axes(ha(i));
    plot(Cmap_DecelFwd(GoodBetas_DecelFwd(i),:));hold on;plot(sum(Speed_flow(:,1200:end),1)/5);
end

Onset_encoders=[];Threshold=0.4;
idx_speed=Decel_encode_idx(find(idxKmeans_DecelBack==GoodBetas_DecelBack));
mean_temp=mean(ZS2((idx_speed),1:end),1);
ZS_temp=ZS2((idx_speed),1:end);
corr_temp=zeros(size(ZS_temp,1),1);
parfor i=1:length(idx_speed)
	temp=corrcoef(mean_temp, ZS_temp(i,:));
	corr_temp(i)=temp(1,2);
end
Onset_encoders(1).idx=idx_speed(corr_temp>=Threshold);
Onset_encoders(1).mean=mean(ZS2(Onset_encoders.idx,:),1);
Onset_encoders(1).std=std(ZS2(Onset_encoders.idx,:),1,1);

idx_speed=Decel_FWD_encode_idx(find(idxKmeans_DecelFwd==GoodBetas_DecelFwd));
mean_temp=mean(ZS2((idx_speed),1:end),1);
ZS_temp=ZS2((idx_speed),1:end);
corr_temp=zeros(size(ZS_temp,1),1);
parfor i=1:length(idx_speed)
	temp=corrcoef(mean_temp, ZS_temp(i,:));
	corr_temp(i)=temp(1,2);
end
Onset_encoders(2).idx=idx_speed(corr_temp>=Threshold);
Onset_encoders(2).mean=mean(ZS2(Onset_encoders(2).idx,:),1);
Onset_encoders(2).std=std(ZS2(Onset_encoders(2).idx,:),1,1);

Fighandle=figure;
set(Fighandle, 'Position', [100, 100, 800, 800]);
idx_speed=Onset_encoders(1).idx;
imagesc(ZS2(idx_speed(randperm(length(idx_speed))),1:end),[-2 5]);colormap hot;

Fighandle=figure;
set(Fighandle, 'Position', [100, 100, 1500, 800]);
temp=mean(ZS2((idx_speed),1:end),1);
std_temp=std(ZS2((idx_speed),1:end),1,1);
H=shadedErrorBar(x, temp, std_temp);axis([0 475 -2 6]);title(length(idx_speed));hold on;
plot(x,(Speed_flow(1,:)/2)-2.01,'Color','m');hold on;
plot(x,(Speed_flow(2,:)/2)-2.01,'Color','g');hold off;

Fighandle=figure;
set(Fighandle, 'Position', [100, 100, 800, 800]);
idx_speed=Onset_encoders(2).idx;
imagesc(ZS2(idx_speed(randperm(length(idx_speed))),1:end),[-2 5]);colormap hot;

Fighandle=figure;
set(Fighandle, 'Position', [100, 100, 1500, 800]);
temp=mean(ZS2((idx_speed),1:end),1);
std_temp=std(ZS2((idx_speed),1:end),1,1);
H=shadedErrorBar(x, temp, std_temp);axis([0 475 -2 6]);title(length(idx_speed));hold on;
plot(x,(Speed_flow(1,:)/2)-2.01,'Color','m');hold on;
plot(x,(Speed_flow(2,:)/2)-2.01,'Color','g');hold off;

%Integrators

coefficients={}; %%%to make the coefficients variable that we will use. Regression coefficients represent the mean change in the response variable for one unit of change in the predictor variable while holding other predictors in the model constant.
for idx=1:length(model_predict)%%% to make a variable the size of ModelResults
    coef=[model_predict(idx).coef];%%%% and then put in another variable the coef field from ModelResults
    temp=coef.Properties.RowNames;temp=regexp(temp,'x(\d+)','tokens');%%%to take the name of the rows of the coef variable
    if ~isempty(temp)%%% if temp is not empty...
        %temp=[temp{:}];temp=[temp{:}];temp=[temp{:}];%temp=str2num(temp);
        for coef_idx=2:height(coef)%%%take the number of rows from coef, except the first one(i think because is the intercept)
            %if coef.pValue(coef_idx)<0.05%%%to select the coef that are bellow the p value we want, in this case 0.05
                coefficients{idx,str2num(temp{coef_idx}{1}{1})}=coef.Estimate(coef_idx); %%%to make an array the size of idx,10 with the coefficient values that were significant
            %end
        end
    end
end
idxempty=cellfun('isempty',coefficients); %%%to make a variable with where we will aply in every cell the isempty function wich will help us find the empty places
coefficients(idxempty)={0}; %%% and put a 0 in the places where we found that there were empty cells
clearvars idxempty idx coef_idx coef  %%%clear variables
coefficients=cell2mat(coefficients); %%%to make a matrix of the coefficients array
mean_coef=mean(coefficients,1);
std_coef=std(coefficients,1,1);


coef_rsq=coefficients(idx_rsq_all,:);
Coef_clust=zeros(8,size(ZS2,2));
Coef_ClustData=struct();
for i=1:8
    Coef_clust(i,:)=mean(ZS2(coefficients(:,i)>(mean_coef(i)+2*std_coef(i)),:),1);
    Coef_ClustData(i).idx=find(coefficients(:,i)>(mean_coef(i)+2*std_coef(i)));
end

Fighandle=figure;
set(Fighandle, 'Position', [100, 100, 1300, 900]);
counter=1;xplot=4;yplot=4;
ha=tight_subplot(yplot,xplot);
for i=1:8
    idx_temp=Coef_ClustData(i).idx;
    axes(ha(counter));
    temp=mean(ZS2((idx_temp),:),1);
    std_temp=std(ZS2((idx_temp),:),1,1);    
    idx_speed=find(sum(Speed_flow,1)>0);
    axes(ha(counter));
    H=shadedErrorBar(x, temp, std_temp);axis([0 300 -2 5]);title(num2str((i)));hold on;
    plot(x,Speed_flow(1,:),'Color','m');hold on;
    plot(x,Speed_flow(2,:),'Color','g');hold off;
    axes(ha(counter+1));
    imagesc(ZS2((idx_temp),:),[-2 5]);colormap hot;
    counter=counter+2;
end

%bidirectionnal Onset
figure;
plot(BasicClustData(1).Cmap(BasicClustData(1).GoodBetas(5),:)');
hold on;i=17;idx_temp=find(idxKmeans_basic==GoodBetas_basic_select(i));temp=mean(ZS2(idx_rsq_basic(idx_temp),:),1);plot(temp);
hold on; i=13;idx_temp=find(idxKmeans_NewFlow==GoodBetas_NewFlow_select(i));plot(mean(ZS2(idx_rsq_NewFlow(idx_temp),:),1));
hold on;i=6;idx_temp=find(idxKmeans_basic==GoodBetas_basic_select(i));temp=mean(ZS2(idx_rsq_basic(idx_temp),:),1);plot(temp);
hold on;i=11;idx_temp=find(idxKmeans_predict==GoodBetas_predict_select(i));plot(mean(ZS2(idx_rsq_predict(idx_temp),:),1));

%bidirectionnal ON
figure;
plot(BasicClustData(1).Cmap(BasicClustData(1).GoodBetas(6),:)');hold on;plot(BasicClustData(2).Cmap(BasicClustData(2).GoodBetas(5),:)');

%FWD Onset
figure;
plot(BasicClustData(3).Cmap(BasicClustData(3).GoodBetas(4),:)');hold on; i=2;idx_temp=find(idxKmeans_NewFlow==GoodBetas_NewFlow_select(i));plot(mean(ZS2(idx_rsq_NewFlow(idx_temp),:),1));
hold on; i=9;idx_temp=find(idxKmeans_NewFlow==GoodBetas_NewFlow_select(i));plot(mean(ZS2(idx_rsq_NewFlow(idx_temp),:),1));
hold on;i=1;idx_temp=find(idxKmeans_basic==GoodBetas_basic_select(i));temp=mean(ZS2(idx_rsq_basic(idx_temp),:),1);plot(temp);

%FWD ON
figure;
plot(BasicClustData(4).Cmap(BasicClustData(4).GoodBetas(5),:)');hold on;

%BWD Onset
figure;
plot(BasicClustData(6).Cmap(BasicClustData(6).GoodBetas(1),:)');hold on;
i=14;idx_temp=find(idxKmeans_basic==GoodBetas_basic_select(i));temp=mean(ZS2(idx_rsq_basic(idx_temp),:),1);plot(temp);
hold on; i=10;idx_temp=find(idxKmeans_NewFlow==GoodBetas_NewFlow_select(i));plot(mean(ZS2(idx_rsq_NewFlow(idx_temp),:),1));
hold on;i=6;;idx_temp=find(idxKmeans_predict==GoodBetas_predict_select(i));plot(mean(ZS2(idx_rsq_predict(idx_temp),:),1));

%BWD ON
figure;
plot(BasicClustData(7).Cmap(BasicClustData(7).GoodBetas(2),:)');hold on;

%BWD integrator
figure;
i=16;idx_temp=find(idxKmeans_basic==GoodBetas_basic_select(i));temp=mean(ZS2(idx_rsq_basic(idx_temp),:),1);plot(temp);
hold on; i=12;idx_temp=find(idxKmeans_NewFlow==GoodBetas_NewFlow_select(i));plot(mean(ZS2(idx_rsq_NewFlow(idx_temp),:),1));
hold on;i=7;idx_temp=find(idxKmeans_predict==GoodBetas_predict_select(i));plot(mean(ZS2(idx_rsq_predict(idx_temp),:),1));
hold on;
    plot(Speed_flow(1,:),'Color','m');hold on;
    plot(Speed_flow(2,:),'Color','g');hold off;
