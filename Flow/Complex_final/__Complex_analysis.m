load('Complex_Noise.mat');
options = statset('UseParallel',1); [idxKmeans_ZS Cmap_ZS]=kmeans(ZS2,50,'Options',options,'Distance','cityblock','Replicates',3,'MaxIter',1000,'Display','final');

NewFlow=zeros(8,size(ZS2,2));
GCaMP6=[-0.104392135015146,1.69644104899772,5.13796058542217,8.27886020152244,10.3756715204800,11.8173714529814,12.2425184714093,10.8571417354877,8.80831829681196,6.91339112244670,5.46959264663869,4.30868766622567,3.42533619066766,2.75378443486879,2.18017250852183,1.72816235135824,1.32732537295463,1.00684435500268,0.730210038304555,0.530242444093118,0.362253250339685,0.227668255288566,0.0869242416152502,0.000718266708050853,-0.0828334873368325]';
GCaMP6s=interp(GCaMP6,2);
back=    [56 256 557 1006 1106 1466 1827 2086]-10; %Infuse
back_off=[106 356 706 1056 1206 1516 1926 2176]-10;
fwd=    [156 407 757 856 1256 1316 1526 1626 1986 2236]-10; %Withdraw
fwd_off=[206 506 806 956 1306 1366 1576 1776 2036 2286]-10;
for i=1:length(back)
NewFlow(1,back(i):back(i)+size(GCaMP6,1)-1)=GCaMP6';
NewFlow(2,back(i):back(i)+size(GCaMP6s,1)-1)=GCaMP6s';
NewFlow(3,back_off(i):back_off(i)+size(GCaMP6,1)-1)=GCaMP6';
NewFlow(7,back(i):back_off(i))=1;
end
for i=1:length(fwd)
NewFlow(4,fwd(i):fwd(i)+size(GCaMP6,1)-1)=GCaMP6';
NewFlow(5,fwd(i):fwd(i)+size(GCaMP6s,1)-1)=GCaMP6s';
NewFlow(6,fwd_off(i):fwd_off(i)+size(GCaMP6,1)-1)=GCaMP6';
NewFlow(8,fwd(i):fwd_off(i))=1;
end
clearvars GCaMP6 back back_off fwd fwd_off;

Threshold=0.2;
[Model_ZS,GoodBetas_ZS]=Test_Regress(Cmap_ZS,NewFlow,idxKmeans_ZS,Threshold);

parfor i=1:size(ZS2,1)
    mdl=fitlm(NewFlow',ZS2(i,:));    
    model_NewFlow(i).coef=mdl.Coefficients;        
    model_NewFlow(i).rsquared=mdl.Rsquared.Adjusted;
end
clearvars slow_fast slow_fast_fwd i

NewFlow_basic=zeros(8,size(ZS2,2));
load('D:\Pictures\processed\Flow\BasicClusters.mat');
GCaMP6=[-0.104392135015146,1.69644104899772,5.13796058542217,8.27886020152244,10.3756715204800,11.8173714529814,12.2425184714093,10.8571417354877,8.80831829681196,6.91339112244670,5.46959264663869,4.30868766622567,3.42533619066766,2.75378443486879,2.18017250852183,1.72816235135824,1.32732537295463,1.00684435500268,0.730210038304555,0.530242444093118,0.362253250339685,0.227668255288566,0.0869242416152502,0.000718266708050853,-0.0828334873368325]';
GCaMP6s=interp(GCaMP6,2);
back=    [56 256 557 1006 1106 1466 1827 2086]-10; %Infuse
back_off=[106 356 706 1056 1206 1516 1926 2176]-10;
fwd=    [156 407 757 856 1256 1316 1526 1626 1986 2236]-10; %Withdraw
fwd_off=[206 506 806 956 1306 1366 1576 1776 2036 2286]-10;
for i=1:length(back)
    NewFlow_basic(1,back(i):back(i)+100)=Basic_Clusters(6,156:256)+abs(min(Basic_Clusters(6,156:256)));
    NewFlow_basic(2,back(i):back(i)+100)=Basic_Clusters(7,156:256)+abs(min(Basic_Clusters(7,156:256)));
    NewFlow_basic(3,back_off(i):back_off(i)+size(GCaMP6,1)-1)=GCaMP6';
    NewFlow_basic(7,back(i):back(i)+100)=Basic_Clusters(8,156:256)+abs(min(Basic_Clusters(8,156:256)));
end
for i=1:length(fwd)
    NewFlow_basic(4,fwd(i):fwd(i)+100)=Basic_Clusters(6,156:256)+abs(min(Basic_Clusters(6,156:256)));
    NewFlow_basic(5,fwd(i):fwd(i)+100)=Basic_Clusters(7,156:256)+abs(min(Basic_Clusters(7,156:256)));
    NewFlow_basic(6,fwd_off(i):fwd_off(i)+size(GCaMP6,1)-1)=GCaMP6';
    NewFlow_basic(8,fwd(i):fwd(i)+100)=Basic_Clusters(8,156:256)+abs(min(Basic_Clusters(8,156:256)));
end
clearvars GCaMP6 GCaMP6s back back_off fwd fwd_off fwd_basic fwd_long fwd_three NewFlow2

PredictFlow=zeros(8,size(ZS2,2));
GCaMP6=[-0.104392135015146,1.69644104899772,5.13796058542217,8.27886020152244,10.3756715204800,11.8173714529814,12.2425184714093,10.8571417354877,8.80831829681196,6.91339112244670,5.46959264663869,4.30868766622567,3.42533619066766,2.75378443486879,2.18017250852183,1.72816235135824,1.32732537295463,1.00684435500268,0.730210038304555,0.530242444093118,0.362253250339685,0.227668255288566,0.0869242416152502,0.000718266708050853,-0.0828334873368325]';
GCaMP6s=interp(GCaMP6,2);
back=      [56 1006 1466]-10; %Infuse
back_long= [256 1106 1827 2086]-10;
back_three=557-10;
back_off=[106 356 706 1056 1206 1516 1926 2176]-10;
fwd=    [156 757 1256 1316 1526 1986 2236]-10; %Withdraw
fwd_long= [407 856]-10;
fwd_three=1626-10;
fwd_off=[206 506 806 956 1306 1366 1576 1776 2036 2286]-10;
for i=1:length(back)
    PredictFlow(1,back(i):back(i)+100)=Basic_Clusters(6,156:256)+abs(min(Basic_Clusters(6,156:256)));
    PredictFlow(2,back(i):back(i)+100)=Basic_Clusters(7,156:256)+abs(min(Basic_Clusters(7,156:256)));    
    PredictFlow(7,back(i):back(i)+100)=Basic_Clusters(8,156:256)+abs(min(Basic_Clusters(8,156:256)));
end
for i=1:length(back_off)
    PredictFlow(3,back_off(i):back_off(i)+size(GCaMP6,1)-1)=GCaMP6';
end
for i=1:length(back_long)
    PredictFlow(1,back_long(i):back_long(i)+201)=interp(Basic_Clusters(6,156:256)+abs(min(Basic_Clusters(6,156:256))),2);
    PredictFlow(2,back_long(i):back_long(i)+201)=interp(Basic_Clusters(7,156:256)+abs(min(Basic_Clusters(7,156:256))),2);    
    PredictFlow(7,back_long(i):back_long(i)+201)=interp(Basic_Clusters(8,156:256)+abs(min(Basic_Clusters(8,156:256))),2);
end
PredictFlow(1,back_three:back_three+302)=interp(Basic_Clusters(6,156:256)+abs(min(Basic_Clusters(6,156:256))),3);
PredictFlow(2,back_three:back_three+302)=interp(Basic_Clusters(7,156:256)+abs(min(Basic_Clusters(7,156:256))),3);    
PredictFlow(7,back_three:back_three+302)=interp(Basic_Clusters(8,156:256)+abs(min(Basic_Clusters(8,156:256))),3);
for i=1:length(fwd)
    PredictFlow(4,fwd(i):fwd(i)+100)=Basic_Clusters(6,156:256)+abs(min(Basic_Clusters(6,156:256)));
    PredictFlow(5,fwd(i):fwd(i)+100)=Basic_Clusters(7,156:256)+abs(min(Basic_Clusters(7,156:256)));
    PredictFlow(6,fwd_off(i):fwd_off(i)+size(GCaMP6,1)-1)=GCaMP6';
    PredictFlow(8,fwd(i):fwd(i)+100)=Basic_Clusters(8,156:256)+abs(min(Basic_Clusters(8,156:256)));
end
for i=1:length(fwd_long)
    PredictFlow(4,fwd_long(i):fwd_long(i)+100)=Basic_Clusters(6,156:256)+abs(min(Basic_Clusters(6,156:256)));
    PredictFlow(5,fwd_long(i):fwd_long(i)+100)=Basic_Clusters(6,156:256)+abs(min(Basic_Clusters(7,156:256)));
    PredictFlow(8,fwd(i):fwd(i)+201)=interp(Basic_Clusters(8,156:256)+abs(min(Basic_Clusters(8,156:256))),2);
end
PredictFlow(4,fwd_three:fwd_three+302)=interp(Basic_Clusters(6,156:256)+abs(min(Basic_Clusters(6,156:256))),3);
PredictFlow(5,fwd_three:fwd_three+302)=interp(Basic_Clusters(7,156:256)+abs(min(Basic_Clusters(7,156:256))),3);    
PredictFlow(8,fwd_three:fwd_three+302)=interp(Basic_Clusters(8,156:256)+abs(min(Basic_Clusters(8,156:256))),3);
for i=1:length(fwd_off)
    PredictFlow(6,fwd_off(i):fwd_off(i)+size(GCaMP6,1)-1)=GCaMP6';
end
clearvars GCaMP6 GCaMP6s back back_off fwd fwd_off back_long back_three fwd_long;


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

idx_rsq_NewFlow=find([model_NewFlow.rsquared]>0.05);
idx_rsq_predict=find([model_predict.rsquared]>0.05);
idx_rsq_basic=find([model_basic.rsquared]>0.05);

[idxKmeans_NewFlow Cmap_NewFlow]=kmeans(ZS2(idx_rsq_NewFlow,:),20,'Options',options,'Distance','cityblock','Replicates',3,'MaxIter',1000,'Display','final');
[idxKmeans_predict Cmap_predict]=kmeans(ZS2(idx_rsq_predict,:),20,'Options',options,'Distance','cityblock','Replicates',3,'MaxIter',1000,'Display','final');
[idxKmeans_basic Cmap_basic]=kmeans(ZS2(idx_rsq_basic,:),20,'Options',options,'Distance','cityblock','Replicates',3,'MaxIter',1000,'Display','final');

figure;
subplot(3,1,1);imagesc(Cmap_NewFlow);
subplot(3,1,2);imagesc(Cmap_basic);
subplot(3,1,3);imagesc(Cmap_predict);

Threshold=0.2;
[Model_NewFlow,GoodBetas_NewFlow]=Test_Regress(Cmap_NewFlow,NewFlow,idxKmeans_NewFlow,Threshold);
[Model_basic,GoodBetas_basic]=Test_Regress(Cmap_basic,NewFlow_basic,idxKmeans_basic,Threshold);
[Model_predict,GoodBetas_predict]=Test_Regress(Cmap_predict,PredictFlow,idxKmeans_predict,Threshold);

ZS_rsq=ZS2(idx_rsq_predict,:);

Corr_BasicCLust=zeros(size(Basic_Clusters,1),size(ZS_rsq,1));
start_complex=[155 1005];start_basic=[155 258];duration=90;

parfor i=1:size(ZS_rsq,1)
    corr_temp=[];
    ZS_rsq_temp=[ZS_rsq(i,start_complex(1):start_complex(1)+duration) ZS_rsq(i,start_complex(2):start_complex(2)+duration)];
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
idx_BasicClust_corr=find(MaxCorr_BasiClust>0.5);
BasicClustData=[];
for basic=1:max(MaxCorr_BasiClust_ind)
    idx_BasicClust_temp=find(MaxCorr_BasiClust_ind(idx_BasicClust_corr)==basic);
    BasicClustData(basic).ZS=ZS_rsq(idx_BasicClust_corr(idx_BasicClust_temp),:);
    BasicClustData(basic).mean=mean(BasicClustData(basic).ZS,1);      
end

Fighandle=figure;
set(Fighandle, 'Position', [100, 100, 1300, 900]);
counter=1;xplot=floor(sqrt(max(MaxCorr_BasiClust_ind)));yplot=ceil(max(MaxCorr_BasiClust_ind)/xplot);
for i=1:max(MaxCorr_BasiClust_ind)    
    subplot(xplot,yplot,counter);plot(x,BasicClustData(i).mean);xlim([0 500]); %hold on;plot(x,(double(Flow_profile2)/10)-1)
    counter=counter+1;
end

