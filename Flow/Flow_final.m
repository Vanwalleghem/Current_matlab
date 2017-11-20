figure;
for i=6:8
    plot(Basic_Clusters(i,156:256));
    pause;
end

NewFlow3=zeros(8,size(ZS2,2));
GCaMP6=[-0.104392135015146,1.69644104899772,5.13796058542217,8.27886020152244,10.3756715204800,11.8173714529814,12.2425184714093,10.8571417354877,8.80831829681196,6.91339112244670,5.46959264663869,4.30868766622567,3.42533619066766,2.75378443486879,2.18017250852183,1.72816235135824,1.32732537295463,1.00684435500268,0.730210038304555,0.530242444093118,0.362253250339685,0.227668255288566,0.0869242416152502,0.000718266708050853,-0.0828334873368325]';
%GCaMP6s=[0.000256990000000000;0.00850739000000000;0.0654158300000000;0.0784609000000000;0.0764130100000000;0.0665958600000000;0.0579028900000000;0.0467942900000000;0.0232079800000000;0.0144564400000000;0.00695772000000000;0.00526551000000000;0.00299500000000000;0.00198520000000000;0.00128512000000000;0.00134175000000000;0.000403170000000000;0];
GCaMP6=[-0.104392135015146,1.69644104899772,5.13796058542217,8.27886020152244,10.3756715204800,11.8173714529814,12.2425184714093,10.8571417354877,8.80831829681196,6.91339112244670,5.46959264663869,4.30868766622567,3.42533619066766,2.75378443486879,2.18017250852183,1.72816235135824,1.32732537295463,1.00684435500268,0.730210038304555,0.530242444093118,0.362253250339685,0.227668255288566,0.0869242416152502,0.000718266708050853,-0.0828334873368325]';
GCaMP6s=interp(GCaMP6,2);
back=    [56 256 557 1006 1106 1466 1827 2086]-10; %Infuse
back_off=[106 356 706 1056 1206 1516 1926 2176]-10;
fwd=    [156 407 757 856 1256 1316 1526 1626 1986 2236]-10; %Withdraw
fwd_off=[206 506 806 956 1306 1366 1576 1776 2036 2286]-10;
for i=1:length(back)
    NewFlow3(1,back(i):back(i)+100)=Basic_Clusters(6,156:256)+abs(min(Basic_Clusters(6,156:256)));
    NewFlow3(2,back(i):back(i)+100)=Basic_Clusters(7,156:256)+abs(min(Basic_Clusters(7,156:256)));
    NewFlow3(3,back_off(i):back_off(i)+size(GCaMP6,1)-1)=GCaMP6';
    NewFlow3(7,back(i):back(i)+100)=Basic_Clusters(8,156:256)+abs(min(Basic_Clusters(8,156:256)));
end
for i=1:length(fwd)
    NewFlow3(4,fwd(i):fwd(i)+100)=Basic_Clusters(6,156:256)+abs(min(Basic_Clusters(6,156:256)));
    NewFlow3(5,fwd(i):fwd(i)+100)=Basic_Clusters(7,156:256)+abs(min(Basic_Clusters(7,156:256)));
    NewFlow3(6,fwd_off(i):fwd_off(i)+size(GCaMP6,1)-1)=GCaMP6';
    NewFlow3(8,fwd(i):fwd(i)+100)=Basic_Clusters(8,156:256)+abs(min(Basic_Clusters(8,156:256)));
end
clearvars GCaMP6 GCaMP6s back back_off fwd fwd_off;

parfor i=1:size(ZS2,1)
    mdl=stepwiselm(NewFlow3',ZS2(i,:),'Upper','linear','Intercept',false,'Criterion','bic','verbose',0);
    model_ZS_detrend(i).coef=mdl.Coefficients;
    model_ZS_detrend(i).MSE=mdl.MSE;
    model_ZS_detrend(i).Fitted=mdl.Fitted;
    model_ZS_detrend(i).rsquared=mdl.Rsquared.Adjusted;
end
clearvars slow_fast slow_fast_fwd i

options = statset('UseParallel',1); [idxKmeans_ZS2b Cmap_ZS2b]=kmeans(ZS2_rsq,30,'Options',options,'Distance','cityblock','Replicates',10,'MaxIter',2000,'Display','final');
[Model_ZS2b,GoodBetas_ZS2b]=Test_Regress(Cmap_ZS2b,NewFlow3,idxKmeans_ZS2b,0.2);

coefficients={};
for idx=1:length(model_ZS_detrend)
    coef=[model_ZS_detrend(idx).coef];
    temp=coef.Properties.RowNames;temp=regexp(temp,'x(\d+)','tokens');
    if ~isempty(temp)
        %temp=[temp{:}];temp=[temp{:}];temp=[temp{:}];%temp=str2num(temp);
        for coef_idx=1:height(coef)
            if coef.pValue(coef_idx)<0.05
                coefficients{idx,str2num(temp{coef_idx}{1}{1})}=coef.Estimate(coef_idx);
            end
        end
    end
end
idxempty=cellfun('isempty',coefficients);
coefficients(idxempty)={0};
clearvars idxempty idx coef_idx coef
coefficients=cell2mat(coefficients);
coefficients_all=coefficients;
coefficients=coefficients(idx_rsq_ZS2,:);

options = statset('UseParallel',1); [idxKmeans_coef Cmap_coef]=kmeans(coefficients,10,'Options',options,'Distance','cityblock','Replicates',10,'MaxIter',2000,'Display','final');

temp=zeros(10,size(ZS2,2));
for i=1:10
   temp(i,:)=mean(ZS2_rsq(find(idxKmeans_coef==i),:),1);
end

figure;
for i=1:10
   plot(mean(ZS2_rsq(find(idxKmeans_coef==i),:),1));hold on;plot((Flow_profile2/10)-1);hold off;pause;
end

test=ZS2_rsq;
test(ZS2_rsq<0.1)=0;
options = statset('UseParallel',1); [idxKmeans_test Cmap_test]=kmeans(test,30,'Options',options,'Distance','cityblock','Replicates',10,'MaxIter',2000,'Display','final');
[Model_test,GoodBetas_test]=Test_Regress(Cmap_test,NewFlow3,idxKmeans_test,0.2);

Spikes=zeros(size(ZS2));
parfor i=1:size(ZS2,1)
    [~, Spikes(i,:),~]=deconvolveCa(ZS2(i,:));
end

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
    %PredictFlow(1,fwd_long(i):fwd_long(i)+201)=interp(Basic_Clusters(6,156:256)+abs(min(Basic_Clusters(6,156:256))),2);
    %PredictFlow(2,fwd_long(i):fwd_long(i)+201)=interp(Basic_Clusters(7,156:256)+abs(min(Basic_Clusters(7,156:256))),2);    
    %PredictFlow(7,fwd_long(i):fwd_long(i)+201)=interp(Basic_Clusters(8,156:256)+abs(min(Basic_Clusters(8,156:256))),2);
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
    mdl=stepwiselm(PredictFlow',ZS2(i,:),'Upper','linear','Intercept',false,'Criterion','bic','verbose',0);
    model_ZS_pretend(i).coef=mdl.Coefficients;
    model_ZS_pretend(i).MSE=mdl.MSE;
    model_ZS_pretend(i).Fitted=mdl.Fitted;
    model_ZS_pretend(i).rsquared=mdl.Rsquared.Adjusted;
end
clearvars slow_fast slow_fast_fwd i

idx_rsq_predict=find([model_ZS_pretend.rsquared]>0.1);
Fighandle=figure;
set(Fighandle, 'Position', [100, 100, 1400, 900]);
imagesc(ZS2(idx_rsq_predict,:),[0 5]);colormap hot
ZS2_predict=ZS2(idx_rsq_predict,:);

options = statset('UseParallel',1); [idxKmeans_ZS2_rsq Cmap_ZS2_rsq]=kmeans(ZS2_predict,30,'Options',options,'Distance','cityblock','Replicates',10,'MaxIter',2000,'Display','final');

MatFiles_short=dir('*matlab*.mat');
name=strcat(MatFiles_short(1).name);
Calcium_short=load(name, 'DenoisedTraces');
Calcium_short=Calcium_short.DenoisedTraces;
MatFiles_short(1).number=size(Calcium_short,1);
Spikes_short=load(name, 'Spikes');
Spikes_short=Spikes_short.Spikes;
Noise_short=load(name, 'Noise');
Noise_short=Noise_short.Noise;
Fitness=load(name, 'idx_components');
Fitness=Fitness.idx_components+1;
GoodCalcium_short=Calcium_short(Fitness,:);
GoodSpikes_short=Spikes_short(Fitness,:);
GoodNoise_short=Noise_short(Fitness,:);
MatFiles_short(1).GoodNumber=length(Fitness);
MatFiles_short(1).GC=GoodCalcium_short;
for i = 2:length(MatFiles_short)
    name=strcat(MatFiles_short(i).name);
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
    GC=C(F,:);
    GS=S(F,:);
    GN=N(F,:);
    Noise_short=vertcat(Noise_short,N);
    Calcium_short=vertcat(Calcium_short,C);
    Spikes_short=vertcat(Spikes_short,S);
    Fitness=horzcat(Fitness,F);
    GoodCalcium_short=vertcat(GoodCalcium_short,GC);
    GoodSpikes_short=vertcat(GoodSpikes_short,GS);
    GoodNoise_short=vertcat(GoodNoise_short,GN);
    MatFiles_short(i).number=size(Calcium_short,1);
    MatFiles_short(i).GoodNumber=MatFiles_short(i-1).GoodNumber+length(F);
    MatFiles_short(i).GC=GC;
end
clearvars GC C S F N name i GS GN;
ZS_short=zscore(GoodCalcium_short,1,2);
ZS2_short=detrend(ZS_short(:,20:end)')';

parfor i=1:size(ZS2_short,1)
    mdl=stepwiselm(PredictFlow(:,1:1981)',ZS2_short(i,:),'Upper','linear','Intercept',false,'Criterion','bic','verbose',0);
    model_ZS_pretend_short(i).coef=mdl.Coefficients;
    model_ZS_pretend_short(i).MSE=mdl.MSE;
    model_ZS_pretend_short(i).Fitted=mdl.Fitted;
    model_ZS_pretend_short(i).rsquared=mdl.Rsquared.Adjusted;
end
clearvars slow_fast slow_fast_fwd i
idx_rsq_short=find([model_ZS_pretend_short.rsquared]>0.05);
ZS2_short_rsq=ZS2_short(idx_rsq_short,:);

coefficients_predict={};
for idx=1:length(model_ZS_pretend)
    coef=[model_ZS_pretend(idx).coef];
    temp=coef.Properties.RowNames;temp=regexp(temp,'x(\d+)','tokens');
    if ~isempty(temp)
        %temp=[temp{:}];temp=[temp{:}];temp=[temp{:}];%temp=str2num(temp);
        for coef_idx=1:height(coef)
            if coef.pValue(coef_idx)<0.05
                coefficients_predict{idx,str2num(temp{coef_idx}{1}{1})}=coef.Estimate(coef_idx);
            end
        end
    end
end
idxempty=cellfun('isempty',coefficients_predict);
coefficients_predict(idxempty)={0};
clearvars idxempty idx coef_idx coef
coefficients_predict=cell2mat(coefficients_predict);
coefficients_all_pred=coefficients_predict;
coefficients_predict=coefficients_all_pred(idx_rsq_predict,:);

[~, max_coef]=max(coefficients_predict,[],2);

MeanByCoef=zeros(max(max_coef),size(ZS2,2));
for i=1:max(max_coef)
    MeanByCoef(i,:)=mean(ZS2_predict(find(max_coef==i),:),1);
end

figure;
for i=1:max(max_coef)
    plot(MeanByCoef(i,:));pause
end


PredictFlow=zeros(8,size(ZS2,2));
back=      [56 1006 1466]-10; %Infuse
back_long= [256 1106 1827 2086]-10;
back_three=557-10;
back_off=[106 356 706 1056 1206 1516 1926 2176]-10;
back_basic=57;
fwd=    [156 757 1256 1316 1526 1986 2236]-10; %Withdraw
fwd_long= [407 856]-10;
fwd_three=1626-10;
fwd_off=[206 506 806 956 1306 1366 1576 1776 2036 2286]-10;
fwd_basic=157;
for i=1:length(back)
    PredictFlow(1,back(i):back(i)+100)=Basic_Clusters(3,back_basic:back_basic+100)+abs(min(Basic_Clusters(3,back_basic:back_basic+100)));
    PredictFlow(2,back(i):back(i)+100)=Basic_Clusters(4,back_basic:back_basic+100)+abs(min(Basic_Clusters(4,back_basic:back_basic+100)));    
    PredictFlow(3,back(i):back(i)+100)=Basic_Clusters(5,back_basic:back_basic+100)+abs(min(Basic_Clusters(5,back_basic:back_basic+100)));
    PredictFlow(4,back(i):back(i)+100)=Basic_Clusters(8,fwd_basic:fwd_basic+100)+abs(min(Basic_Clusters(8,fwd_basic:fwd_basic+100)));
end
for i=1:length(back_long)
    PredictFlow(1,back_long(i):back_long(i)+201)=interp(Basic_Clusters(3,back_basic:back_basic+100)+abs(min(Basic_Clusters(3,back_basic:back_basic+100))),2);
    PredictFlow(2,back_long(i):back_long(i)+201)=interp(Basic_Clusters(4,back_basic:back_basic+100)+abs(min(Basic_Clusters(4,back_basic:back_basic+100))),2);    
    PredictFlow(3,back_long(i):back_long(i)+201)=interp(Basic_Clusters(5,back_basic:back_basic+100)+abs(min(Basic_Clusters(5,back_basic:back_basic+100))),2);
    PredictFlow(4,back_long(i):back_long(i)+201)=interp(Basic_Clusters(8,fwd_basic:fwd_basic+100)+abs(min(Basic_Clusters(8,fwd_basic:fwd_basic+100))),2);
end
PredictFlow(1,back_three:back_three+302)=interp(Basic_Clusters(3,back_basic:back_basic+100)+abs(min(Basic_Clusters(3,back_basic:back_basic+100))),3);
PredictFlow(2,back_three:back_three+302)=interp(Basic_Clusters(4,back_basic:back_basic+100)+abs(min(Basic_Clusters(4,back_basic:back_basic+100))),3);    
PredictFlow(3,back_three:back_three+302)=interp(Basic_Clusters(5,back_basic:back_basic+100)+abs(min(Basic_Clusters(5,back_basic:back_basic+100))),3);
PredictFlow(4,back_three:back_three+302)=interp(Basic_Clusters(8,fwd_basic:fwd_basic+100)+abs(min(Basic_Clusters(8,fwd_basic:fwd_basic+100))),3);
for i=1:length(fwd)
    PredictFlow(5,fwd(i):fwd(i)+100)=Basic_Clusters(6,156:256)+abs(min(Basic_Clusters(6,156:256)));
    PredictFlow(6,fwd(i):fwd(i)+100)=Basic_Clusters(7,156:256)+abs(min(Basic_Clusters(7,156:256)));    
    PredictFlow(7,fwd(i):fwd(i)+100)=Basic_Clusters(8,156:256)+abs(min(Basic_Clusters(8,156:256)));
    PredictFlow(8,fwd(i):fwd(i)+100)=Basic_Clusters(5,back_basic:back_basic+100)+abs(min(Basic_Clusters(5,back_basic:back_basic+100)));
end
for i=1:length(fwd_long)
    PredictFlow(5,fwd_long(i):fwd_long(i)+201)=interp(Basic_Clusters(6,156:256)+abs(min(Basic_Clusters(6,156:256))),2);
    PredictFlow(6,fwd_long(i):fwd_long(i)+201)=interp(Basic_Clusters(7,156:256)+abs(min(Basic_Clusters(7,156:256))),2);
    PredictFlow(7,fwd_long(i):fwd_long(i)+201)=interp(Basic_Clusters(8,156:256)+abs(min(Basic_Clusters(8,156:256))),2);
    PredictFlow(8,fwd_long(i):fwd_long(i)+201)=interp(Basic_Clusters(5,back_basic:back_basic+100)+abs(min(Basic_Clusters(5,back_basic:back_basic+100))),2);
end
PredictFlow(5,fwd_three:fwd_three+302)=interp(Basic_Clusters(6,156:256)+abs(min(Basic_Clusters(6,156:256))),3);
PredictFlow(6,fwd_three:fwd_three+302)=interp(Basic_Clusters(7,156:256)+abs(min(Basic_Clusters(7,156:256))),3);    
PredictFlow(7,fwd_three:fwd_three+302)=interp(Basic_Clusters(8,156:256)+abs(min(Basic_Clusters(8,156:256))),3);
PredictFlow(8,fwd_three:fwd_three+302)=interp(Basic_Clusters(5,back_basic:back_basic+100)+abs(min(Basic_Clusters(5,back_basic:back_basic+100))),3);

clearvars GCaMP6 GCaMP6s back back_off fwd fwd_off back_long back_three fwd_long;

parfor i=1:size(ZS2,1)
    mdl=stepwiselm(PredictFlow',ZS2(i,:),'Upper','linear','Intercept',true,'Criterion','bic','verbose',0);
    model_ZS_predict(i).coef=mdl.Coefficients;
    model_ZS_predict(i).rsquared=mdl.Rsquared.Adjusted;
end
clearvars slow_fast slow_fast_fwd i
idx_rsq_ZS2b=find([model_ZS_predict.rsquared]>0.05);
figure;imagesc(ZS2(idx_rsq_ZS2b,:), [-0.5 4]);colormap hot;

figure;
for i=1:30
    plot(Cmap_ZS2_pred(i,:));hold on;plot(NewFlow3);pause
end