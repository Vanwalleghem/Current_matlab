MatFiles=dir('*matlab*.mat');
name=strcat(MatFiles(1).name);
Calcium=load(name, 'DenoisedTraces');
Calcium=Calcium.DenoisedTraces;
MatFiles(1).number=size(Calcium,1);
Noise=load(name, 'Noise');
Noise=Noise.Noise;
Fitness=load(name, 'idx_components');
Fitness=Fitness.idx_components+1;
GoodCalcium=Calcium(Fitness,1:2500);
GoodNoise=Noise(Fitness,1:2500);
MatFiles(1).GoodNumber=length(Fitness);
for i = 2:length(MatFiles)
    name=strcat(MatFiles(i).name);
    C=load(name, 'DenoisedTraces');
    C=C.DenoisedTraces;
%     if i==3
%         C=[C(:,1) C(:,1) C(:,1:58)];
%     end  
    N=load(name, 'Noise');
    N=N.Noise;
    F=load(name, 'idx_components');
    F=F.idx_components+1;
    GC=C(F,1:2500);
    GN=N(F,1:2500);
    GoodNoise=vertcat(GoodNoise,GN);
    
    GoodCalcium=vertcat(GoodCalcium,GC);

    MatFiles(i).number=size(Calcium,1);
    MatFiles(i).GoodNumber=MatFiles(i-1).GoodNumber+length(F);

end
clearvars GC C S F N name i GS;

ZS2=zscore(GoodCalcium+GoodNoise,1,2);
x = linspace(0.2,size(ZS2,2)/5,size(ZS2,2));y = linspace(1,size(ZS2,1),size(ZS2,1));
Fighandle=figure;
set(Fighandle, 'Position', [100, 100, 1200, 900]);
imagesc(x,y,ZS2(randperm(size(ZS2,1)),:),[-0.5 5]);colormap hot;set(gca,'YTickLabel',[]);

options = statset('UseParallel',1); [idxKmeans Cmap]=kmeans(ZS2,20,'Options',options,'Distance','cityblock','Replicates',5,'MaxIter',1000,'Display','final');

NewFlow=zeros(6,size(ZS2,2));
GCaMP6=[-0.104392135015146,1.69644104899772,5.13796058542217,8.27886020152244,10.3756715204800,11.8173714529814,12.2425184714093,10.8571417354877,8.80831829681196,6.91339112244670,5.46959264663869,4.30868766622567,3.42533619066766,2.75378443486879,2.18017250852183,1.72816235135824,1.32732537295463,1.00684435500268,0.730210038304555,0.530242444093118,0.362253250339685,0.227668255288566,0.0869242416152502,0.000718266708050853,-0.0828334873368325]';
%GCaMP6=[0.000256990000000000;0.00850739000000000;0.0654158300000000;0.0784609000000000;0.0764130100000000;0.0665958600000000;0.0579028900000000;0.0467942900000000;0.0232079800000000;0.0144564400000000;0.00695772000000000;0.00526551000000000;0.00299500000000000;0.00198520000000000;0.00128512000000000;0.00134175000000000;0.000403170000000000;0];
back=    [56 256 556 1006 1106 1466 1826 2086]; %Infuse
back_off=[106 356 706 1056 1206 1516 1926 2176];
fwd=    [156 406 756 856 1256 1316 1526 1626 1986 2236]; %Withdraw
fwd_off=[206 506 806 956 1306 1366 1576 1776 2036 2286] ;
for i=1:length(back)
    NewFlow(1,back(i):back(i)+size(GCaMP6,1)-1)=GCaMP6';
    NewFlow(2,back_off(i):back_off(i)+size(GCaMP6,1)-1)=GCaMP6';
    NewFlow(5,back(i):back_off(i))=1;
end
for i=1:length(fwd)
    NewFlow(3,fwd(i):fwd(i)+size(GCaMP6,1)-1)=GCaMP6';
    NewFlow(4,fwd_off(i):fwd_off(i)+size(GCaMP6,1)-1)=GCaMP6';
    NewFlow(6,fwd(i):fwd_off(i))=1;
end
clearvars GCaMP6 back back_off fwd fwd_off;

[Model_ZS2,GoodBetas_ZS2]=Test_Regress(Cmap_rsq3,NewFlow,idxKmeans,0.4);

PredictFlow=zeros(8,size(ZS2,2));
back=      [56 1006 1466]; %Infuse
back_long= [256 1106 1827 2086];
back_three=557;
back_off=[106 356 706 1056 1206 1516 1926 2176];
back_basic=57;
fwd=    [156 757 1256 1316 1526 1986 2236]; %Withdraw
fwd_long= [407 856];
fwd_three=1626;
fwd_off=[206 506 806 956 1306 1366 1576 1776 2036 2286];
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

Lin_reg_ZS2=[];
parfor i=1:length(ZS2)
    %mdl=stepwiselm(NewFlow',ZS2(i,:),'linear','Criterion','adjrsquared','Upper','linear','Verbose',0);
    mdl=fitlm(NewFlow',ZS2(i,:));    
    Lin_reg_ZS2(i).coef=mdl.Coefficients;
    %Lin_reg_ZS2(i).MSE=mdl.MSE;
    %Lin_reg_ZS2(i).Fitted=mdl.Fitted;
    Lin_reg_ZS2(i).rsquared=mdl.Rsquared.Adjusted; 
    mdl=fitlm(PredictFlow',ZS2(i,:));
    model_ZS_predict(i).coef=mdl.Coefficients;
    model_ZS_predict(i).rsquared=mdl.Rsquared.Adjusted;
end


idx_rsq_ZS=find([Lin_reg_ZS2.rsquared]>0.1);
Fighandle=figure;
set(Fighandle, 'Position', [100, 100, 1400, 900]);
imagesc(ZS(idx_rsq_ZS,:),[0 5]);colormap hot
temp=ZS(idx_rsq_ZS,:);