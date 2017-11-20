MatFiles=dir('*analysis_matlab.mat');
name=strcat(MatFiles(1).name);
Calcium=load(name, 'DenoisedTraces');
Calcium=Calcium.DenoisedTraces;
MatFiles(1).number=size(Calcium,1);
Spikes=load(name, 'Spikes');
Spikes=Spikes.Spikes;
Noise=load(name, 'Noise');
Noise=Noise.Noise;
Fitness=load(name, 'idx_components');
Fitness=Fitness.idx_components+1;
GoodCalcium=Calcium(Fitness,:);
GoodSpikes=Spikes(Fitness,:);
MatFiles(1).GoodNumber=length(Fitness);
MatFiles(1).GC=GoodCalcium;
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
    GC=C(F,:);
    GS=S(F,:);
    Noise=vertcat(Noise,N);
    Calcium=vertcat(Calcium,C);
    Spikes=vertcat(Spikes,S);
    Fitness=horzcat(Fitness,F);
    GoodCalcium=vertcat(GoodCalcium,GC);
    GoodSpikes=vertcat(GoodSpikes,GS);
    MatFiles(i).number=size(Calcium,1);
    MatFiles(i).GoodNumber=MatFiles(i-1).GoodNumber+length(F);
    MatFiles(i).GC=GC;
end
clearvars GC C S F N name i GS;

ZS=zscore(GoodCalcium,1,2);
x = linspace(0.2,size(ZS,2)/5,size(ZS,2));y = linspace(1,size(ZS,1),size(ZS,1));
Fighandle=figure;
set(Fighandle, 'Position', [100, 100, 1200, 900]);
imagesc(x,y,ZS,[0 5]);colormap hot;set(gca,'YTickLabel',[]);

Numbers=[1 [MatFiles.GoodNumber]];
counter=1;
idx_Plane=nan(length(GoodCalcium),1);
idx_Fish=nan(length(GoodCalcium),1);
name=strcat(MatFiles(1).name);
[Plane,~]=regexp(name,'\d+_(\d+)_','tokens','match');Plane=str2num(Plane{1}{1});
% [Plane,~]=regexp(name,'\d\D(\d+)um','tokens','match');Plane=str2num(Plane{1}{1});
[Fish,~]=regexp(name,'(\d+)_\d+_','tokens','match');Fish=str2num(Fish{1}{1});
% [Fish,~]=regexp(name,'(\d)\D\d+um','tokens','match');Fish=str2num(Fish{1}{1});
% idx_Plane(1:Numbers(2))=Plane;
% idx_Fish(1:Numbers(2))=Fish;
for i=1:length(MatFiles)
	%[Fish,~]=regexp(files{i},'(\d+)_','tokens','match');Fish=str2num(Fish{1}{1});
    name=strcat(MatFiles(i).name);
    if findstr(MatFiles(i).name,'2planes')
        [Plane,~]=regexp(name,'f\d-(\d+)um_','tokens','match');Plane=str2num(Plane{1}{1});
        [Fish,~]=regexp(name,'f(\d)-\d+um_','tokens','match');Fish=str2num(Fish{1}{1});
    else
        [Plane,~]=regexp(name,'\d+_(\d+)_','tokens','match');Plane=str2num(Plane{1}{1});
        [Fish,~]=regexp(name,'(\d+)_\d+_','tokens','match');Fish=str2num(Fish{1}{1});
    end
    %[Plane,~]=regexp(name,'\d+_(\d+)_','tokens','match');Plane=str2num(Plane{1}{1});
    %[Fish,~]=regexp(name,'(\d+)_\d+_','tokens','match');Fish=str2num(Fish{1}{1});
   
    idx_Plane(Numbers(i):Numbers(i+1))=Plane;
    idx_Fish(Numbers(i):Numbers(i+1))=Fish;
end
clearvars i Fish Plane name counter

flow=zeros(6,700);
GCaMP6=[-0.104392135015146,1.69644104899772,5.13796058542217,8.27886020152244,10.3756715204800,11.8173714529814,12.2425184714093,10.8571417354877,8.80831829681196,6.91339112244670,5.46959264663869,4.30868766622567,3.42533619066766,2.75378443486879,2.18017250852183,1.72816235135824,1.32732537295463,1.00684435500268,0.730210038304555,0.530242444093118,0.362253250339685,0.227668255288566,0.0869242416152502,0.000718266708050853,-0.0828334873368325]';
%GCaMP6=[0.000256990000000000;0.00850739000000000;0.0654158300000000;0.0784609000000000;0.0764130100000000;0.0665958600000000;0.0579028900000000;0.0467942900000000;0.0232079800000000;0.0144564400000000;0.00695772000000000;0.00526551000000000;0.00299500000000000;0.00198520000000000;0.00128512000000000;0.00134175000000000;0.000403170000000000;0];
back=[57 257 457];
back_off=[106 306 506];
fwd=[157 357 557];
fwd_off=[207 407 607];
flow(1,back(1):back(1)+size(GCaMP6,1)-1)=GCaMP6';
flow(1,back(2):back(2)+size(GCaMP6,1)-1)=GCaMP6';
flow(1,back(3):back(3)+size(GCaMP6,1)-1)=GCaMP6';
flow(2,back_off(1):back_off(1)+size(GCaMP6,1)-1)=GCaMP6';
flow(2,back_off(2):back_off(2)+size(GCaMP6,1)-1)=GCaMP6';
flow(2,back_off(3):back_off(3)+size(GCaMP6,1)-1)=GCaMP6';
flow(3,fwd(1):fwd(1)+size(GCaMP6,1)-1)=GCaMP6';
flow(3,fwd(2):fwd(2)+size(GCaMP6,1)-1)=GCaMP6';
flow(3,fwd(3):fwd(3)+size(GCaMP6,1)-1)=GCaMP6';
flow(4,fwd_off(1):fwd_off(1)+size(GCaMP6,1)-1)=GCaMP6';
flow(4,fwd_off(2):fwd_off(2)+size(GCaMP6,1)-1)=GCaMP6';
flow(4,fwd_off(3):fwd_off(3)+size(GCaMP6,1)-1)=GCaMP6';
flow(5,back(1):back(1)+43)=1;
flow(5,back(2):back(2)+43)=1;
flow(5,back(3):back(3)+43)=1;
flow(6,fwd(1):fwd(1)+43)=1;
flow(6,fwd(2):fwd(2)+43)=1;
flow(6,fwd(3):fwd(3)+43)=1;
clearvars GCaMP6 back back_off fwd fwd_off;

flow=repmat(flow,3,3);
flow(1:6,701:end)=0;
flow(7:12,1401:end)=0;flow(7:12,1:701)=0;
flow(13:end,1:1401)=0;

Model_LR=[];
parfor i=1:size(ZS,1)
    mdl=stepwiselm(flow',ZS(i,:),'Upper','linear','Intercept',false,'Criterion','bic','verbose',0);
    Model_LR(i).coef=mdl.Coefficients;
    Model_LR(i).MSE=mdl.MSE;
    %Model_LR(i).Fitted=mdl.Fitted;
    Model_LR(i).rsquared=mdl.Rsquared.Adjusted;
end

idx_rsq_sort=find([Model_LR.rsquared]>0.1 & [Model_LR.rsquared]<1);
ZS_sort=ZS(idx_rsq_sort,:);
figure;imagesc(ZS_sort,[-1 4]);colormap hot