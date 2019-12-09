%of 4 (may need to interpolate it or make a new one.
spike=[0,1.69644104899772,5.13796058542217,8.27886020152244,10.3756715204800,11.8173714529814,12.2425184714093,10.8571417354877,8.80831829681196,6.91339112244670,5.46959264663869,4.30868766622567,3.42533619066766,2.75378443486879,2.18017250852183,1.72816235135824,1.32732537295463,1.00684435500268,0.730210038304555,0.530242444093118,0.362253250339685,0.227668255288566,0.0869242416152502,0.000718266708050853,-0.0828334873368325]';
spike=spike/max(spike); %Normalize the spike so the regression coefficients are easier to "read"
framerate=2;

NewFlow=zeros(4,size(ZS,2));
fwd=[0 40 100 190 210 282 374 404]; %Infuse
fwd_off=[10 60 130 200 230 292 394 424];
back=[20 70 140 160 240 252 294 314 394 434]; %Withdraw
back_off=[30 90 150 180 250 262 304 344 404 444];
framerate=2.1;start=50;
back= round(start*framerate+ back*framerate); %Infuse
back_off=round(start*framerate+ back_off*framerate);
fwd=  round(start*framerate+ fwd*framerate); %Withdraw
fwd_off= round(start*framerate+ fwd_off*framerate);
for i=1:length(back)
NewFlow(1,back(i):back(i)+size(spike,1)-1)=spike';
NewFlow(2,back(i):back_off(i))=1;
end
for i=1:length(fwd)
NewFlow(3,fwd(i):fwd(i)+size(spike,1)-1)=spike';
NewFlow(4,fwd(i):fwd_off(i))=1;
end
clearvars GCaMP6 back back_off fwd fwd_off;

StimLength=1200;
x = linspace(0,StimLength/framerate,StimLength);


Stimuli=zeros(4,size(ZS,2));
fwd2=     [0  40 100 192 208 290 364 409]; %Infuse
fwd_off2=[10 60 134 202 228 300 384 429];
back2=     [22 69 142 165 250 262 294 314 384 430]; %Withdraw
back_off2=[32 89 152 185 260 272 304 344 394 440];
framerate=2.1;start=50;
back2= round(start*framerate+ back2*framerate); %Infuse
back_off2=round(start*framerate+ back_off2*framerate);
fwd2=  round(start*framerate+ fwd2*framerate); %Withdraw
fwd_off2= round(start*framerate+ fwd_off2*framerate);
for i=1:length(back2)
Stimuli(1,back2(i):back2(i)+size(spike,1)-1)=spike';
Stimuli(2,back_off2(i):back_off2(i)+size(spike,1)-1)=spike';
Stimuli(3,back2(i):back_off2(i))=1;
end
for i=1:length(fwd2)
Stimuli(4,fwd2(i):fwd2(i)+size(spike,1)-1)=spike';
Stimuli(5,fwd_off2(i):fwd_off2(i)+size(spike,1)-1)=spike';
Stimuli(6,fwd2(i):fwd_off2(i))=1;
end
clearvars GCaMP6 back back_off fwd fwd_off;
Stimuli(7,:)=linspace(0,1,size(Stimuli,2));


%Our basic analysis starts with a Linear Regression
%You can choose to only use part of the timeseries
LinearModelStimuli=[];
parfor i=1:size(ZS,1)    
    mdl=fitlm(Stimuli',ZS(i,:));
    LinearModelStimuli(i).coef=mdl.Coefficients;    
    LinearModelStimuli(i).rsquared=mdl.Rsquared.Adjusted;
end

rsq_all=[LinearModelStimuli.rsquared];
figure;histogram(rsq_all);
ThresholdR2=0.1
idx_rsq=find([LinearModelStimuli.rsquared]>ThresholdR2);
ZS_rsq=ZS(idx_rsq,:);

%The alternative/follow up on linear regression is clustering, we usually
%use k-means, before or after filtering, as you prefer
options = statset('UseParallel',1); %parallelize the replicates
[idxKmeans_ZS Cmap_ZS]=kmeans(ZS_rsq,20,'Options',options,'Distance','correlation','Replicates',3,'MaxIter',1000,'Display','final');
figure;
imagesc(Cmap_ZS);
