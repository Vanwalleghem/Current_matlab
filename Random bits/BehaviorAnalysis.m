RawData=zeros(5,60000);
for i=1:5
    RawData(i,:) = csvread(strcat('D:\Pictures\processed\Flow\Behavior\Fish',num2str(i+1),'_backSub.csv'),1,1);
end


% AlignedData=zscore(RawData,1,2);
% for i=[1 3 4 5]
%     [temp,temp2]=alignsignals(AlignedData(2,:),AlignedData(i,:),200);
%     AlignedData(i,:)=temp2(1:60000);
% end
load('_Complex_finalish_v2.mat', 'Speed_flow')
temp=timeseries(sum(Speed_flow,1)',[0.2:0.2:500]);
temp2=resample(temp3,[0.2:0.01:500]);

temp={};
temp{1}=resample(timeseries(Speed_flow(1,:)',[0.2:0.2:500]),[0.2:0.01:500]);
temp{2}=resample(timeseries(Speed_flow(2,:)',[0.2:0.2:500]),[0.2:0.01:500]);

Flow_resampled=temp{1}.data-temp{2}.data;
Flow_resampled=[zeros(1020,1); Flow_resampled];
Flow_resampled(length(Flow_resampled):60000)=0;

Flow=zeros(2,60000);
Flow(1,:)=csvread(strcat('D:\Pictures\processed\Flow\Behavior\Flow_Fish3.csv'),1,1);
Flow(2,:)=csvread(strcat('D:\Pictures\processed\Flow\Behavior\Flow_Fish5.csv'),1,1);
Flow(1,:)=zscore(lowpass(Flow(1,:),5,100),1,2);
Flow(2,:)=zscore(lowpass(Flow(2,:),5,100),1,2);

temp2=[zeros(1,1020) temp2.data'];
temp2(length(temp2):60000)=0;

figure;plot(zscore(Flow,1,2)');hold on;plot(temp2);

AlignedData=zscore(RawData,1,2);
for i=1:5
    [temp,temp2]=alignsignals(Flow(2,:),AlignedData(i,:),200);
    AlignedData(i,:)=temp2(1:60000);
end

HighPassData=AlignedData;
LowPassData=AlignedData;
for i=1:5
    HighPassData(i,:)=highpass(AlignedData(i,:),5,100);
    LowPassData(i,:)=lowpass(AlignedData(i,:),5,100);
end

Baseline_mvt={};
Speed_1_mvt={};
Speed_2_mvt={};
for i=1:5
    [Speed_1_mvt{i,1},Speed_1_mvt{i,2}]=findpeaks(HighPassData(i,temp2==1),'MinPeakProminence',0.5,'MinPeakDistance',5);
    [Speed_2_mvt{i,1},Speed_2_mvt{i,2}]=findpeaks(HighPassData(i,temp2==2),'MinPeakProminence',0.5,'MinPeakDistance',5);
    [Baseline_mvt{i,1},Baseline_mvt{i,2}]=findpeaks(HighPassData(i,temp2==0),'MinPeakProminence',0.5,'MinPeakDistance',5);
end

[sum(temp2==1) sum(temp2==2) sum(temp2==0)]

PrismTemp=zeros(5,3);
temp=cellfun(@length,Speed_1_mvt);PrismTemp(:,2)=temp(:,1)/(sum(temp2==1) / 100);
temp=cellfun(@length,Speed_2_mvt);PrismTemp(:,3)=temp(:,1)/(sum(temp2==2) / 100);
temp=cellfun(@length,Baseline_mvt);PrismTemp(:,1)=temp(:,1)/(sum(temp2==0) / 100);

tail_flicks_perS=zeros(5,3000);
for i=1:5    
    [~,temp]=findpeaks(HighPassData(i,:),'MinPeakProminence',0.5,'MinPeakDistance',5);
	tail_flicks_perS(i,:)=histcounts(temp,[0:20:60000]);
end
figure;plot(tail_flicks_perS');hold on;plot([zeros(1,1020/20) sum(Speed_flow,1)]);

Baseline_mvt_complex={};
Speed_1_mvt_forward={};
Speed_2_mvt_forward={};
Speed_1_mvt_backward={};
Speed_2_mvt_backward={};
for i=1:5
    [Speed_1_mvt_forward{i,1},Speed_1_mvt_forward{i,2}]=findpeaks(HighPassData(i,Flow_resampled==1),'MinPeakProminence',0.5,'MinPeakDistance',5);
    [Speed_2_mvt_forward{i,1},Speed_2_mvt_forward{i,2}]=findpeaks(HighPassData(i,Flow_resampled==2),'MinPeakProminence',0.5,'MinPeakDistance',5);
    [Baseline_mvt_complex{i,1},Baseline_mvt_complex{i,2}]=findpeaks(HighPassData(i,Flow_resampled==0),'MinPeakProminence',0.5,'MinPeakDistance',5);
    [Speed_1_mvt_backward{i,1},Speed_1_mvt_backward{i,2}]=findpeaks(HighPassData(i,Flow_resampled==-1),'MinPeakProminence',0.5,'MinPeakDistance',5);
    [Speed_2_mvt_backward{i,1},Speed_2_mvt_backward{i,2}]=findpeaks(HighPassData(i,Flow_resampled==-2),'MinPeakProminence',0.5,'MinPeakDistance',5);
end


PrismTemp=zeros(5,5);
temp=cellfun(@length,Speed_1_mvt_forward);PrismTemp(:,2)=temp(:,1)/(sum(Flow_resampled==1) / 100);
temp=cellfun(@length,Speed_2_mvt_forward);PrismTemp(:,3)=temp(:,1)/(sum(Flow_resampled==2) / 100);
temp=cellfun(@length,Baseline_mvt_complex);PrismTemp(:,1)=temp(:,1)/(sum(Flow_resampled==0) / 100);
temp=cellfun(@length,Speed_1_mvt_backward);PrismTemp(:,4)=temp(:,1)/(sum(Flow_resampled==-1) / 100);
temp=cellfun(@length,Speed_2_mvt_backward);PrismTemp(:,5)=temp(:,1)/(sum(Flow_resampled==-2) / 100);

Fighandle=figure;
set(Fighandle, 'Position', [100, 100, 1300, 900]);
ha=tight_subplot(5,1);
for i=1:5
    axes(ha(i));
    plot(tail_flicks_perS(i,:));hold on;plot([zeros(2,1020/20) Speed_flow]');
end

Flow_resampled_ind=zeros(3,size(Flow_resampled,1));
Flow_resampled_ind(1,Flow_resampled==1)=1;
Flow_resampled_ind(2,Flow_resampled==2)=1;
Flow_resampled_ind(3,Flow_resampled==0)=1;
Flow_resampled_ind(4,Flow_resampled==-1)=1;
Flow_resampled_ind(5,Flow_resampled==-2)=1;

Corr_Motor_flow=[];
for i=1:5
    Corr_Motor_flow(i,:)=pdist2(HighPassData(i,:),Flow_resampled_ind,'correlation');
end
Corr_Motor_flow=1-Corr_Motor_flow;
