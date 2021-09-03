%of 4 (may need to interpolate it or make a new one.
spike=[0,1.69644104899772,5.13796058542217,8.27886020152244,10.3756715204800,11.8173714529814,12.2425184714093,10.8571417354877,8.80831829681196,6.91339112244670,5.46959264663869,4.30868766622567,3.42533619066766,2.75378443486879,2.18017250852183,1.72816235135824,1.32732537295463,1.00684435500268,0.730210038304555,0.530242444093118,0.362253250339685,0.227668255288566,0.0869242416152502,0.000718266708050853,-0.0828334873368325]';
spike=spike/max(spike); %Normalize the spike so the regression coefficients are easier to "read"

test=F(iscell(:,2)>0.8,:);
DF=DeltaF2(test,50,11);

Y = datasample(DF,500);
SpikeAvg=[];
for i=1:size(Y,1)
    [~,temp]=findpeaks(Y(i,:),'MinPeakProminence',0.05,'MinPeakDistance',20);
    SpikeTemp=zeros(length(temp),13);
    temp(temp<5)=[];
    temp(temp>1450)=[];
    for frame=1:length(temp)
        SpikeTemp(frame,:)=Y(i,temp(frame)-3:temp(frame)+9)/max(Y(i,temp(frame)-3:temp(frame)+9));
    end
    SpikeAvg=cat(1,SpikeAvg,SpikeTemp);
end


figure;plot(spike(1:2:end))
test=mean(SpikeAvg,1)-min(mean(SpikeAvg,1));
test=test/max(test);
hold on;plot(test);


