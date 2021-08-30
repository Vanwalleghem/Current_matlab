Onset_stack=[ZS2(idx_biOnset(randperm(length(idx_biOnset))),:);
ZS2(idx_bwdOnset(randperm(length(idx_bwdOnset))),:);
ZS2(idx_fwdOnset(randperm(length(idx_fwdOnset))),:)];

figure;imagesc(detrend(Onset_stack')',[-0.5 5]);colormap hot

BasicOnset=struct();
Threshold=0.4;
for i=1:3
    switch i
        case 1
            idx_temp=idx_biOnset;
            i
        case 2
            idx_temp=idx_bwdOnset;
            i
        case 3
            idx_temp=idx_fwdOnset;
            i
    end
    mean_temp=mean(ZS2((idx_temp),:),1);
    ZS_temp=ZS2((idx_temp),1:end);
    corr_temp=zeros(size(ZS_temp,1),1);
    parfor i=1:length(idx_temp)
        temp=corrcoef(mean_temp, ZS_temp(i,:));
        corr_temp(i)=temp(1,2);
    end
    switch i
        case 1
            BasicOnset(i).idx=idx_biOnset(corr_temp>=Threshold);
        case 2
            BasicOnset(i).idx=idx_bwdOnset(corr_temp>=Threshold);
        case 3
            BasicOnset(i).idx=idx_fwdOnset(corr_temp>=Threshold);
    end
    
    BasicOnset(i).mean=mean(ZS2(BasicOnset(i).idx,:),1);
    BasicOnset(i).std=std(ZS2(BasicOnset(i).idx,:),1);
end

figure;
ha=tight_subplot(1,3);
for i=1:3
    axes(ha(i));
    plot(BasicOnset(i).mean);hold on;plot(Speed_flow(1,:)/2);hold on;plot(-Speed_flow(2,:)/2);ylim([-1 5]);
end

ROI_temp=[];
for i=1:3    
    idx_speed=BasicOnset(i).idx;
    ROI_speed=ROI_rotated(intersect(idx_NotF9orF7,idx_speed),:);
    IsInBrainRegion=ismember(round(ROI_speed),Zbrain_AllMask,'rows');
    ROI_speed(:,4)=i;
    ROI_temp=vertcat(ROI_temp,ROI_speed(IsInBrainRegion,:));
end
csvwrite(strcat('BasicOnset_encoding_','BiFWDRev','_ROIs.csv'),ROI_temp);

Fighandle=figure;
set(Fighandle, 'Position', [10,10, 600, 1400]);

colors_onset = [0         0    1.0000         
    0    0.7000    0.2000
    1.0000         0         0];

Fighandle=figure;
set(Fighandle, 'Position', [10, 10, 800, 400]);
ha=tight_subplot(2,3,0.01);
for i=1:3
    axes(ha(i));
    plot(BasicOnset(i).mean(550:700)-min(BasicOnset(i).mean(550:700)),'Color',colors_onset(i,:),'LineWidth',4);ylim([-0.1 3.5]);
    hold on;a=area(Speed_flow(1,550:700),'FaceColor',[0.5 0.5 0.5]);a.FaceAlpha = 0.3;a.LineStyle='None';a.BaseLine.LineStyle='None';xlim([-5 155]);set(gca,'xtick',[]);set(gca,'ytick',[]);
    axes(ha(i+3));
    plot(BasicOnset(i).mean(1620:1770)-min(BasicOnset(i).mean(1620:1770)),'Color',colors_onset(i,:),'LineWidth',4);ylim([-0.1 3.5]);
    hold on;a=area(Speed_flow(2,1620:1770),'FaceColor',[0.5 0.5 0.5]);a.FaceAlpha = 0.3;a.LineStyle='None';a.BaseLine.LineStyle='None';xlim([-5 155]);set(gca,'xtick',[]);set(gca,'ytick',[]);
end
print(Fighandle,strcat('D:\Pictures\processed\Flow\Complex_final\Figure\OnsetEncoding_complexb.svg'),'-dsvg','-r0');


