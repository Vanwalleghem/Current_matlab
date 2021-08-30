%% speed encoding per fish

idx_temp=Speed_encoders(2).idx;
ZS_temp=ZS2(idx_temp,:);
fish_temp=unique(idx_Fish(idx_temp));
ToPlot_slow=nan(length(idx_time_temp),length(unique(idx_Fish)),61,4);

idx_time_temp=[2 4 6 7];
for i=1:length(idx_time_temp)    
    for fish_nb=1:length(fish_temp)
        ToPlot_slow(i,fish_nb,:,1)=mean(ZS_temp(idx_Fish(idx_temp)==fish_temp(fish_nb),back(idx_time_temp(i)):back(idx_time_temp(i))+60));
    end
end
idx_time_temp=[1 5 8];
for i=1:length(idx_time_temp)    
    for fish_nb=1:length(fish_temp)
        ToPlot_slow(i,fish_nb,:,2)=mean(ZS_temp(idx_Fish(idx_temp)==fish_temp(fish_nb),back(idx_time_temp(i)):back(idx_time_temp(i))+60));
    end
end

idx_temp=Speed_encoders(2).idx;
ZS_temp=ZS2(idx_temp,:);
fish_temp=unique(idx_Fish(idx_temp));

idx_time_temp=[1 4 5 7];
for i=1:length(idx_time_temp)    
    for fish_nb=1:length(fish_temp)
        ToPlot_slow(i,fish_nb,:,3)=mean(ZS_temp(idx_Fish(idx_temp)==fish_temp(fish_nb),fwd(idx_time_temp(i)):fwd(idx_time_temp(i))+60));
    end
end
idx_time_temp=[2 3 9 10];
for i=1:length(idx_time_temp)    
    for fish_nb=1:length(fish_temp)
        ToPlot_slow(i,fish_nb,:,4)=mean(ZS_temp(idx_Fish(idx_temp)==fish_temp(fish_nb),fwd(idx_time_temp(i)):fwd(idx_time_temp(i))+60));
    end
end

Fighandle=figure;
set(Fighandle, 'Position', [10, 10, 800, 800]);
ha = tight_subplot(2,2,[.01 .01],[.01 .01],[.01 .01]);
ToPlot_temp2=squeeze(nanmean(ToPlot_slow,1));
ToPlot_temp2=squeeze(nanmean(ToPlot_temp2,1));
for i=1:4
    axes(ha(i));    
    plot(ToPlot_temp2(:,i));ylim([-0.5 3]);
end


idx_temp=Speed_encoders(1).idx;
ZS_temp=ZS2(idx_temp,:);
fish_temp=unique(idx_Fish(idx_temp));
ToPlot_fast=nan(length(idx_time_temp),length(unique(idx_Fish)),61,4);

idx_time_temp=[2 4 6 7];
for i=1:length(idx_time_temp)    
    for fish_nb=1:length(fish_temp)
        if sum(idx_Fish(idx_temp)==fish_temp(fish_nb))>5
            ToPlot_fast(i,fish_nb,:,1)=mean(ZS_temp(idx_Fish(idx_temp)==fish_temp(fish_nb),back(idx_time_temp(i)):back(idx_time_temp(i))+60));
        end
    end
end
idx_time_temp=[1 5 8];
for i=1:length(idx_time_temp)    
    for fish_nb=1:length(fish_temp)
        if sum(idx_Fish(idx_temp)==fish_temp(fish_nb))>5
        ToPlot_fast(i,fish_nb,:,2)=mean(ZS_temp(idx_Fish(idx_temp)==fish_temp(fish_nb),back(idx_time_temp(i)):back(idx_time_temp(i))+60));
        end
    end
end

idx_temp=Speed_encoders(1).idx;
ZS_temp=ZS2(idx_temp,:);
fish_temp=unique(idx_Fish(idx_temp));

idx_time_temp=[1 4 5 7];
for i=1:length(idx_time_temp)    
    for fish_nb=1:length(fish_temp)
        if sum(idx_Fish(idx_temp)==fish_temp(fish_nb))>5
        ToPlot_fast(i,fish_nb,:,3)=mean(ZS_temp(idx_Fish(idx_temp)==fish_temp(fish_nb),fwd(idx_time_temp(i)):fwd(idx_time_temp(i))+60));
        end
    end
end
idx_time_temp=[2 3 9 10];
for i=1:length(idx_time_temp)    
    for fish_nb=1:length(fish_temp)
        if sum(idx_Fish(idx_temp)==fish_temp(fish_nb))>5
        ToPlot_fast(i,fish_nb,:,4)=mean(ZS_temp(idx_Fish(idx_temp)==fish_temp(fish_nb),fwd(idx_time_temp(i)):fwd(idx_time_temp(i))+60));
        end
    end
end

Fighandle=figure;
set(Fighandle, 'Position', [10, 10, 800, 800]);
ha = tight_subplot(2,2,[.01 .01],[.01 .01],[.01 .01]);
ToPlot_temp2=squeeze(nanmean(ToPlot_fast,1));
ToPlot_temp2=squeeze(nanmean(ToPlot_temp2,1));
for i=1:4
    axes(ha(i));    
    plot(ToPlot_temp2(:,i));ylim([-0.5 3]);
end


%% final figure

Fighandle=figure;
set(Fighandle, 'Position', [10, 10, 800, 400]);
ha = tight_subplot(1,2,[.01 .01],[.01 .01],[.01 .01]);
axes(ha(2));
ToPlot_temp2=squeeze(nanmean(ToPlot_fast,1));
ToPlot_temp2=squeeze(nanmean(ToPlot_temp2,1));
plot(ToPlot_temp2(:,1)-min(ToPlot_temp2(:,1)),'r','LineWidth',3);ylim([-0.1 3]);hold on;
ToPlot_temp2=squeeze(nanmean(ToPlot_slow,1));
ToPlot_temp2=squeeze(nanmean(ToPlot_temp2,1));
plot(ToPlot_temp2(:,1)-min(ToPlot_temp2(:,1)),'b','LineWidth',3);ylim([-0.1 3]);hold off;
%set(gca,'visible','off');
set(gca,'xtick',[]);set(gca,'ytick',[]);
axes(ha(1));
ToPlot_temp2=squeeze(nanmean(ToPlot_fast,1));
ToPlot_temp2=squeeze(nanmean(ToPlot_temp2,1));
plot(ToPlot_temp2(:,2)-min(ToPlot_temp2(:,2)),'r','LineWidth',3);ylim([-0.1 3]);hold on;
ToPlot_temp2=squeeze(nanmean(ToPlot_slow,1));
ToPlot_temp2=squeeze(nanmean(ToPlot_temp2,1));
plot(ToPlot_temp2(:,2)-min(ToPlot_temp2(:,2)),'b','LineWidth',3);ylim([-0.1 3]);hold off;
%set(gca,'visible','off');
set(gca,'xtick',[]);set(gca,'ytick',[]);
print(Fighandle,strcat('D:\Pictures\processed\Flow\Complex_final\Figure\SpeedEncoding_complex.svg'),'-dsvg','-r0');
