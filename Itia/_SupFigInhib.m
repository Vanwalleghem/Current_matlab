colors=[0.14 1 0.14; 0.7 0.4 1;0 0.6 0.6];
Fighandle=figure;
set(Fighandle, 'Position', [100, 100, 1000, 600]);
ylimits=([5000 7500 ;3000 5500]);
%ha = tight_subplot(2,1,[.01 .01],[.01 .01],[.01 .01]);
for j=1:2
    %axes(ha(j));
    subplot(1,2,j);
    plot(RAW(:,j),'color',colors(j,:),'LineWidth',3);ylim(ylimits(j,:));xlim([0 760]);%set(gca,'XTickLabel',[]);set(gca,'YTickLabel',[]);
    %h=refline(0);h.LineStyle='- -';h.Color=[0.5 0.5 0.5];
    for i=1:5
        rectangle('EdgeColor','none','FaceColor',[0.3, 0.3, 0.3, 0.4-i*0.06],'Position',[40+120*(i-1) 0 4 8000]);
        rectangle('EdgeColor','none','FaceColor',[0.3, 0.3, 0.3, 0.4-i*0.06],'Position',[80+120*(i-1) 0 4 8000]);
        rectangle('EdgeColor','none','FaceColor',[0.3, 0.3, 0.3, 0.4-i*0.06],'Position',[120+120*(i-1) 0 4 8000]);
    end    
end