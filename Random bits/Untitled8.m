Fighandle=figure;
set(Fighandle, 'Position', [10,10, 600, 1400]);
plot(Zbrain_brainMask2D(:,2),1400-Zbrain_brainMask2D(:,1),'k');hold on;
scatter(Node_all(:,2),Node_all(:,1),90,'k','filled');hold on;
scatter(Node_all(:,2),Node_all(:,1),55,P,'filled');colormap hot;hold off;axis([0 600 0 1400]);camroll(180);caxis([0 1]);
set(gca,'Visible','off')
print(Fighandle,strcat('D:\Pictures\processed\Flow\Basic\Figure\Graph Theory\','Participation.tif'),'-dtiff','-r0');
colorbar;caxis([0 1]);
print(Fighandle,strcat('D:\Pictures\processed\Flow\Basic\Figure\Graph Theory\','Participation_cbar.tif'),'-dtiff','-r0');


Fighandle=figure;
graph_temp=graph(Correlation_bin_prop);
set(Fighandle, 'Position', [10,10, 600, 1400]);
plot(Zbrain_brainMask2D(:,2),1400-Zbrain_brainMask2D(:,1),'k');hold on;
scatter(Node_all(:,2),Node_all(:,1),(degree(graph_temp)/5)+20,'k','filled');hold on;
LWidths = 1*graph_temp.Edges.Weight/max(graph_temp.Edges.Weight);
LWidths(~isfinite(LWidths))=0.01;
LColors = ones(length(LWidths),3);LColors=LColors-LWidths;
plot(graph_temp,'EdgeColor',LColors,'NodeColor',color_temp,'LineStyle','-','LineWidth',LWidths,'Marker','o','MarkerSize',degree(graph_temp)/5,'XData',Node_all(:,2),'YData',Node_all(:,1),'NodeLabel',char.empty(100,0));hold off;axis([0 600 0 1400]);camroll(180);
set(gca,'Visible','off')
print(Fighandle,strcat('D:\Pictures\processed\Flow\Basic\Figure\Graph Theory\','25PercGraph.tif'),'-dtiff','-r0');

[RdBu]=cbrewer('div','RdBu',101);
[Bu]=cbrewer('seq','Blues',101);
Fighandle=figure;
set(Fighandle, 'Position', [10,10, 1200, 1200]);
imagesc(Correlation_bin_prop,[0 1]);
colormap(Bu);
set(gca,'Visible','off')
print(Fighandle,strcat('D:\Pictures\processed\Flow\Basic\Figure\Graph Theory\','25PercCorrMatrix.tif'),'-dtiff','-r0');
colorbar;
print(Fighandle,strcat('D:\Pictures\processed\Flow\Basic\Figure\Graph Theory\','25PercCorrMatrix_cbar.tif'),'-dtiff','-r0');

Fighandle=figure;
set(Fighandle, 'Position', [10,10, 1200, 1200]);
imagesc(Correlation_bin_prop,[0 1]);
rectangle('FaceColor','r','Position',[10 -1 10 0.2]);
rectangle('FaceColor','r','Position',[50 -1 10 0.2]);
rectangle('FaceColor','r','Position',[90 -1 10 0.2]);
rectangle('FaceColor','b','Position',[30 -1 10 0.2]);
rectangle('FaceColor','b','Position',[70 -1 10 0.2]);
rectangle('FaceColor','b','Position',[110 -1 10 0.2]);