NullModel=zeros(size(Mean_allNodes_perFish));
for i=1:size(NullModel,1);
    NullModel(i,:,:)=AAFT(squeeze(nanmean(Mean_allNodes_perFish(i,:,:),2)),13)';
end

CorrelationMatrices_NullModel=nan(size(NullModel,1),size(NullModel,1),size(NullModel,2));
for fish_nb=1:size(NullModel,2)
    CorrelationMatrices_NullModel(:,:,fish_nb)=squareform(pdist(squeeze(NullModel(:,fish_nb,:)),'correlation'));
end
CorrelationMatrixNull=squeeze(nanmean(CorrelationMatrices_NullModel,3));

figure;imagesc(CorrelationMatrixNull,[0 2]);

Nodes_graph_Null=graph(CorrelationMatrixNull,'upper','omitselfloops');
Nodes_Null.density=degree(Nodes_graph_Null);

color_temp=zeros(size(Node_all));
for i=1:8
    idx_temp=find(Node_ID==i);
    for node_nb=1:length(idx_temp)
        color_temp(idx_temp(node_nb),:)=colors(i,:)/256;
    end
end

CorrelationMatrixNull=squeeze(nanmean(CorrelationMatrices_NullModel,3));
CorrelationMatrixNull=1-CorrelationMatrixNull;
CorrelationMatrixNull=weight_conversion(CorrelationMatrixNull,'autofix');
CorrelationMatrixNull=weight_conversion(CorrelationMatrixNull,'normalize');
figure;imagesc(CorrelationMatrixNull,[0 1]);
Correlation_bin_prop_Null=threshold_proportional(CorrelationMatrixNull,0.25);
figure;imagesc(Correlation_bin_prop_Null,[0 1]);
Correlation_bin_Null=threshold_absolute(CorrelationMatrixNull,0.75);
figure;imagesc(Correlation_bin_Null,[0 1]);

temp=strengths_und(Correlation_bin_prop_Null);

Strength_perClusterAndThr=zeros(8,10);
edges=[0.05:0.05:0.5];
for thr_nb=1:10
    graphTemp=threshold_proportional(CorrelationMatrixNull,edges(thr_nb));
    DensTemp=strengths_und(graphTemp);
    for i=1:8
        Strength_perClusterAndThr(i,thr_nb)=mean(DensTemp(Node_ID==i));
    end
end

Density_perClusterAndThr=zeros(8,10);
edges=[0.05:0.05:0.5];
for thr_nb=1:10
    graphTemp=threshold_proportional(CorrelationMatrixNull,edges(thr_nb));
    graphTemp=graph(graphTemp,'upper','omitselfloops');
    DensTemp=degree(graphTemp);
    for i=1:8
        Density_perClusterAndThr(i,thr_nb)=mean(DensTemp(Node_ID==i));
    end     
end

Fighandle=figure;
set(Fighandle, 'Position', [10,10, 1200, 1200]);
imagesc(CorrelationMatrixNull,[0 1]);colormap(Bu);
for i=1:8
    temp=find(Node_ID==i);    
    rectangle('FaceColor',colors(i,:)/256,'Position',[temp(1)-1 -5 temp(1)-1+length(temp) 5]);
    rectangle('FaceColor',colors(i,:)/256,'Position',[-5 temp(1)-1 5 temp(1)-1+length(temp)]);    
end
axis([-5 length(Node_ID) -5 length(Node_ID)]);
set(gca,'Visible','off')
print(Fighandle,strcat('D:\Pictures\processed\Flow\Basic\Figure\Graph Theory\','CorrMatrix_Null.tif'),'-dtiff','-r0');

% Fighandle=figure;
% graph_temp=graph(Correlation_bin_prop_Null);
% set(Fighandle, 'Position', [10,10, 600, 1400]);
% plot(Zbrain_brainMask2D(:,2),1400-Zbrain_brainMask2D(:,1),'k');hold on;
% scatter(Node_all(:,2),Node_all(:,1),(degree(graph_temp)/6)+10,'k','filled');hold on;
% %plot(Nodes_graph,'NodeColor','k','LineStyle','None','Marker','o','MarkerSize',log(Nodes.density)+3,'XData',Node_all(:,2),'YData',Node_all(:,1),'NodeLabel',char.empty(100,0));hold on;
% LWidths = 1*graph_temp.Edges.Weight/max(graph_temp.Edges.Weight);
% LWidths(~isfinite(LWidths))=0.01;
% LColors = ones(length(LWidths),3);LColors=LColors-LWidths;
% plot(graph_temp,'EdgeColor',LColors,'NodeColor',color_temp,'LineStyle','-','LineWidth',LWidths,'Marker','o','MarkerSize',degree(graph_temp)/6,'XData',Node_all(:,2),'YData',Node_all(:,1),'NodeLabel',char.empty(100,0));hold off;
% 
% [N,E] = rentian_scaling_3d(Correlation_bin_prop_Null>0,Node_all,10000,1e-6);
% figure; loglog(N,E,'.');
% [b,stats] = robustfit(log10(N),log10(E));[b(2) stats.se(2)]
% 
% P_Null=participation_coef(Correlation_bin_prop_Null,Node_ID);
% Fighandle=figure;
% set(Fighandle, 'Position', [1000,1000, 600, 1400]);
% plot(Zbrain_brainMask2D(:,2),1400-Zbrain_brainMask2D(:,1),'k');hold on;
% scatter(Node_all(:,2),Node_all(:,1),55,'k','filled');hold on;
% scatter(Node_all(:,2),Node_all(:,1),45,P_Null,'filled');colormap hot;hold off;

CorrelationMatrixNull=squeeze(nanmean(CorrelationMatrices_NullModel,3));
CorrelationMatrixNull=1-CorrelationMatrixNull;
CorrelationMatrixNull=abs(CorrelationMatrixNull);
CorrelationMatrixNull=weight_conversion(CorrelationMatrixNull,'autofix');
%CorrelationMatrix=weight_conversion(CorrelationMatrix,'normalize');
figure;imagesc(CorrelationMatrixNull,[0 1]);

Correlation_bin_prop_Null=threshold_absolute(CorrelationMatrixNull,0.6);
figure;imagesc(Correlation_bin_prop_Null,[0 1]);

Fighandle=figure;
graph_temp=graph(Correlation_bin_prop_Null);
set(Fighandle, 'Position', [10,10, 600, 1400]);
plot(Zbrain_brainMask2D(:,2),1400-Zbrain_brainMask2D(:,1),'k');hold on;
scatter(Node_all(:,2),Node_all(:,1),150,'k','filled');hold on;
LWidths = 1*graph_temp.Edges.Weight/max(graph_temp.Edges.Weight);
LWidths(~isfinite(LWidths))=0.01;
LWidths(LWidths==0)=0.01;
LColors = ones(length(LWidths),3);LColors=LColors-LWidths;
plot(graph_temp,'EdgeColor',LColors,'NodeColor',color_temp,'LineStyle','-','LineWidth',LWidths,'Marker','o','MarkerSize',10,'XData',Node_all(:,2),'YData',Node_all(:,1),'NodeLabel',char.empty(100,0));hold off;axis([0 600 0 1400]);camroll(180);
ax = gca;
outerpos = ax.OuterPosition;
ti = ax.TightInset; 
left = outerpos(1) + ti(1);
bottom = outerpos(2) + ti(2);
ax_width = outerpos(3) - ti(1) - ti(3);
ax_height = outerpos(4) - ti(2) - ti(4);
ax.Position = [left bottom ax_width ax_height];
set(gca,'Visible','off')
print(Fighandle,strcat('D:\Pictures\processed\Flow\Basic\Figure\Graph Theory\','Edges_06_Null.tif'),'-dtiff','-r0');


P=participation_coef(Correlation_bin_prop,Node_ID);
Fighandle=figure;Multiplier=300;Baseline=0.01;
set(Fighandle, 'Position', [10,10, 600, 1400]);
plot(Zbrain_brainMask2D(:,2),1400-Zbrain_brainMask2D(:,1),'k');hold on;
scatter(Node_all(:,2),Node_all(:,1),20+(P+Baseline)*Multiplier,[0.2 0.2 0.2],'filled');
scatter(Node_all(:,2),Node_all(:,1),(P+Baseline)*Multiplier,color_temp,'filled');hold on;axis([0 600 0 1400]);camroll(180);caxis([0 1]);
ax = gca;
outerpos = ax.OuterPosition;
ti = ax.TightInset; 
left = outerpos(1) + ti(1);
bottom = outerpos(2) + ti(2);
ax_width = outerpos(3) - ti(1) - ti(3);
ax_height = outerpos(4) - ti(2) - ti(4);
ax.Position = [left bottom ax_width ax_height];
set(gca,'Visible','off')
print(Fighandle,strcat('D:\Pictures\processed\Flow\Basic\Figure\Graph Theory\','Participation_size06_Null.tif'),'-dtiff','-r0');
