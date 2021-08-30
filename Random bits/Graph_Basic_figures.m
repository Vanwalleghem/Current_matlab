color_temp=zeros(length(Node_ID),3);
for i=1:length(Node_ID)
    color_temp(i,:)=colors(Node_ID(i),:);
end
color_temp=color_temp/256;

%Move Rev_ON pLLG 10px
Node_all_REVON_backup=Node_all(107,:);
Node_all(107,2)=Node_all(107,2)+5;
Node_all(107,1)=Node_all(107,1)-5;

Fighandle=figure;
set(Fighandle, 'Position', [10,10, 600, 1400]);
plot(Zbrain_brainMask2D(:,2),1400-Zbrain_brainMask2D(:,1),'k');hold on;
scatter(Node_all(:,2),Node_all(:,1),150,'k','filled');hold on;
scatter(Node_all(:,2),Node_all(:,1),130,color_temp,'filled');colormap hot;hold off;axis([0 600 0 1400]);camroll(180);caxis([0 1]);
set(gca,'Visible','off');
ax = gca;
outerpos = ax.OuterPosition;
ti = ax.TightInset; 
left = outerpos(1) + ti(1);
bottom = outerpos(2) + ti(2);
ax_width = outerpos(3) - ti(1) - ti(3);
ax_height = outerpos(4) - ti(2) - ti(4);
ax.Position = [left bottom ax_width ax_height];
print(Fighandle,strcat('D:\Pictures\processed\Flow\Basic\Figure\Graph Theory\','Nodes.tif'),'-dtiff','-r0');

CorrelationMatrix=squeeze(nanmean(CorrelationMatrices_perFish,3));
CorrelationMatrix=1-CorrelationMatrix;
CorrelationMatrix=abs(CorrelationMatrix);
CorrelationMatrix=weight_conversion(CorrelationMatrix,'autofix');
%CorrelationMatrix=weight_conversion(CorrelationMatrix,'normalize');
figure;imagesc(CorrelationMatrix,[0 1]);
% Correlation_bin_prop=threshold_proportional(CorrelationMatrix,0.25);
Correlation_bin_prop=threshold_absolute(CorrelationMatrix,0.6);

Fighandle=figure;
set(Fighandle, 'Position', [10,10, 1200, 1200]);
imagesc(CorrelationMatrix,[0 1]);colormap(Bu);
for i=1:8
    temp=find(Node_ID==i);    
    rectangle('FaceColor',colors(i,:)/256,'Position',[temp(1)-1 -5 temp(1)-1+length(temp) 5]);
    rectangle('FaceColor',colors(i,:)/256,'Position',[-5 temp(1)-1 5 temp(1)-1+length(temp)]);    
end
axis([-5 length(Node_ID) -5 length(Node_ID)]);
set(gca,'Visible','off')
print(Fighandle,strcat('D:\Pictures\processed\Flow\Basic\Figure\Graph Theory\','CorrMatrix.tif'),'-dtiff','-r0');


Fighandle=figure;
graph_temp=graph(Correlation_bin_prop);
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
print(Fighandle,strcat('D:\Pictures\processed\Flow\Basic\Figure\Graph Theory\','Edges_06.tif'),'-dtiff','-r0');


P=participation_coef(Correlation_bin_prop,Node_ID);
Fighandle=figure;Multiplier=300;Baseline=0.01;
set(Fighandle, 'Position', [10,10, 600, 1400]);
plot(Zbrain_brainMask2D(:,2),1400-Zbrain_brainMask2D(:,1),'k');hold on;
scatter(Node_all(:,2),Node_all(:,1),20+(P+Baseline)*Multiplier,[0.2 0.2 0.2],'filled');
scatter(Node_all(:,2),Node_all(:,1),(P+Baseline)*Multiplier,color_temp,'filled');hold on;axis([0 600 0 1400]);camroll(180);caxis([0 1]);
counter=0;
for degree_sz = 0:0.1:0.6
    scatter(550,1100-7*counter,(degree_sz+Baseline)*Multiplier,'k','LineWidth',2);hold on;
    counter=counter+5;
end
ax = gca;
outerpos = ax.OuterPosition;
ti = ax.TightInset; 
left = outerpos(1) + ti(1);
bottom = outerpos(2) + ti(2);
ax_width = outerpos(3) - ti(1) - ti(3);
ax_height = outerpos(4) - ti(2) - ti(4);
ax.Position = [left bottom ax_width ax_height];
set(gca,'Visible','off')
print(Fighandle,strcat('D:\Pictures\processed\Flow\Basic\Figure\Graph Theory\','Participation_size06b.tif'),'-dtiff','-r0');

Particip_perFish=nan(8,13);
for i=1:13
    CorrelationMatrix=squeeze(nanmean(CorrelationMatrices_perFish,3));
    CorrelationMatrix=1-CorrelationMatrix;
    CorrelationMatrix=abs(CorrelationMatrix);
    CorrelationMatrix=weight_conversion(CorrelationMatrix,'autofix');
    %CorrelationMatrix=weight_conversion(CorrelationMatrix,'normalize');
    Correlation_bin_prop=threshold_absolute(CorrelationMatrix,0.6);    
    P=participation_coef(Correlation_bin_prop,Node_ID);
end


% Fighandle=figure;
% set(Fighandle, 'Position', [10,10, 600, 1400]);
% plot(Zbrain_brainMask2D(:,2),1400-Zbrain_brainMask2D(:,1),'k');hold on;
% scatter(Node_all(:,2),Node_all(:,1),90,color_temp,'filled');hold on;
% scatter(Node_all(:,2),Node_all(:,1),55,P,'filled');colormap gray;hold off;axis([0 600 0 1400]);camroll(180);caxis([0 1]);
% set(gca,'Visible','off')
% print(Fighandle,strcat('D:\Pictures\processed\Flow\Basic\Figure\Graph Theory\','Participation_greydot.tif'),'-dtiff','-r0');

Fighandle=figure;
graph_temp=graph(Correlation_bin_prop);
set(Fighandle, 'Position', [10,10, 600, 1400]);
plot(Zbrain_brainMask2D(:,2),1400-Zbrain_brainMask2D(:,1),'k');hold on;
scatter(Node_all(:,2),Node_all(:,1),20+(degree(graph_temp))*6,'k','filled');
scatter(Node_all(:,2),Node_all(:,1),5+(degree(graph_temp))*6,color_temp,'filled');hold on;axis([0 600 0 1400]);camroll(180);caxis([0 1]);
for degree_sz = 0:5:40
    scatter(550,1100-7*degree_sz,5+degree_sz*6,'k','LineWidth',2);hold on;
end
ax = gca;
outerpos = ax.OuterPosition;
ti = ax.TightInset; 
left = outerpos(1) + ti(1);
bottom = outerpos(2) + ti(2);
ax_width = outerpos(3) - ti(1) - ti(3);
ax_height = outerpos(4) - ti(2) - ti(4);
ax.Position = [left bottom ax_width ax_height];
set(gca,'Visible','off')
print(Fighandle,strcat('D:\Pictures\processed\Flow\Basic\Figure\Graph Theory\','Degree_size06b.tif'),'-dtiff','-r0');

% Fighandle=figure;
% set(Fighandle, 'Position', [10,10, 600, 1400]);
% plot(Zbrain_brainMask2D(:,2),1400-Zbrain_brainMask2D(:,1),'k');hold on;
% scatter(Node_all(:,2),Node_all(:,1),90,'k','filled');hold on;
% scatter(Node_all(:,2),Node_all(:,1),55,P,'filled');colormap hot;hold off;axis([0 600 0 1400]);camroll(180);caxis([0 1]);
% set(gca,'Visible','off')
% print(Fighandle,strcat('D:\Pictures\processed\Flow\Basic\Figure\Graph Theory\','Participation.tif'),'-dtiff','-r0');
% colorbar;caxis([0 1]);
% print(Fighandle,strcat('D:\Pictures\processed\Flow\Basic\Figure\Graph Theory\','Participation_cbar.tif'),'-dtiff','-r0');

Particip_perFish=nan(13,8);
for fish_nb=1:13
    CorrTemp=squeeze(CorrelationMatrices_minusFish(:,:,fish_nb));
    CorrTemp=1-CorrTemp;
    CorrTemp=abs(CorrTemp);
    CorrTemp=weight_conversion(CorrTemp,'autofix');    
    %CorrTemp=threshold_absolute(CorrTemp,0.5);    
    CorrTemp=threshold_proportional(CorrTemp,0.25);    
    P=participation_coef(CorrTemp,Node_ID);
    for cluster_id=1:8
        Particip_perFish(fish_nb,cluster_id)=mean(P(Node_ID==cluster_id));
    end
end

Degree_perFish=nan(13,8);
for fish_nb=1:13
    CorrTemp=squeeze(CorrelationMatrices_minusFish(:,:,fish_nb));
    CorrTemp=1-CorrTemp;
    CorrTemp=abs(CorrTemp);
    CorrTemp=weight_conversion(CorrTemp,'autofix');    
    CorrTemp=threshold_absolute(CorrTemp,0.5);    
    graph_temp=graph(CorrTemp);
    P=degree(graph_temp);
    for cluster_id=1:8
        Degree_perFish(fish_nb,cluster_id)=mean(P(Node_ID==cluster_id));
    end
end

% Fighandle=figure;
% graph_temp=graph(Correlation_bin_prop);
% set(Fighandle, 'Position', [10,10, 600, 1400]);
% plot(Zbrain_brainMask2D(:,2),1400-Zbrain_brainMask2D(:,1),'k');hold on;
% scatter(Node_all(:,2),Node_all(:,1),(degree(graph_temp)/3)+20,'k','filled');hold on;
% LWidths = 1*graph_temp.Edges.Weight/max(graph_temp.Edges.Weight);
% LWidths(~isfinite(LWidths))=0.01;
% LWidths(LWidths==0)=0.01;
% LColors = ones(length(LWidths),3);LColors=LColors-LWidths;
% plot(graph_temp,'EdgeColor',LColors,'NodeColor',color_temp,'LineStyle','-','LineWidth',LWidths,'Marker','o','MarkerSize',5+degree(graph_temp)/3,'XData',Node_all(:,2),'YData',Node_all(:,1),'NodeLabel',char.empty(100,0));hold off;axis([0 600 0 1400]);camroll(180);
% set(gca,'Visible','off')
% print(Fighandle,strcat('D:\Pictures\processed\Flow\Basic\Figure\Graph Theory\','25PercGraph.tif'),'-dtiff','-r0');

[RdBu]=cbrewer('div','RdBu',101);
[Bu]=cbrewer('seq','Blues',101);
% Fighandle=figure;
% set(Fighandle, 'Position', [10,10, 1200, 1200]);
% imagesc(Correlation_bin_prop,[0 1]);
% colormap(Bu);
% set(gca,'Visible','off')
% print(Fighandle,strcat('D:\Pictures\processed\Flow\Basic\Figure\Graph Theory\','25PercCorrMatrix.tif'),'-dtiff','-r0');
% colorbar;
% print(Fighandle,strcat('D:\Pictures\processed\Flow\Basic\Figure\Graph Theory\','25PercCorrMatrix_cbar.tif'),'-dtiff','-r0');

Fighandle=figure;
set(Fighandle, 'Position', [10,10, 1200, 1200]);
imagesc(Correlation_bin_prop,[0 1]);colormap(Bu);
for i=1:8
    temp=find(Node_ID==i);    
    rectangle('FaceColor',colors(i,:)/256,'Position',[temp(1)-1 -5 temp(1)-1+length(temp) 5]);
    rectangle('FaceColor',colors(i,:)/256,'Position',[-5 temp(1)-1 5 temp(1)-1+length(temp)]);    
end
axis([-5 length(Node_ID) -5 length(Node_ID)]);
set(gca,'Visible','off')
print(Fighandle,strcat('D:\Pictures\processed\Flow\Basic\Figure\Graph Theory\','CorrMatrix06.tif'),'-dtiff','-r0');

% Correlation_bin_prop=threshold_proportional(CorrelationMatrix,0.1);
% Correlation_bin_prop=threshold_absolute(CorrelationMatrix,0.7);
Correlation_matrix_binary=Correlation_bin_prop>0;

EdgesColorMap =Correlation_bin_prop(Correlation_matrix_binary>0);
EdgesColorMap=repmat([0.5 0.5 0.5],length(EdgesColorMap),1).*(1-EdgesColorMap);
NodeColorsFct=zeros(size(Node_all));
for i=1:8
    NodeColorsFct(Node_ID==i,:)=repmat(colors(i,:)/256,sum(Node_ID==i),1);
end
label_temp=cell(1,length(NodeColorsFct));
for i=1:size(Node_all,1)
    %label_temp{i}=RegionList{brain_nb_temp(i)};
    label_temp{i}=[' '];
end

Fighandle=figure;
set(Fighandle, 'Position', [50,50, 1600, 1600]);
h=circularGraph(Correlation_bin_prop,'ColorMapNode',NodeColorsFct,'ColorMapEdges',EdgesColorMap,'Label',label_temp);
%h=circularGraph(Correlation_matrix_binary,'ColorMapNode',color_temp,'Label',label_temp);
positions=[];
for i=1:length(h.Node)
    positions(i,:)=h.Node(i).Position;
end
hold on;
scatter(positions(:,1),positions(:,2),150,'k','filled');
hold on;
scatter(positions(:,1),positions(:,2),130,NodeColorsFct,'filled');
hold off;
print(Fighandle,strcat('D:\Pictures\processed\Flow\Basic\Figure\Graph Theory\','CircularPerFct06.tif'),'-dtiff','-r0');
close all;

[B,sort_brain]=sort(Node_selectRegion);

CorrelationMatrices_perFish_sortedBrain=nan(size(Mean_allNodes_perFish,1),size(Mean_allNodes_perFish,1),size(Mean_allNodes_perFish,2));
for fish_nb=1:size(Mean_allNodes_perFish,2)
    CorrelationMatrices_perFish_sortedBrain(:,:,fish_nb)=squareform(pdist(squeeze(Mean_allNodes_perFish(sort_brain,fish_nb,:)),'correlation'));
end

CorrelationMatrix_sortedBrain=squeeze(nanmean(CorrelationMatrices_perFish_sortedBrain,3));
CorrelationMatrix_sortedBrain=1-CorrelationMatrix_sortedBrain;
CorrelationMatrix_sortedBrain=weight_conversion(CorrelationMatrix_sortedBrain,'autofix');
% CorrelationMatrix_sortedBrain=weight_conversion(CorrelationMatrix_sortedBrain,'normalize');
figure;imagesc(CorrelationMatrix_sortedBrain,[0 1]);
Correlation_bin_prop_sortedBrain=threshold_proportional(CorrelationMatrix_sortedBrain,0.1);
Correlation_bin_prop_sortedBrain=threshold_absolute(CorrelationMatrix_sortedBrain,0.6);

Correlation_matrix_binary=Correlation_bin_prop_sortedBrain>0;

EdgesColorMap =Correlation_bin_prop_sortedBrain(Correlation_matrix_binary>0);
EdgesColorMap=repmat([0.5 0.5 0.5],length(EdgesColorMap),1).*(1-EdgesColorMap);

label_temp=cell(1,length(NodeColorsFct));
for i=1:size(B,2)
    %label_temp{i}=RegionList{brain_nb_temp(i)};
    label_temp{i}=RegionList_select2{B(i)};
end

Fighandle=figure;
set(Fighandle, 'Position', [50,50, 1600, 1600]);
h=circularGraph(Correlation_matrix_binary,'ColorMapNode',NodeColorsFct(sort_brain,:),'ColorMapEdges',EdgesColorMap,'Label',label_temp);
%h=circularGraph(Correlation_matrix_binary,'ColorMapNode',color_temp,'Label',label_temp);
positions=[];
for i=1:length(h.Node)
    positions(i,:)=h.Node(i).Position;
end
hold on;
scatter(positions(:,1),positions(:,2),120,'k','filled');
hold on;
scatter(positions(:,1),positions(:,2),100,NodeColorsFct(sort_brain,:),'filled');
hold off;
print(Fighandle,strcat('D:\Pictures\processed\Flow\Basic\Figure\Graph Theory\','CircularPerBrain06.tif'),'-dtiff','-r0');
close all;

Particip_perFish=nan(13,8);
for fish_nb=1:13
    CorrTemp=squeeze(CorrelationMatrices_minusFish(:,:,fish_nb));
    CorrTemp=1-CorrTemp;
    CorrTemp=abs(CorrTemp);
    CorrTemp=weight_conversion(CorrTemp,'autofix');    
    CorrTemp=threshold_absolute(CorrTemp,0.6);    
    %CorrTemp=threshold_proportional(CorrTemp,0.25);    
    P=participation_coef(CorrTemp,Node_ID);
    for cluster_id=1:8
        Particip_perFish(fish_nb,cluster_id)=mean(P(Node_ID==cluster_id));
    end
end

Degree_perFish=nan(13,8);
for fish_nb=1:13
    CorrTemp=squeeze(CorrelationMatrices_minusFish(:,:,fish_nb));
    CorrTemp=1-CorrTemp;
    CorrTemp=abs(CorrTemp);
    CorrTemp=weight_conversion(CorrTemp,'autofix');    
    CorrTemp=threshold_absolute(CorrTemp,0.6);    
    graph_temp=graph(CorrTemp);
    P=degree(graph_temp);
    for cluster_id=1:8
        Degree_perFish(fish_nb,cluster_id)=mean(P(Node_ID==cluster_id));
    end
end


Degree_perFish_perCluster=[];
for thr_nb = [0.9:-0.05:0.5]
    Degree_perFish=nan(8,13);
    for fish_nb=1:13
        
        CorrTemp=squeeze(CorrelationMatrices_minusFish(:,:,fish_nb));
        CorrTemp=1-CorrTemp;
        CorrTemp=abs(CorrTemp);
        CorrTemp=weight_conversion(CorrTemp,'autofix');
        CorrTemp=threshold_absolute(CorrTemp,thr_nb);
        graph_temp=graph(CorrTemp);
        P=degree(graph_temp);
        for cluster_id=1:8
            Degree_perFish(cluster_id,fish_nb)=mean(P(Node_ID==cluster_id));
        end
    end
    Degree_perFish_perCluster=horzcat(Degree_perFish_perCluster,Degree_perFish);
end

Part_perFish_perCluster=[];
for thr_nb = [0.9:-0.05:0.5]
    Degree_perFish=nan(8,13);
    for fish_nb=1:13
        
        CorrTemp=squeeze(CorrelationMatrices_minusFish(:,:,fish_nb));
        CorrTemp=1-CorrTemp;
        CorrTemp=abs(CorrTemp);
        CorrTemp=weight_conversion(CorrTemp,'autofix');
        CorrTemp=threshold_absolute(CorrTemp,thr_nb);        
        P=participation_coef(CorrTemp,Node_ID);
        for cluster_id=1:8
            Degree_perFish(cluster_id,fish_nb)=mean(P(Node_ID==cluster_id));
        end
    end
    Part_perFish_perCluster=horzcat(Part_perFish_perCluster,Degree_perFish);
end
