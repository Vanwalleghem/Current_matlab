Fighandle=figure;
set(Fighandle, 'Position', [100, 100, 1400, 900]);x = linspace(0.2,size(Cmap,2)/5,size(Cmap,2));
counter=1;counter2=1;xplot=floor(sqrt(length(GoodBetas_select)));yplot=ceil(length(GoodBetas_select)/xplot);
for i=1:length(GoodBetas_select)  
    if counter==3
        counter=counter+1;
    end    
    subplot(3,3,counter);plot(x,mean(ZS(idxKmeans_ZS_goodmembers==(GoodBetas_select(i)),:),1),'color',colors(counter2,:)/256);axis([0 131 -1 4]);rectangle('FaceColor','r','Position',[11 -1 10 0.25]);rectangle('FaceColor','r','Position',[51 -1 10 0.25]);rectangle('FaceColor','r','Position',[91 -1 10 0.25]);rectangle('FaceColor','b','Position',[31 -1 10 0.25]);rectangle('FaceColor','b','Position',[71 -1 10 0.25]);rectangle('FaceColor','b','Position',[111 -1 10 0.25]);
    counter=counter+1;
    counter2=counter2+1;
end

options = statset('UseParallel',1); %parallelize the replicates
Nodes_coords={};
for i=1:length(GoodBetas_select)  
    idx_temp=find(idxKmeans_ZS_goodmembers==GoodBetas_select(i));
    ROI_temp=ROI_fish(idx_temp,:);
    Fish_temp=idx_Fish(idx_temp);
    moduleN=10+ceil(length(ROI_temp)/20);
    test_fish=1;
    while moduleN>1 & test_fish==1
        [idxKmeans_ROI Cmap_ROI]=kmeans(ROI_temp,moduleN,'Options',options,'Replicates',3,'MaxIter',1000);
        for roi_nb=1:max(idxKmeans_ROI)
            idx_roi=find(idxKmeans_ROI==roi_nb);
            Fish_sum=histcounts(categorical(Fish_temp(idx_roi),Fish_list));
            if sum(Fish_sum>10)<3
                moduleN=moduleN-1;
                break
            else
                test_fish=0;
            end
        end
    end
    Nodes_coords{i,1}=Cmap_ROI;
    Nodes_coords{i,2}=idxKmeans_ROI;
end
clearvars i moduleN test_fish idx_roi idx_temp ROI_temp idxKmeans_ROI Cmap_ROI

Node_all=[];
Mean_allNodes=[];
Mean_allNodes_perFish=[];
Node_ID=[];
counter=1;
for i=1:length(Nodes_coords)
    Node_all=cat(1,Node_all,Nodes_coords{i,1});
    idx_temp=find(idxKmeans_ZS_goodmembers==GoodBetas_select(i));
    Fish_temp=idx_Fish(idx_temp);
    idxKmeans_ROI=Nodes_coords{i,2};
    for roi_nb=1:max(idxKmeans_ROI)
        idx_roi=find(idxKmeans_ROI==roi_nb);
        Mean_allNodes(counter,:)=mean(ZS(idx_temp(idx_roi),:),1);
        Fish_temp2=Fish_temp(idx_roi);
        for fish_nb=1:length(Fish_list)       
            idx_temp2=find(Fish_temp2==Fish_list(fish_nb));
            if idx_temp2
                Mean_allNodes_perFish(counter,fish_nb,:)=mean(ZS(idx_temp(idx_roi(idx_temp2)),:),1);
            else
                Mean_allNodes_perFish(counter,fish_nb,:)=nan(1,size(ZS,2));
            end
        end
        counter=counter+1;
    end
    Node_ID=cat(1,Node_ID,ones(size(Nodes_coords{i,1},1),1)*i);
end
clearvars i moduleN test_fish idx_roi idx_temp ROI_temp idxKmeans_ROI Cmap_ROI Fish_temp Fish_temp2 idx_temp2

Node_all(find(round(Node_all(:,1))==109),2)=270;

CorrelationMatrices_perFish=nan(size(Mean_allNodes_perFish,1),size(Mean_allNodes_perFish,1),size(Mean_allNodes_perFish,2));
for fish_nb=1:size(Mean_allNodes_perFish,2)
    CorrelationMatrices_perFish(:,:,fish_nb)=squareform(pdist(squeeze(Mean_allNodes_perFish(:,fish_nb,:)),'correlation'));
end
CorrelationMatrix=squeeze(nanmean(CorrelationMatrices_perFish,3));
Nodes_graph=graph(CorrelationMatrix,'upper','omitselfloops');
Nodes.density=degree(Nodes_graph);

color_temp=zeros(size(Node_all));
for i=1:8
    idx_temp=find(Node_ID==i);
    for node_nb=1:length(idx_temp)
        color_temp(idx_temp(node_nb),:)=colors(i,:)/256;
    end
end

CorrelationMatrix=squeeze(nanmean(CorrelationMatrices_perFish,3));
CorrelationMatrix=abs(1-CorrelationMatrix);
CorrelationMatrix=weight_conversion(CorrelationMatrix,'autofix');
CorrelationMatrix=weight_conversion(CorrelationMatrix,'normalize');
figure;imagesc(CorrelationMatrix,[0 1]);
Correlation_bin_prop=threshold_proportional(CorrelationMatrix,0.25);
figure;imagesc(Correlation_bin_prop,[0 1]);
Correlation_bin=threshold_absolute(CorrelationMatrix,0.75);
figure;imagesc(Correlation_bin,[0 1]);

temp=strengths_und(Correlation_bin_prop);

Strength_perClusterAndThr=zeros(8,10);
edges=[0.05:0.05:0.5];
for thr_nb=1:10
    graphTemp=threshold_proportional(CorrelationMatrix,edges(thr_nb));
    DensTemp=strengths_und(graphTemp);
    for i=1:8
        Strength_perClusterAndThr(i,thr_nb)=mean(DensTemp(Node_ID==i));
    end
end

Density_perClusterAndThr=zeros(8,10);
edges=[0.05:0.05:0.5];
for thr_nb=1:10
    graphTemp=threshold_proportional(CorrelationMatrix,edges(thr_nb));
    graphTemp=graph(graphTemp,'upper','omitselfloops');
    DensTemp=degree(graphTemp);
    for i=1:8
        Density_perClusterAndThr(i,thr_nb)=mean(DensTemp(Node_ID==i));
    end     
end

Correlation_bin_prop=threshold_proportional(CorrelationMatrix,0.25);

Fighandle=figure;
graph_temp=graph(Correlation_bin_prop);
set(Fighandle, 'Position', [10,10, 600, 1400]);
plot(Zbrain_brainMask2D(:,2),1400-Zbrain_brainMask2D(:,1),'k');hold on;
scatter(Node_all(:,2),Node_all(:,1),(degree(graph_temp)/6)+10,'k','filled');hold on;
%plot(Nodes_graph,'NodeColor','k','LineStyle','None','Marker','o','MarkerSize',log(Nodes.density)+3,'XData',Node_all(:,2),'YData',Node_all(:,1),'NodeLabel',char.empty(100,0));hold on;
LWidths = 1*graph_temp.Edges.Weight/max(graph_temp.Edges.Weight);
LWidths(~isfinite(LWidths))=0.01;
LColors = ones(length(LWidths),3);LColors=LColors-LWidths;
plot(graph_temp,'EdgeColor',LColors,'NodeColor',color_temp,'LineStyle','-','LineWidth',LWidths,'Marker','o','MarkerSize',degree(graph_temp)/6,'XData',Node_all(:,2),'YData',Node_all(:,1),'NodeLabel',char.empty(100,0));hold off;



temp=gtom(Correlation_bin_prop,3);

[N,E] = rentian_scaling_3d(CorrelationMatrix>0.75,Node_all,5000,1e-6);
figure; loglog(E,N,'.');
[b,stats] = robustfit(log10(N),log10(E));[b(2) stats.se(2)]


Fighandle=figure;
set(Fighandle, 'Position', [1000,1000, 600, 1400]);
plot(Zbrain_brainMask2D(:,2),1400-Zbrain_brainMask2D(:,1),'k');hold on;
scatter(Node_all(:,2),Node_all(:,1),log(Nodes.density)+3,'k','filled');hold on;
%plot(Nodes_graph,'NodeColor','k','LineStyle','None','Marker','o','MarkerSize',log(Nodes.density)+3,'XData',Node_all(:,2),'YData',Node_all(:,1),'NodeLabel',char.empty(100,0));hold on;
LWidths = 1*Nodes_graph.Edges.Weight/max(Nodes_graph.Edges.Weight);
LWidths(~isfinite(LWidths))=0.01;
LColors = ones(length(LWidths),3);LColors=LColors-LWidths;
plot(Nodes_graph,'EdgeColor',LColors,'NodeColor',color_temp,'LineStyle','-','LineWidth',LWidths,'Marker','o','MarkerSize',log(Nodes.density),'XData',Node_all(:,2),'YData',Node_all(:,1),'NodeLabel',char.empty(100,0));hold off;

Nodes.graphThr=graph(CorrelationMatrix<0.35,'upper','omitselfloops');
Nodes.densityThr=degree(Nodes.graphThr);

Fighandle=figure;
set(Fighandle, 'Position', [1000,1000, 600, 1400]);
%plot(Zbrain_brainMask2D(:,2),Zbrain_brainMask2D(:,1),'k');hold on;
scatter(Node_all(:,2),Node_all(:,1),(Nodes.densityThr/5)+3,'k','filled');hold on;
%plot(Nodes_graph,'NodeColor','k','LineStyle','None','Marker','o','MarkerSize',log(Nodes.density)+3,'XData',Node_all(:,2),'YData',Node_all(:,1),'NodeLabel',char.empty(100,0));hold on;
plot(Nodes.graphThr,'EdgeColor',[0.2 0.2 0.2],'NodeColor',color_temp,'LineStyle','-','LineWidth',0.5,'Marker','o','MarkerSize',(Nodes.densityThr/5)+3,'XData',Node_all(:,2),'YData',Node_all(:,1),'NodeLabel',char.empty(100,0));hold off;

Density_perClusterAndThr=zeros(8,10);
edges=[0.5:-0.05:0.05];
for thr_nb=1:10
    graphTemp=graph(CorrelationMatrix<edges(thr_nb),'upper','omitselfloops');
    DensTemp=degree(graphTemp);
    for i=1:8
        Density_perClusterAndThr(i,thr_nb)=mean(DensTemp(Node_ID==i));
    end
end



Fighandle=figure;
set(Fighandle, 'Position', [10,10, 600, 1400]);
plot(Zbrain_brainMask2D(:,2),1400-Zbrain_brainMask2D(:,1),'k');hold on;
scatter(Node_all(:,2),Node_all(:,1),90,'k','filled');hold on;
scatter(Node_all(:,2),Node_all(:,1),55,P,'filled');colormap hot;hold off;axis([0 600 0 1400]);camroll(180);colorbar;caxis([0 1]);
set(gca,'Visible','off')
print(Fighandle,strcat('D:\Pictures\processed\Flow\Basic\Figure\Graph Theory\','Participation.tif'),'-dtiff','-r0');

Fighandle=figure;
graph_temp=graph(Correlation_bin_prop);
set(Fighandle, 'Position', [10,10, 600, 1400]);
plot(Zbrain_brainMask2D(:,2),1400-Zbrain_brainMask2D(:,1),'k');hold on;
scatter(Node_all(:,2),Node_all(:,1),(degree(graph_temp)/5)+20,'k','filled');hold on;
LWidths = 1*graph_temp.Edges.Weight/max(graph_temp.Edges.Weight);
LWidths(~isfinite(LWidths))=0.01;
LColors = ones(length(LWidths),3);LColors=LColors-LWidths;
plot(graph_temp,'EdgeColor',LColors,'NodeColor',color_temp,'LineStyle','-','LineWidth',LWidths,'Marker','o','MarkerSize',degree(graph_temp)/5,'XData',Node_all(:,2),'YData',Node_all(:,1),'NodeLabel',char.empty(100,0));hold off;axis([0 600 0 1400]);camroll(180);colorbar;
set(gca,'Visible','off')
print(Fighandle,strcat('D:\Pictures\processed\Flow\Basic\Figure\Graph Theory\','25PercGraph.tif'),'-dtiff','-r0');