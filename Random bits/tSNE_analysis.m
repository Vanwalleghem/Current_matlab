tSNE_Flow=[];tSNE_Flow = tsne_perp(ZS,[],2,100,10:10:100);

Fighandle=figure;
set(Fighandle, 'Position', [100, 100, 1400, 900]);
xplot=floor(sqrt(length(tSNE_Flow)));yplot=ceil(length(tSNE_Flow)/xplot);IDX=[];counter=1;
for i=1:length(tSNE_Flow)
    %[IDX{i}, isnoise]=DBSCAN(tSNE_Flow{i},10,20);
    subplot(xplot,yplot,i);gscatter(tSNE_Flow{i}(:,1),tSNE_Flow{i}(:,2))%,IDX{i})    
    counter=counter+1;
end

Fighandle=figure;
set(Fighandle, 'Position', [100, 100, 1400, 900]);
xplot=floor(sqrt(30));yplot=ceil(30/xplot);
for i=1:30
    [IDX_temp{i}, isnoise]=DBSCAN(tSNE_Flow{1},i,20);
    subplot(xplot,yplot,i);gscatter(tSNE_Flow{1}(:,1),tSNE_Flow{2}(:,2),IDX_temp{i})
end

Fighandle=figure;
set(Fighandle, 'Position', [100, 100, 1400, 900]);
xplot=floor(sqrt(30));yplot=ceil(30/xplot);
for i=1:30
    [IDX_min{i}, isnoise]=DBSCAN(tSNE_Flow{1},12,i);
    subplot(xplot,yplot,i);gscatter(tSNE_Flow{1}(:,1),tSNE_Flow{2}(:,2),IDX_min{i},[],[],[],'off')
end

for i=1:30
    temp(i)=max(IDX_min{i});    
end
figure;plot(temp);

Fighandle=figure;
set(Fighandle, 'Position', [100, 100, 1400, 900]);
xplot=floor(sqrt(30));yplot=ceil(30/xplot);
for i=1:30    
    subplot(xplot,yplot,i);gscatter(tSNE_Flow{2}(:,1),tSNE_Flow{2}(:,2),IDX_temp{i},[],[],[],'off')
end

Fighandle=figure;
set(Fighandle, 'Position', [100, 100, 1400, 900]);
gscatter(tSNE_Flow{1}(:,1),tSNE_Flow{1}(:,2),IDX_min{24})


DBClusters=zeros(max(IDX_temp{8}),size(ZS,2));
for i=1:max(IDX_temp{8})
    idx=find(IDX_temp{8}==i);
    DBClusters(i,:)=mean(ZS2(idx,:),1);
end
figure;imagesc(DBClusters);


Fighandle=figure;
set(Fighandle, 'Position', [100, 100, 1400, 900]);
xplot=floor(sqrt(length(tSNE_Flow)));yplot=ceil(length(tSNE_Flow)/xplot);
for i=1:length(tSNE_Flow)    
    subplot(xplot,yplot,i);gscatter(tSNE_Flow{i}(:,1),tSNE_Flow{i}(:,2),IDX{i},[],[],[],'off')    
    counter=counter+1;
end