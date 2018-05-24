load('D:\temp\allfisheyes.mat')
idx_ELO=find(sum(ismember(namefish,'ELO'),2)==3);
ZS_eyes=zscore(matrixeyes,1,2);

x = linspace(0.2,size(ZS_eyes,2)/5,size(ZS_eyes,2));
Fighandle=figure;
set(Fighandle, 'Position', [0, 0, 1280, 1024]);
set(findall(Fighandle,'type','text'),'fontSize',12,'fontWeight','bold','FontName','Arial')
xplot=2;yplot=2;
temp_ELO=ZS_eyes(idx_ELO,:);
temp_ERO=ZS_eyes(find(sum(ismember(namefish,'ERO'),2)==3),:);
subplot(xplot,yplot,1);
imagesc(temp_ELO,[-5 5]);
subplot(xplot,yplot,2);
temp=mean(temp_ELO,1);std_temp=std(temp_ELO,1,1);
shadedErrorBar(x, temp, std_temp);axis([0 150 -6 6]);
subplot(xplot,yplot,3);
imagesc(temp_ERO,[-5 5]);
subplot(xplot,yplot,4);
temp=mean(temp_ERO,1);std_temp=std(temp_ERO,1,1);
shadedErrorBar(x, temp, std_temp);axis([0 150 -6 6]);

Fighandle=figure;
set(Fighandle, 'Position', [0, 0, 1280, 1024]);
set(findall(Fighandle,'type','text'),'fontSize',12,'fontWeight','bold','FontName','Arial')
xplot=2;yplot=2;
temp_ELO=matrixeyes(idx_ELO([1 3 4 5 6]),:);
temp_ERO=matrixeyes(find(sum(ismember(namefish,'ERO'),2)==3),:);
subplot(xplot,yplot,1);
imagesc(temp_ELO,[-5 5]);
subplot(xplot,yplot,2);
temp=mean(temp_ELO,1);std_temp=std(temp_ELO,1,1);
shadedErrorBar(x, temp, std_temp);axis([0 150 -6 6]);
subplot(xplot,yplot,3);
imagesc(temp_ERO,[-5 5]);
subplot(xplot,yplot,4);
temp=mean(temp_ERO,1);std_temp=std(temp_ERO,1,1);
shadedErrorBar(x, temp, std_temp);axis([0 150 -6 6]);