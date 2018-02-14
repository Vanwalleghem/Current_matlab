figure;
counter=1;xplot=length(ItiaList);yplot=4;counter2=1;
Cluster{1,1}=[2,3];Cluster{1,2}=[1];
Cluster{2,1}=[1,3];Cluster{2,2}=[2];
Cluster{3,1}=[1];Cluster{3,2}=[2];
Cluster{4,1}=[1];
Cluster{5,1}=[2 3];Cluster{5,2}=[1];
Cluster{6,1}=[2 3];Cluster{6,2}=[1];
Cluster{7,1}=[1];Cluster{7,2}=[2];
Cluster{8,1}=[1];
Cluster{9,1}=[2 3];Cluster{9,2}=[1];Cluster{9,3}=[4];
Cluster{10,1}=[1];

colors=[0 1 0;1 0 0; 1 0 0.5];

for i=1:length(ItiaList)
    regionName=ItiaList{i};
    idx_temp=LinReg.(regionName).KmeansIdx_select;
    idx_modif=zeros(size(idx_temp));
    GoodBet_temp=LinReg.(regionName).GoodBetas;
    LinReg.(regionName).GB_select=[1,2,3];
    for j=1:3
        for k=1:length(Cluster{i,j})
            if k
                idx_modif(idx_temp==GoodBet_temp(Cluster{i,j}(k)))=j;
            end
        end
    end
    LinReg.(regionName).idx_sorted=idx_modif;    
end

figure;
counter=1;xplot=length(ItiaList);yplot=3;counter2=1;
for i=1:length(ItiaList)
    regionName=ItiaList{i};
    idx_temp=LinReg.(regionName).idx_sorted;
    GoodBet_temp=LinReg.(regionName).GB_select;
    counter=counter2;
    for j=1:3
        idx_temp2=find(idx_temp==GoodBet_temp(j));
        if idx_temp2
            subplot(xplot,yplot,counter+(j-1));plot(mean(LinReg.(regionName).ZS_AVG(idx_temp2,:),1));title(num2str(length(idx_temp2)));
        end
    end
    counter2=counter2+yplot;
end


%Building the AVG across 3 presentations
for j=1:length(ItiaList)
    regionName=ItiaList{j};
    ZS2=LinReg.(regionName).ZS_rsq;
    ZS_AVG2=zeros(size(ZS2,1),246);
    parfor idx_ZS=1:size(ZS2,1)
        start=30;
        AVG=[];
        for i=1:3
            AVG(i,:)=ZS2(idx_ZS,start:start+40);
            start=start+40;
        end
        AVG=mean(AVG,1);
        AVG=AVG-mean(AVG(1:10));%-min(AVG);
        temp=[];
        for j=2:6
            for i=1:3
                temp(i,:)=ZS2(idx_ZS,start:start+40);
                start=start+40;
            end
            temp=mean(temp,1);
            temp=temp-mean(temp(1:10));%-min(temp);
            AVG=[AVG temp];
        end
        ZS_AVG2(idx_ZS,:)=AVG;
    end
    LinReg.(regionName).ZS_AVG=ZS_AVG2;
end

%Kmeans of AVG
options = statset('UseParallel',1); 
for j=[2 5 8]%1:length(ItiaList)
    regionName=ItiaList{j};    
    if j==6 | j==10
        [idxKmeans Cmap]=kmeans(LinReg.(regionName).ZS_AVG,15,'Distance','cityblock','Options',options,'Replicates',3,'MaxIter',1000,'Display','final');
    elseif j==2 | j==5 | j==8        
        [idxKmeans Cmap]=kmeans(LinReg.(regionName).ZS_AVG,5,'Distance','cityblock','Options',options,'Replicates',3,'MaxIter',1000,'Display','final');
    else
        [idxKmeans Cmap]=kmeans(LinReg.(regionName).ZS_AVG,5,'Distance','cityblock','Options',options,'Replicates',3,'MaxIter',1000,'Display','final');
    end
    LinReg.(regionName).KmeansCenter_AVG=Cmap;
    LinReg.(regionName).KmeansIdx_AVG=idxKmeans;
end

%Get the good ones
for j=1:[2 5 8]%length(ItiaList)
    regionName=ItiaList{j};
    [~,LinReg.(regionName).GoodBetas_AVG]=Test_Regress( LinReg.(regionName).KmeansCenter_AVG,Stimuli_AVG,LinReg.(regionName).KmeansIdx_AVG,0.5);
end

figure;
counter=1;xplot=length(ItiaList);yplot=6;counter2=1;
for i=1:length(ItiaList)
    regionName=ItiaList{i};
    idx_temp=LinReg.(regionName).KmeansIdx_AVG;
    GoodBet_temp=LinReg.(regionName).GoodBetas_AVG;
    counter=counter2;
    for j=1:length(GoodBet_temp)
        idx_temp2=find(idx_temp==GoodBet_temp(j));
        subplot(xplot,yplot,counter+(j-1));plot(mean(LinReg.(regionName).ZS_AVG(idx_temp2,:),1));title(num2str(length(idx_temp2)));
    end
    counter2=counter2+yplot;
end


%Correlate back to Kmeans to remove crap
Threshold=0.5;
for i=1:length(ItiaList)
    regionName=ItiaList{i};
    temp_ZS=LinReg.(regionName).ZS_AVG;
    idx_temp=LinReg.(regionName).KmeansIdx_AVG;
    GoodBet_temp=LinReg.(regionName).GoodBetas_AVG;
    for j=1:length(GoodBet_temp)
        idx_g=find(idx_temp==GoodBet_temp(j));
        temp_g=temp_ZS(idx_g,:);
        corr_temp=zeros(size(idx_g));
        for k=1:length(idx_g)        
            temp=corrcoef(LinReg.(regionName).KmeansCenter_AVG(GoodBet_temp(j),:), temp_g(k,:));
            corr_temp(k)=temp(1,2);
        end        
        idx_temp(idx_g(find(corr_temp<=Threshold)))=0;
    end
    LinReg.(regionName).KmeansIdx_select_AVG=idx_temp;
end
clearvars i j k temp ans Cmap coeff coefficients explained i Hindbrain_Mask idx_temp idxKmeans ish_nb

Threshold=0.3;
i=3;
regionName=ItiaList{i};
temp_ZS=LinReg.(regionName).ZS_AVG;
idx_temp=LinReg.(regionName).KmeansIdx_AVG;
GoodBet_temp=LinReg.(regionName).GoodBetas_AVG;
for j=1:length(GoodBet_temp)
    idx_g=find(idx_temp==GoodBet_temp(j));
    temp_g=temp_ZS(idx_g,:);
    corr_temp=zeros(size(idx_g));
    for k=1:length(idx_g)
        temp=corrcoef(LinReg.(regionName).KmeansCenter_AVG(GoodBet_temp(j),:), temp_g(k,:));
        corr_temp(k)=temp(1,2);
    end
    idx_temp(idx_g(find(corr_temp<=Threshold)))=0;    
end
LinReg.(regionName).KmeansIdx_select_AVG=idx_temp;

clearvars i j k temp ans Cmap coeff coefficients explained i Hindbrain_Mask idx_temp idxKmeans ish_nb


Threshold=0.5;
i=1;
regionName=ItiaList{i};
temp_ZS=LinReg.(regionName).ZS_AVG;
idx_temp=LinReg.(regionName).KmeansIdx_AVG;
GoodBet_temp=LinReg.(regionName).GoodBetas_AVG;
idx_g=find(idx_temp==GoodBet_temp(2) | idx_temp==GoodBet_temp(3));
temp_g=temp_ZS(idx_g,:);
corr_temp=zeros(size(idx_g));
for k=1:length(idx_g)
    temp=corrcoef(LinReg.(regionName).KmeansCenter_AVG(GoodBet_temp(2),:), temp_g(k,:));
    corr_temp(k)=temp(1,2);
end
idx_temp(idx_g(find(corr_temp<=Threshold)))=0;
idx_temp(idx_g(find(corr_temp>Threshold)))=GoodBet_temp(2);
LinReg.(regionName).KmeansIdx_select_AVG=idx_temp;

i=3;
regionName=ItiaList{i};
temp_ZS=LinReg.(regionName).ZS_AVG;
idx_temp=LinReg.(regionName).KmeansIdx_AVG;
GoodBet_temp=LinReg.(regionName).GoodBetas_AVG;
idx_g=find(idx_temp==GoodBet_temp(1) | idx_temp==GoodBet_temp(2));
temp_g=temp_ZS(idx_g,:);
corr_temp=zeros(size(idx_g));
for k=1:length(idx_g)
    temp=corrcoef(LinReg.(regionName).KmeansCenter_AVG(GoodBet_temp(1),:), temp_g(k,:));
    corr_temp(k)=temp(1,2);
end
idx_temp(idx_g(find(corr_temp<=Threshold)))=0;
idx_temp(idx_g(find(corr_temp>Threshold)))=GoodBet_temp(1);
LinReg.(regionName).KmeansIdx_select_AVG=idx_temp;

i=6;
regionName=ItiaList{i};
temp_ZS=LinReg.(regionName).ZS_AVG;
idx_temp=LinReg.(regionName).KmeansIdx_AVG;
GoodBet_temp=LinReg.(regionName).GoodBetas_AVG;
idx_g=find(idx_temp==GoodBet_temp(1) | idx_temp==GoodBet_temp(2) | idx_temp==GoodBet_temp(3) | idx_temp==GoodBet_temp(4));
temp_g=temp_ZS(idx_g,:);
corr_temp=zeros(size(idx_g));
for k=1:length(idx_g)
    temp=corrcoef(LinReg.(regionName).KmeansCenter_AVG(GoodBet_temp(1),:), temp_g(k,:));
    corr_temp(k)=temp(1,2);
end
idx_temp(idx_g(find(corr_temp<=Threshold)))=0;
idx_temp(idx_g(find(corr_temp>Threshold)))=GoodBet_temp(4);
LinReg.(regionName).KmeansIdx_select_AVG=idx_temp;

i=9;
regionName=ItiaList{i};
temp_ZS=LinReg.(regionName).ZS_AVG;
idx_temp=LinReg.(regionName).KmeansIdx_AVG;
GoodBet_temp=LinReg.(regionName).GoodBetas_AVG;
idx_g=find(idx_temp==GoodBet_temp(1) | idx_temp==GoodBet_temp(2));
temp_g=temp_ZS(idx_g,:);
corr_temp=zeros(size(idx_g));
for k=1:length(idx_g)
    temp=corrcoef(LinReg.(regionName).KmeansCenter_AVG(GoodBet_temp(2),:), temp_g(k,:));
    corr_temp(k)=temp(1,2);
end
idx_temp(idx_g(find(corr_temp<=Threshold)))=0;
idx_temp(idx_g(find(corr_temp>Threshold)))=GoodBet_temp(2);
LinReg.(regionName).KmeansIdx_select_AVG=idx_temp;

i=10;
regionName=ItiaList{i};
temp_ZS=LinReg.(regionName).ZS_AVG;
idx_temp=LinReg.(regionName).KmeansIdx_AVG;
GoodBet_temp=LinReg.(regionName).GoodBetas_AVG;
idx_g=find(idx_temp==GoodBet_temp(2) | idx_temp==GoodBet_temp(4) | idx_temp==GoodBet_temp(6)  | idx_temp==GoodBet_temp(1));
temp_g=temp_ZS(idx_g,:);
corr_temp=zeros(size(idx_g));
for k=1:length(idx_g)
    temp=corrcoef(LinReg.(regionName).KmeansCenter_AVG(GoodBet_temp(4),:), temp_g(k,:));
    corr_temp(k)=temp(1,2);
end
idx_temp(idx_g(find(corr_temp<=Threshold)))=0;
idx_temp(idx_g(find(corr_temp>Threshold)))=GoodBet_temp(4);
LinReg.(regionName).KmeansIdx_select_AVG=idx_temp;

figure;
counter=1;xplot=length(ItiaList);yplot=6;counter2=1;
for i=1:length(ItiaList)
    regionName=ItiaList{i};
    idx_temp=LinReg.(regionName).KmeansIdx_select_AVG;
    GoodBet_temp=LinReg.(regionName).GoodBetas_AVG;
    counter=counter2;
    for j=1:length(GoodBet_temp)
        idx_temp2=find(idx_temp==GoodBet_temp(j));
        subplot(xplot,yplot,counter+(j-1));plot(mean(LinReg.(regionName).ZS_AVG(idx_temp2,:),1));title(num2str(length(idx_temp2)));
    end
    counter2=counter2+yplot;
end

colors={};
i=1;colors{i}=[0 0.8 0; 0.8 0 0.8];regionName=ItiaList{i};LinReg.(regionName).GoodBetas_AVG_final=LinReg.(regionName).GoodBetas_AVG([1 2]);
i=2;colors{i}=[0 0.8 0; 0.8 0 0.8];regionName=ItiaList{i};LinReg.(regionName).GoodBetas_AVG_final=LinReg.(regionName).GoodBetas_AVG([1 2]);
i=3;colors{i}=[0 0.8 0; 0.8 0 0.8];regionName=ItiaList{i};LinReg.(regionName).GoodBetas_AVG_final=LinReg.(regionName).GoodBetas_AVG([3 1]);
i=4;colors{i}=[0 0.8 0];regionName=ItiaList{i};LinReg.(regionName).GoodBetas_AVG_final=LinReg.(regionName).GoodBetas_AVG([1]);
i=5;colors{i}=[0.8 0 0.8;0.2 1 0.8];regionName=ItiaList{i};LinReg.(regionName).GoodBetas_AVG_final=LinReg.(regionName).GoodBetas_AVG([1 2]);
%i=6;colors{i}=[0 0.8 0; 0.8 0 0.8; 0.2 1 0.8];regionName=ItiaList{i};LinReg.(regionName).GoodBetas_AVG_final=LinReg.(regionName).GoodBetas_AVG([1 5 4]);
i=6;colors{i}=[0 0.8 0; 0.8 0 0.8; 0.2 1 0.8];regionName=ItiaList{i};LinReg.(regionName).GoodBetas_AVG_final=LinReg.(regionName).GoodBetas_AVG([1 5 4]);
i=7;colors{i}=[0 0.8 0; 0.8 0 0.8];regionName=ItiaList{i};LinReg.(regionName).GoodBetas_AVG_final=LinReg.(regionName).GoodBetas_AVG([1 2]);
i=8;colors{i}=[0 0.8 0];regionName=ItiaList{i};LinReg.(regionName).GoodBetas_AVG_final=LinReg.(regionName).GoodBetas_AVG([1]);
i=9;colors{i}=[0 0.8 0; 0.8 0 0.8;0.2 1 0.8];regionName=ItiaList{i};LinReg.(regionName).GoodBetas_AVG_final=LinReg.(regionName).GoodBetas_AVG([3 2 1]);
i=10;colors{i}=[0 0.8 0; 0.8 0 0.8];regionName=ItiaList{i};LinReg.(regionName).GoodBetas_AVG_final=LinReg.(regionName).GoodBetas_AVG([4 3]);

figure;
counter=1;xplot=length(ItiaList);yplot=3;counter2=1;
for i=1:length(ItiaList)
    regionName=ItiaList{i};
    idx_temp=LinReg.(regionName).KmeansIdx_select_AVG;
    GoodBet_temp=LinReg.(regionName).GoodBetas_AVG_final;
    counter=counter2;
    if i==5
        counter=counter+1;
    end
    for j=1:length(GoodBet_temp)
        idx_temp2=find(idx_temp==GoodBet_temp(j));        
        subplot(xplot,yplot,counter+(j-1));plot(mean(LinReg.(regionName).ZS_AVG(idx_temp2,:),1),'color',colors{i}(j,:));title(num2str(length(idx_temp2)));
    end
    counter2=counter2+yplot;
end

