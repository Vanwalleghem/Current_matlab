MatFiles=dir('*_output_analysis_matlab.mat');
name=strcat(MatFiles(1).name);
Calcium=load(name, 'DenoisedTraces');
Calcium=Calcium.DenoisedTraces;
MatFiles(1).number=size(Calcium,1);
%Spikes=load(name, 'Spikes');
%Spikes=Spikes.Spikes;
Noise=load(name, 'Noise');
Noise=Noise.Noise;
Fitness=load(name, 'idx_components');
Fitness=Fitness.idx_components+1;
GoodCalcium=Calcium(Fitness,:);%should be (Fitness,:)
%GoodSpikes=Spikes(Fitness,:);
MatFiles(1).GoodNumber=length(Fitness);
MatFiles(1).GC=GoodCalcium;
MatFiles(1).GN=Noise(Fitness,:);
%Cn=load(name, 'Correlation_image');
%MatFiles(1).Cn=Cn.Correlation_image;
for i = 2:length(MatFiles)
name=strcat(MatFiles(i).name);
C=load(name, 'DenoisedTraces');
C=C.DenoisedTraces;
%     if i==3
%         C=[C(:,1) C(:,1) C(:,1:58)];
%     end
%S=load(name, 'Spikes');
%S=S.Spikes;
N=load(name, 'Noise');
N=N.Noise;
F=load(name, 'idx_components');
F=F.idx_components+1;
GC=C(F,:);
GN=N(F,:);
%GS=S(F,:);
%Cn=load(name, 'Correlation_image');
%Noise=vertcat(Noise,N);
%Calcium=vertcat(Calcium,C);
%Spikes=vertcat(Spikes,S);
%Fitness=horzcat(Fitness,F);
%GoodCalcium=vertcat(GoodCalcium,GC);
%GoodSpikes=vertcat(GoodSpikes,GS);
MatFiles(i).number=size(Calcium,1);
MatFiles(i).GoodNumber=MatFiles(i-1).GoodNumber+length(F);
MatFiles(i).GC=GC;
MatFiles(i).GN=GN;
%MatFiles(i).Cn=Cn.Correlation_image;
end
clearvars GC C S F N name i GS DF Cn Calcium GN Noise Fitness;

i=1;GoodCalcium=MatFiles(i).GC;
for i=2:length(MatFiles)	
    GoodCalcium=vertcat(GoodCalcium,MatFiles(i).GC);
end

ZS=zscore(GoodCalcium,1,2);

% ZS_diff = zeros(size(ZS));
% parfor i=1:size(ZS,1)
%     temp = TVRegDiff(ZS(i,:), 50, 1e-1,[],'large',1e-8,[],0,0);
%     ZS_diff(i,:)=temp;
% end

options = statset('UseParallel',1); [idxKmeans_ZS Cmap_ZS]=kmeans(ZS,30,'Options',options,'Distance','cityblock','Replicates',3,'MaxIter',1000,'Display','final');
%options = statset('UseParallel',1); [idxKmeans_ZS2 Cmap_ZS2]=kmeans(ZS,5,'Options',options,'Distance','cityblock','Replicates',5,'MaxIter',1000,'Display','final');
% options = statset('UseParallel',1); [idxKmeans_ZS_diff Cmap_ZS_diff]=kmeans(ZS_diff,20,'Options',options,'Distance','cityblock','Replicates',5,'MaxIter',1000,'Display','final');

Stimuli=zeros(1,size(GoodCalcium,2));
GCaMP6=[0,1.69644104899772,5.13796058542217,8.27886020152244,10.3756715204800,11.8173714529814,12.2425184714093,10.8571417354877,8.80831829681196,6.91339112244670,5.46959264663869,4.30868766622567,3.42533619066766,2.75378443486879,2.18017250852183,1.72816235135824,1.32732537295463,1.00684435500268,0.730210038304555,0.530242444093118,0.362253250339685,0.227668255288566,0.0869242416152502,0.000718266708050853,-0.0828334873368325]';
idxStart=60;
for i=1:5
    Stimuli(1,(idxStart+(i-1)*40):(idxStart+(i-1)*40)+size(GCaMP6,1)-1)=GCaMP6;
end

Correlation=[];
parfor i=1:size(ZS,1)
    temp=corrcoef(Stimuli,ZS(i,:));    
    Correlation(i)=temp(1,2);
end

figure;
for i=[0.2:0.1:1]
    plot(mean(ZS(find(Correlation>i),:),1));hold on;
end

figure;
for i=[-0.2:-0.1:-1]
    plot(mean(ZS(find(Correlation<i),:),1));hold on;
end

figure;imagesc(ZS(find(Correlation>0.2),:));


rsquare_loom=[ModelResults.rsquared];
idx_rsq=find(rsquare_loom>0.2 & rsquare_loom<1);
figure;
imagesc(ZS(idx_rsq,:));

figure;
for i=1:size(Cmap_ZS,1)
    plot(Cmap_ZS(i,:));pause
end

Numbers=[1 [MatFiles.GoodNumber]];
counter=1;
idx_Plane=nan(length(ZS),1);
idx_Position=nan(length(ZS),1);
idx_Fish=nan(length(ZS),1);
name=strcat(MatFiles(1).name);
for i=1:length(MatFiles)	
    name=strcat(MatFiles(i).name);
    [Plane,~]=regexp(name,'_(\d+)_','tokens','match');Plane=str2num(Plane{1}{1});
    if strfind(name,'CLO')
        Fish=1;
    elseif strfind(name,'ELO')
        Fish=2;
    elseif strfind(name,'ERO')
        Fish=3;
    else
        Fish=4;
    end
    idx_Plane(Numbers(i):Numbers(i+1))=Plane;
    idx_Position(Numbers(i):Numbers(i+1))=Fish;
    [Fish,~]=regexp(name,'Fish201709(\d+)_','tokens','match');Fish=str2num(Fish{1}{1});
    idx_Fish(Numbers(i):Numbers(i+1))=Fish;
end
clearvars i Fish Plane name counter

GoodBetas=21;
rows=length(GoodBetas);
counter=1;
Fighandle=figure;
set(Fighandle, 'Position', [100, 100, 1300, 900]);
for i=GoodBetas
    idx_temp=find(idxKmeans_ZS==i);
    subplot(rows,5,counter);plot(mean(ZS(idx_temp,:),1));
    subplot(rows,5,counter+1);imagesc(ZS(idx_temp,:),[0 3]);
    subplot(rows,5,counter+2);histogram(idx_Plane(idx_temp));
    subplot(rows,5,counter+3);histogram(idx_Position(idx_temp),[0.5:1:4.5]);
    subplot(rows,5,counter+4);histogram(idx_Fish(idx_temp));
    counter=counter+5;
end

figure;histogram(idx_Position);

counter=1;
Fighandle=figure;
set(Fighandle, 'Position', [100, 100, 1300, 900]);
for i=GoodBetas
    idx_temp=find(idxKmeans_ZS==i);
    idx_temp_ELO=find(idx_Position(idx_temp)==1);
    idx_temp_ILO=find(idx_Position(idx_temp)==2);
    subplot(rows,1,counter);plot(mean(ZS(idx_temp(idx_temp_ELO),:),1));ylim([-1.5 3.5])
    hold on;plot(mean(ZS(idx_temp(idx_temp_ILO),:),1)); legend('ELO','ILO')
    counter=counter+1;
end


All_ROIs=[];
ROIs_idx=[];
i=1;
name=strcat(MatFiles(i).name);
Rs=load(name, 'ROIs');
Rs=Rs.ROIs;
F=load(name, 'idx_components');
F=F.idx_components+1;
Rs=Rs(:,F);
All_ROIs{1}=Rs;
ROIs_idx(1)=size(Rs,2);
for i = 2:length(MatFiles)
    name=strcat(MatFiles(i).name);
    Rs=load(name, 'ROIs');
    Rs=Rs.ROIs;
    F=load(name, 'idx_components');
    F=F.idx_components+1;
    Rs=Rs(:,F);
    All_ROIs{i}=Rs;
    ROIs_idx(i)=ROIs_idx(i-1)+length(F);    
end
clearvars GC C S F N name i;

Numbers=[0 [ROIs_idx]];
temp=[];
counter=1;
for i=GoodBetas
    temp{counter}=find(idxKmeans_ZS==i);    
    counter=counter+1;    
end

Start=min(cellfun(@min, temp));Start=find(Numbers<Start,1,'last');
filename=MatFiles(Start).name;

%colors = [1,0,0;0,1,0;0,0,1;1,0.600000000000000,0;1,0,0.7];
colors = distinguishable_colors(length(GoodBetas),[1 1 1; 0 0 0]);
% colors = [0         0    1.0000
%          0    0.5000    1.0000
%     1.0000         0         0
%     1.0000    0.1034    0.7241
%     1.0000    0.5000    0.3000
%          0    0.7000    0.2000
%     0.5000    0.5000         0
%          0    0.5000    0.5000];
colors = colors*256;
for idx=Start:length(MatFiles)
    filename=MatFiles(idx).name;
    ROIsNb=[];ClusterNb=[];
    %for k = 1 : length(temp)
    for k = 1 : length(temp)
        tempROIsNb=find([temp{k}]<=Numbers(idx+1));
        if tempROIsNb            
            ROIsNb=[ROIsNb; temp{k}(tempROIsNb)];
            temp{k}(tempROIsNb)=[];
            ClusterNb=[ClusterNb ; repmat(k,length(tempROIsNb),1)];
        end
    end
    if ROIsNb
        imagename=regexp(filename,'_output','split');
        %imagename=regexp(imagename,'_output_analysis_matlab2.mat','split');
        imagename=strcat(imagename{1},'_mean.tif');
        image=double(imread(imagename));image=image/max(max(image));image=image*128;
        image=uint8(image);
        image2=zeros(size(image(:,:,1)));
        image3=repmat(image,1,1,3);
        ROIs=All_ROIs{idx};       
        ROIsNb=ROIsNb-Numbers(idx);
        ROIs=ROIs(:,ROIsNb);
        for k = 1 : size(ROIs,2)
            image2=zeros(size(image(:,:,1)));
            ROI=full(reshape(ROIs(:,k),size(image,1),size(image,2)));            
            image2=image2+ROI;image2=(image2/max(max(image2)));image2=uint8(image2);
            for j=1:3
                image3(:,:,j)=image3(:,:,j)+image2*colors(ClusterNb(k),j);
            end
        end
        %image3(:,:,3)=image;
        name=strcat('_Kmeans_',imagename(4:end));
    imwrite(image3,name,'tif');
    end
    %image3=uint8(image3);

end
clearvars idx i temp tempFileNb fileNb AVG_files filename image counter Numbers image2 image3 k ROI ROIs ROIsNb Start tempROIsNb name imagename tempidx Raster



Fighandle=figure;
set(Fighandle, 'Position', [100, 100, 1400, 900]);x = linspace(0.2,size(ZS,2)/5,size(ZS,2));
counter=1;counter2=1;xplot=floor(sqrt(length(GoodBetas)));yplot=ceil(length(GoodBetas)/xplot);
for i=GoodBetas
%     if counter==3
%         counter=counter+1;
%     end
    idx_temp=find(idxKmeans_ZS==i);    
    subplot(xplot,yplot,counter);plot(x,mean(ZS(idx_temp,:),1),'color',colors(counter2,:)/256);%axis([0 131 -1 4]);rectangle('FaceColor','r','Position',[11 -1 10 0.25]);rectangle('FaceColor','r','Position',[51 -1 10 0.25]);rectangle('FaceColor','r','Position',[91 -1 10 0.25]);rectangle('FaceColor','b','Position',[31 -1 10 0.25]);rectangle('FaceColor','b','Position',[71 -1 10 0.25]);rectangle('FaceColor','b','Position',[111 -1 10 0.25]);    
    counter=counter+1;
    counter2=counter2+1;
end


GoodBetas_diff=[4 5 8 12 24 28];
rows=length(GoodBetas_diff);
counter=1;
Fighandle=figure;
set(Fighandle, 'Position', [100, 100, 1300, 900]);
for i=GoodBetas_diff
    idx_temp=find(idxKmeans_diff==i);
    subplot(rows,4,counter);plot(mean(ZS(idx_temp,:),1));
    subplot(rows,4,counter+1);imagesc(ZS(idx_temp,:),[0 3]);
    subplot(rows,4,counter+2);histogram(idx_Plane(idx_temp));
    subplot(rows,4,counter+3);histogram(idx_Position(idx_temp));
    counter=counter+4;
end

Numbers=[0 [ROIs_idx]];
temp=[];
counter=1;
for i=GoodBetas_diff
    temp{counter}=find(idxKmeans_diff==i);    
    counter=counter+1;    
end

Start=min(cellfun(@min, temp));Start=find(Numbers<Start,1,'last');
filename=MatFiles(Start).name;

%colors = [1,0,0;0,1,0;0,0,1;1,0.600000000000000,0;1,0,0.7];
colors = distinguishable_colors(length(GoodBetas),[1 1 1; 0 0 0]);
% colors = [0         0    1.0000
%          0    0.5000    1.0000
%     1.0000         0         0
%     1.0000    0.1034    0.7241
%     1.0000    0.5000    0.3000
%          0    0.7000    0.2000
%     0.5000    0.5000         0
%          0    0.5000    0.5000];
colors = colors*256;
for idx=Start:length(MatFiles)
    filename=MatFiles(idx).name;
    ROIsNb=[];ClusterNb=[];
    %for k = 1 : length(temp)
    for k = 1 : length(temp)
        tempROIsNb=find([temp{k}]<=Numbers(idx+1));
        if tempROIsNb            
            ROIsNb=[ROIsNb; temp{k}(tempROIsNb)];
            temp{k}(tempROIsNb)=[];
            ClusterNb=[ClusterNb ; repmat(k,length(tempROIsNb),1)];
        end
    end
    if ROIsNb
        imagename=regexp(filename,'_output_analysis','split');
        %imagename=regexp(imagename,'_output_analysis_matlab2.mat','split');
        imagename=strcat('AVG_',imagename{1},'.tif');
        image=double(imread(imagename));image=image/max(max(image));image=image*128;
        image=uint8(image);
        image2=zeros(size(image(:,:,1)));
        image3=repmat(image,1,1,3);
        ROIs=All_ROIs{idx};       
        ROIsNb=ROIsNb-Numbers(idx);
        ROIs=ROIs(:,ROIsNb);
        for k = 1 : size(ROIs,2)
            image2=zeros(size(image(:,:,1)));
            ROI=full(reshape(ROIs(:,k),size(image,1),size(image,2)));            
            image2=image2+ROI;image2=(image2/max(max(image2)));image2=uint8(image2);
            for j=1:3
                image3(:,:,j)=image3(:,:,j)+image2*colors(ClusterNb(k),j);
            end
        end
        %image3(:,:,3)=image;
        name=strcat('_Kmeans_diff_',imagename(4:end));
    imwrite(image3,name,'tif');
    end
    %image3=uint8(image3);

end
clearvars idx i temp tempFileNb fileNb AVG_files filename image counter Numbers image2 image3 k ROI ROIs ROIsNb Start tempROIsNb name imagename tempidx Raster

Fighandle=figure;
set(Fighandle, 'Position', [100, 100, 1400, 900]);x = linspace(0.2,size(ZS,2)/5,size(ZS,2));
counter=1;counter2=1;xplot=floor(sqrt(length(GoodBetas_diff)));yplot=ceil(length(GoodBetas_diff)/xplot);
for i=GoodBetas_diff
%     if counter==3
%         counter=counter+1;
%     end
    idx_temp=find(idxKmeans_diff==i);    
    subplot(xplot,yplot,counter);plot(x,mean(ZS(idx_temp,:),1),'color',colors(counter2,:)/256);%axis([0 131 -1 4]);rectangle('FaceColor','r','Position',[11 -1 10 0.25]);rectangle('FaceColor','r','Position',[51 -1 10 0.25]);rectangle('FaceColor','r','Position',[91 -1 10 0.25]);rectangle('FaceColor','b','Position',[31 -1 10 0.25]);rectangle('FaceColor','b','Position',[71 -1 10 0.25]);rectangle('FaceColor','b','Position',[111 -1 10 0.25]);    
    counter=counter+1;
    counter2=counter2+1;
end

Opt_trap=zeros(2,size(ZS,2));
ON=[10 45 80];
OFF=[13 48 83];
GCaMP6=[0.267549369506403,1.55457586790809,2.27073744424031,2.38555016082230,2.18676492016293,1.82043261673794,1.43235606073233,1.11018612933621,0.851746802276788,0.630647114302851,0.462591970408579,0.333697437213782,0.242186236930119,0.183255706836320,0.129716047724292,0.0837715067907188,0.0524466692619288,0.0313691326187021,0.0173098198493843,0.00749011066138279,0]';
for i=1:length(ON)
    Opt_trap(1,ON(i):ON(i)+size(GCaMP6,1)-1)=GCaMP6';
    Opt_trap(2,OFF(i):OFF(i)+size(GCaMP6,1)-1)=GCaMP6';
end
clearvars GCaMP6 back back_off fwd fwd_off GCaMP6s;
Opt_trap=Opt_trap(:,1:100);

[Model_ZS,GoodBetas_ZS]=Test_Regress(Cmap_ZS,Opt_trap,idxKmeans_ZS,0.2);

ModelResultsSeg_ZS=[];
parfor i=1:length(ZS)
    %mdl=stepwiselm(Regressor',DataSet(i,:),'linear','Criterion','adjrsquared','Intercept',false,'Upper','interactions','Verbose',0);
    mdl=stepwiselm(Opt_trap',ZS(i,:),'linear','Criterion','adjrsquared','Upper','linear','Verbose',0);
    ModelResultsSeg_ZS(i).coef=mdl.Coefficients;
    ModelResultsSeg_ZS(i).MSE=mdl.MSE;
    ModelResultsSeg_ZS(i).Fitted=mdl.Fitted;
    ModelResultsSeg_ZS(i).rsquared=mdl.Rsquared.Adjusted;
end

rsq_ZS=[ModelResultsSeg_ZS.rsquared];

idx_rsq=find(rsq_ZS>0.3);

for i = 1:length(MatFiles)
    name=strcat(MatFiles(i).name);
    Rs=load(name, 'ROIs');
    Rs=Rs.ROIs;
    F=load(name, 'idx_components');
    F=F.idx_components+1;
    Rs=Rs(:,F);
    MatFiles(i).ROI=Rs;     
end
clearvars GC C S F N name i;

coefficients={};
for idx=1:length(ModelResultsSeg_ZS)
    coef=[ModelResultsSeg_ZS(idx).coef];
    temp=coef.Properties.RowNames;temp=regexp(temp,'x(\d+)','tokens');
    if ~isempty(temp)
        %temp=[temp{:}];temp=[temp{:}];temp=[temp{:}];%temp=str2num(temp);
        for coef_idx=2:height(coef)
            if coef.pValue(coef_idx)<0.05
                coefficients{idx,str2num(temp{coef_idx}{1}{1})}=coef.Estimate(coef_idx);
            end
        end
    end
end
idxempty=cellfun('isempty',coefficients);
coefficients(idxempty)={0};
clearvars idxempty idx coef_idx coef
coefficients=cell2mat(coefficients);

for i=1:length(MatFiles)
    if i==1        
        MatFiles(i).rsquared=rsq_ZS(1:MatFiles(i).GoodNumber);
        MatFiles(i).CoefNb=coefficients(1:MatFiles(i).GoodNumber,:);
    else        
        MatFiles(i).rsquared=rsq_ZS(MatFiles(i-1).GoodNumber+1:MatFiles(i).GoodNumber);
        MatFiles(i).CoefNb=coefficients(MatFiles(i-1).GoodNumber+1:MatFiles(i).GoodNumber,:);
    end
end

colors = [256 0 0; 0 256 0];coloring=jet(20);
for idx=1:length(MatFiles)
    filename=MatFiles(idx).name;    
    %image=MatFiles(idx).Cn;image=image/max(max(image));image=image*128;
    imagename=regexp(filename,'_output_analysis','split');imagename=imagename{1};
    %imagename(ismember(imagename,'_'))=[];
    imagename=strcat('AVG_',imagename,'.tif');
    image=double(imread(imagename));image=image/max(max(image));image=image*128;
    image=uint8(image);
    image2=zeros(size(image(:,:,1)));
    image3=repmat(image,1,1,3);
    %idx_rsquared=find(MatFiles(idx).rsquared>0.05);
    idx_rsquared=find(MatFiles(idx).rsquared>0.3);
    ROIs=MatFiles(idx).ROI;
    ROIs=ROIs(:,idx_rsquared);
    Coefs=MatFiles(idx).CoefNb(idx_rsquared,:);
    for k = 1 : size(ROIs,2)
        image2=zeros(size(image(:,:,1)));
        ROI=full(reshape(ROIs(:,k),size(image,1),size(image,2)));
        image2=image2+ROI;image2=(image2/max(max(image2)));image2=uint8(image2);
        for j=1:3
            image3(:,:,j)=image3(:,:,j)+image2*colors(1,j)*Coefs(k,1)*MatFiles(idx).rsquared(idx_rsquared(k));
            image3(:,:,j)=image3(:,:,j)+image2*colors(2,j)*Coefs(k,2)*MatFiles(idx).rsquared(idx_rsquared(k));
        end
    end
    %image3(:,:,3)=image;
    %name=strcat('_Rsquared',filename(1:end-3),'.tif');
    name=strcat('_RsquaredCoef',imagename(4:end));
    imwrite(image3,name,'tif');
end
%image3=uint8(image3);

