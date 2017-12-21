Stimuli=zeros(6,size(ZS,2));
GCaMP6=[0,1.69644104899772,5.13796058542217,8.27886020152244,10.3756715204800,11.8173714529814,12.2425184714093,10.8571417354877,8.80831829681196,6.91339112244670,5.46959264663869,4.30868766622567,3.42533619066766,2.75378443486879,2.18017250852183,1.72816235135824,1.32732537295463,1.00684435500268,0.730210038304555,0.530242444093118,0.362253250339685,0.227668255288566,0.0869242416152502,0.000718266708050853,-0.0828334873368325]';
idxStart=58;
counter=0;
for i=1:18
    if mod(i,3)==1
        counter=counter+1;
    end
    Stimuli(counter,(idxStart+(i-1)*40):(idxStart+(i-1)*40)+size(GCaMP6,1)-1)=GCaMP6;
end

ModelMultipower=[];
parfor i=1:length(ZS)
    %mdl=stepwiselm(Regressor',DataSet(i,:),'linear','Criterion','adjrsquared','Intercept',false,'Upper','interactions','Verbose',0);
    mdl=fitlm(Stimuli',ZS(i,:),'interactions');
    ModelMultipower(i).coef=mdl.Coefficients;
    %ModelResultsSeg_ZS(i).MSE=mdl.MSE;
    %ModelResultsSeg_ZS(i).Fitted=mdl.Fitted;
    ModelMultipower(i).rsquared=mdl.Rsquared.Adjusted;
end

coefficients={};
for idx=1:length(ModelMultipower)
    coef=[ModelMultipower(idx).coef];
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
clearvars idxempty idx coef_idx coef temp
coefficients=cell2mat(coefficients);

ZS_AVG=zeros(size(ZS,1),246);
ZS_STD=zeros(size(ZS,1),246);
parfor idx_ZS=1:size(ZS,1)
    start=58;
    AVG=[];
    for i=1:3
        AVG(i,:)=ZS(idx_ZS,start:start+40);
        start=start+40;
    end
    STD=std(AVG,1,1);
    AVG=mean(AVG,1);
    AVG=AVG-min(AVG);
    j=1;
    for j=2:6
        for i=1:3
            temp(i,:)=ZS(idx_ZS,start:start+40);
            start=start+40;
        end
        temp_STD=std(temp,1,1);
        temp=mean(temp,1);
        temp=temp-min(temp);
        STD=[STD temp_STD];
        AVG=[AVG temp];
    end
    ZS_AVG(idx_ZS,:)=AVG;
    ZS_STD(idx_ZS,:)=STD;
end

AVG_Stimuli=zeros(6,size(ZS_AVG,2));
GCaMP6=[0,1.69644104899772,5.13796058542217,8.27886020152244,10.3756715204800,11.8173714529814,12.2425184714093,10.8571417354877,8.80831829681196,6.91339112244670,5.46959264663869,4.30868766622567,3.42533619066766,2.75378443486879,2.18017250852183,1.72816235135824,1.32732537295463,1.00684435500268,0.730210038304555,0.530242444093118,0.362253250339685,0.227668255288566,0.0869242416152502,0.000718266708050853,-0.0828334873368325]';
idxStart=1;
for i=1:6    
    AVG_Stimuli(i,(idxStart+(i-1)*40):(idxStart+(i-1)*40)+size(GCaMP6,1)-1)=GCaMP6;
end

ModelMultipower_AVG=[];
parfor i=1:length(ZS_AVG)
    %mdl=stepwiselm(Regressor',DataSet(i,:),'linear','Criterion','adjrsquared','Intercept',false,'Upper','interactions','Verbose',0);
    mdl=fitlm(AVG_Stimuli',ZS_AVG(i,:),'interactions');
    ModelMultipower_AVG(i).coef=mdl.Coefficients;
    %ModelResultsSeg_ZS(i).MSE=mdl.MSE;
    %ModelResultsSeg_ZS(i).Fitted=mdl.Fitted;
    ModelMultipower_AVG(i).rsquared=mdl.Rsquared.Adjusted;
end

options = statset('UseParallel',1); [idxKmeans_ZS_AVG_rsq Cmap_ZS_AVG_rsq]=kmeans(ZS_AVG(idx_rsq_AVG,:),10,'Options',options,'Distance','cityblock','Replicates',3,'MaxIter',1000,'Display','final');
options = statset('UseParallel',1); [idxKmeans_ZS_AVG_rsq2 Cmap_ZS_AVG_rsq2]=kmeans(ZS_AVG(idx_rsq_AVG,:),20,'Options',options,'Distance','cityblock','Replicates',3,'MaxIter',1000,'Display','final');
options = statset('UseParallel',1); [idxKmeans_ZS_AVG_rsq3 Cmap_ZS_AVG_rsq3]=kmeans(ZS_AVG(idx_rsq_AVG,:),5,'Options',options,'Distance','cityblock','Replicates',3,'MaxIter',1000,'Display','final');
[Model_ZS,GoodBetas_ZS]=Test_Regress(Cmap_ZS_AVG_rsq2,AVG_Stimuli,idxKmeans_ZS_AVG_rsq2,0.2);
[Model_ZS2,GoodBetas_ZS2]=Test_Regress(Cmap_ZS_AVG_rsq,AVG_Stimuli,idxKmeans_ZS_AVG_rsq,0.2);
[Model_ZS3,GoodBetas_ZS3]=Test_Regress(Cmap_ZS_AVG_rsq3,AVG_Stimuli,idxKmeans_ZS_AVG_rsq3,0.2);

GoodBetas=GoodBetas_ZS([2 3 7 8 10 12 14 18 20]);

figure;
for i=1:20
    idx_temp=find(idxKmeans_ZS_AVG_rsq2==i);
    plot(mean(ZS_AVG(idx_rsq_AVG(idx_temp),:),1));
    pause
end

idxKmeans_final=zeros(size(ZS,1),1);
idxKmeans_final(idx_rsq_AVG)=idxKmeans_ZS_AVG_rsq2;

Numbers=[1 [MatFiles.GoodNumber]];
counter=1;
idx_Plane=nan(length(ZS),1);
idx_Position=nan(length(ZS),1);
idx_Fish=nan(length(ZS),1);
name=strcat(MatFiles(1).name);
for i=1:length(MatFiles)	
    name=strcat(MatFiles(i).name);
    [Plane,~]=regexp(name,'_(\d+)_','tokens','match');Plane=str2num(Plane{1}{1});
    if strfind(name,'CRO')
        Fish=1;
    elseif strfind(name,'ERO')
        Fish=2;
    elseif strfind(name,'ELO')
        Fish=3;    
    elseif strfind(name,'CLO')
        Fish=4;
    end
    idx_Plane(Numbers(i):Numbers(i+1))=Plane;
    idx_Position(Numbers(i):Numbers(i+1))=Fish;
    [Fish,~]=regexp(name,'Fish2017(\d+)_','tokens','match');Fish=str2num(Fish{1}{1});
    idx_Fish(Numbers(i):Numbers(i+1))=Fish;
end
clearvars i Fish Plane name counter

temp=unique(idx_Fish);
for i=1:length(unique(idx_Fish))
    idx_Fish(idx_Fish==temp(i))=i;
end

counter=1;
x = linspace(0.5,size(ZS,2)/2,size(ZS,2));
rows=length(GoodBetas);
counter=1;
Fighandle=figure;
set(Fighandle, 'Position', [100, 100, 1300, 900]);
for i=GoodBetas
    idx_temp=find(idxKmeans_final==i);
    subplot(rows,5,counter);shadedErrorBar(x, mean(ZS(idx_temp,:),1), std(ZS(idx_temp,:),1,1));axis([0 400 -1 5]);
    subplot(rows,5,counter+1);imagesc(ZS(idx_temp,:),[-1 5]); colormap hot
    subplot(rows,5,counter+2);histogram(idx_Plane(idx_temp));
    subplot(rows,5,counter+3);histogram(idx_Position(idx_temp),[0.5:1:3.5]);h = gca;h.XTick=[1 2 3];h.XTickLabel={'CRO','ERO','ELO'};
    subplot(rows,5,counter+4);histogram(idx_Fish(idx_temp));h = gca;h.XTickLabel={'1','2','3','4','5','6','7'};
    counter=counter+5;
end

counter=1;
x = linspace(0.5,size(ZS_AVG,2)/2,size(ZS_AVG,2));
rows=length(GoodBetas);
counter=1;
Fighandle=figure;
set(Fighandle, 'Position', [100, 100, 1300, 900]);
for i=GoodBetas
    idx_temp=find(idxKmeans_final==i);
    subplot(rows,5,counter);shadedErrorBar(x, mean(ZS_AVG(idx_temp,:),1), mean(ZS_STD(idx_temp,:)));axis([0 123 -0.5 5]);
    subplot(rows,5,counter+1);imagesc(ZS_AVG(idx_temp,:),[-0.5 5]); colormap hot
    subplot(rows,5,counter+2);histogram(idx_Plane(idx_temp));
    subplot(rows,5,counter+3);histogram(idx_Position(idx_temp),[0.5:1:3.5]);h = gca;h.XTick=[1 2 3];h.XTickLabel={'CRO','ERO','ELO'};
    subplot(rows,5,counter+4);histogram(idx_Fish(idx_temp));h = gca;h.XTickLabel={'1','2','3','4','5','6','7'};
    counter=counter+5;
end

colors = distinguishable_colors(length(GoodBetas),[1 1 1; 0 0 0]);
colors = colors*256;
Fighandle=figure;
set(Fighandle, 'Position', [100, 100, 500, 1400]);x = linspace(0.5,size(ZS_AVG,2)/2,size(ZS_AVG,2));
counter=1;counter2=1;xplot=1;yplot=length(GoodBetas);%yplot=ceil(length(GoodBetas)/xplot);
for i=GoodBetas
    idx_temp=find(idxKmeans_final==i);
    subplot(xplot,yplot,counter);plot(x,mean(ZS_AVG(idx_temp,:),1),'color',colors(counter2,:)/256);axis([0 123 -0.5 5]);
    start=30;    
    counter=counter+1;
    counter2=counter2+1;
end

All_ROIs=[];
ROIs_idx=[];
for i = 1:length(MatFiles)
    name=strcat(MatFiles(i).name);
    Rs=load(name, 'ROIs');
    Rs=Rs.ROIs;
    F=load(name, 'idx_components');
    F=F.idx_components+1;
    Rs=Rs(:,F);
    All_ROIs{i}=Rs;
    if i==1
        ROIs_idx(i)=length(F);
    else
        ROIs_idx(i)=ROIs_idx(i-1)+length(F);
    end
end
clearvars GC C S F N name i;

Numbers=[0 [ROIs_idx]];
temp=[];
counter=1;
for i=GoodBetas
    temp{counter}=find(idxKmeans_final==i);
    %tempidx=find(idxKmeans==idx);
    %temp{counter}=GoodClusters_goodmembers(counter).idx;
    counter=counter+1;    
end

for idx=1:length(MatFiles)
    filename=MatFiles(idx).name;
    ROIsNb=[];ClusterNb=[];
    %for k = 1 : length(temp)
    for k = 1 : length(temp)
        tempROIsNb=find([temp{k}]<=Numbers(idx+1));
        if tempROIsNb            
            ROIsNb=[ROIsNb ; temp{k}(tempROIsNb)];
            temp{k}(tempROIsNb)=[];
            ClusterNb=[ClusterNb ; repmat(k,length(tempROIsNb),1)];
        end
    end
    if ROIsNb
        imagename=regexp(filename,'_output_analysis','split');
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
