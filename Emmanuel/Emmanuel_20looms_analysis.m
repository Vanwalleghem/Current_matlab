MatFiles=dir('*0lo*_step10*output_analysis_matlab.mat');
name=strcat(MatFiles(1).name);
Calcium=load(name, 'DenoisedTraces');
Calcium=Calcium.DenoisedTraces;
MatFiles(1).number=size(Calcium,1);
Spikes=load(name, 'Spikes');
Spikes=Spikes.Spikes;
Noise=load(name, 'Noise');
Noise=Noise.Noise;
%DF=load(name, 'dFonF');
%DF=DF.dFonF;
Fitness=load(name, 'idx_components');
Fitness=Fitness.idx_components+1;
GoodCalcium=Calcium(Fitness,:);
GoodSpikes=Spikes(Fitness,:);
GoodNoise=Noise(Fitness,:);
%GoodDF=DF(Fitness,:);
MatFiles(1).GoodNumber=length(Fitness);
for i = 2:length(MatFiles)
name=strcat(MatFiles(i).name);
C=load(name, 'DenoisedTraces');
C=C.DenoisedTraces;
%     if i==3
%         C=[C(:,1) C(:,1) C(:,1:58)];
%     end
S=load(name, 'Spikes');
S=S.Spikes;
N=load(name, 'Noise');
N=N.Noise;
F=load(name, 'idx_components');
F=F.idx_components+1;
%D=load(name, 'dFonF');
%D=D.dFonF;
GC=C(F,:);
GS=S(F,:);
%GD=D(F,:);
Noise=vertcat(Noise,N);
GN=N(F,:);
Calcium=vertcat(Calcium,C);
%DF=vertcat(DF,D);
Spikes=vertcat(Spikes,S);
Fitness=horzcat(Fitness,F);
GoodCalcium=vertcat(GoodCalcium,GC);
GoodNoise=vertcat(GoodNoise,GN);
%GoodDF=vertcat(GoodDF,GD);
GoodSpikes=vertcat(GoodSpikes,GS);
MatFiles(i).number=size(Calcium,1);
MatFiles(i).GoodNumber=MatFiles(i-1).GoodNumber+length(F);
end
clearvars GC C S F N name i GS GN N;
%ZS=zscore(detrend(GoodCalcium(:,25:end)')',1,2);
ZS=zscore(GoodCalcium,1,2);
options = statset('UseParallel',1); [idxKmeans_ZS Cmap_ZS]=kmeans(ZS,10,'Options',options,'Distance','cityblock','Replicates',5,'MaxIter',1000,'Display','final');

Stimuli=zeros(1,size(GoodCalcium,2));
GCaMP6=[-0.104392135015146,1.69644104899772,5.13796058542217,8.27886020152244,10.3756715204800,11.8173714529814,12.2425184714093,10.8571417354877,8.80831829681196,6.91339112244670,5.46959264663869,4.30868766622567,3.42533619066766,2.75378443486879,2.18017250852183,1.72816235135824,1.32732537295463,1.00684435500268,0.730210038304555,0.530242444093118,0.362253250339685,0.227668255288566,0.0869242416152502,0.000718266708050853,-0.0828334873368325]';
GCaMP6=interp(GCaMP6,2);
idx=105;
for i=1:20
    Stimuli(1,idx:idx+size(GCaMP6,1)-1)=GCaMP6';
    idx=idx+200;
end
%Stimuli=Stimuli(:,25:end);

Stimuli2=zeros(20,size(GoodCalcium,2));
Stimuli2=zeros(20,4000);
GCaMP6=[-0.104392135015146,1.69644104899772,5.13796058542217,8.27886020152244,10.3756715204800,11.8173714529814,12.2425184714093,10.8571417354877,8.80831829681196,6.91339112244670,5.46959264663869,4.30868766622567,3.42533619066766,2.75378443486879,2.18017250852183,1.72816235135824,1.32732537295463,1.00684435500268,0.730210038304555,0.530242444093118,0.362253250339685,0.227668255288566,0.0869242416152502,0.000718266708050853,-0.0828334873368325]';
GCaMP6=interp(GCaMP6,2);
idx=105;
for i=1:20
    Stimuli2(i,idx:idx+size(GCaMP6,1)-1)=GCaMP6';
    idx=idx+200;
end
Stimuli2=Stimuli2(:,25:end);

options = statset('UseParallel',1); [idxKmeans_ZS Cmap_ZS]=kmeans(ZS,50,'Options',options,'Distance','cityblock','Replicates',5,'MaxIter',1000,'Display','final');
[Model_ZS,GoodBetas_ZS]=Test_Regress(Cmap_ZS,Stimuli,idxKmeans_ZS,0.1);


ZS=zscore(GoodCalcium,1,2);
ModelResults=[];
parfor i=1:size(ZS,1)
    mdl=stepwiselm(Stimuli',ZS(i,:),'linear','Criterion','adjrsquared','Intercept',false,'Upper','interactions','Verbose',0);
    %mdl=stepwiselm(Stimuli',ZS(i,:),'linear','Criterion','adjrsquared','Upper','interactions','Verbose',0);
    ModelResults(i).coef=mdl.Coefficients;
    ModelResults(i).MSE=mdl.MSE;
    ModelResults(i).Fitted=mdl.Fitted;
    ModelResults(i).rsquared=mdl.Rsquared.Adjusted;
end
rsquare_loom=[ModelResults.rsquared];
idx_rsq=find(rsquare_loom>0.025);
figure;
imagesc(ZS(idx_rsq,:));


ModelResults2=[];
ZS2=zscore(GoodCalcium(:,25:end),1,2);
parfor i=1:size(ZS,1)
    %mdl=stepwiselm(Stimuli',ZS(i,:),'linear','Criterion','adjrsquared','Intercept',false,'Upper','interactions','Verbose',0);
    mdl=stepwiselm(Stimuli2',detrend(ZS2(i,:))','linear','Criterion','adjrsquared','Upper','interactions','Verbose',0);
    ModelResults2(i).coef=mdl.Coefficients;
    ModelResults2(i).MSE=mdl.MSE;
    ModelResults2(i).Fitted=mdl.Fitted;
    ModelResults2(i).rsquared=mdl.Rsquared.Adjusted;
end
rsquare_loom2=[ModelResults2.rsquared];
figure;histogram(rsquare_loom2);

idx_rsq2=find(rsquare_loom2>0.15);

figure;imagesc(ZS2(idx_rsq2,:));

Numbers=[1 [MatFiles.GoodNumber]];
counter=1;
idx_Fish=nan(length(GoodCalcium),1);
name=strcat(MatFiles(1).name);
for i=1:length(MatFiles)	
    name=strcat(MatFiles(i).name);
    if ~isempty(regexp(name,'fish#(\d+)_','tokens'))
        [Fish,~]=regexp(name,'fish#(\d+)_','tokens','match');Fish=str2num(Fish{1}{1});        
    elseif ~isempty(regexp(name,'fish#(\d+)6','tokens'))
        [Fish,~]=regexp(name,'fish#(\d+)6','tokens','match');Fish=str2num(Fish{1}{1});
    else
        [Fish,~]=regexp(name,'fish#(\d+)-','tokens','match');Fish=str2num(Fish{1}{1});
    end
    idx_Fish(Numbers(i):Numbers(i+1))=Fish;
end
clearvars i Fish Plane name counter

% figure;histogram(idx_Fish(idx_rsq2));
% 
% Fishes=[133 134 136 137 138 139 140 141];
% counter=1;xplot=floor(sqrt(length(Fishes)));yplot=ceil(length(Fishes)/xplot);
% Fighandle=figure;
% set(Fighandle, 'Position', [100, 100, 1200, 900]);
% 
% for i=[133 134 136 137 138 139 140 141]
%     idx_fish=find(idx_Fish(idx_rsq)==i);
%     subplot(xplot,yplot,counter);plot(mean(ZS(idx_rsq(idx_fish),:),1));
%     counter=counter+1;
% end
% 
for i=1:length(MatFiles)
    if i==1        
        MatFiles(i).rsquared=rsquare_loom(1:MatFiles(i).GoodNumber);
    else        
        MatFiles(i).rsquared=rsquare_loom(MatFiles(i-1).GoodNumber+1:MatFiles(i).GoodNumber);
    end
end


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

% colors = [256 0 0];
% for idx=1:length(MatFiles)
%     filename=MatFiles(idx).name;    
%     %image=MatFiles(idx).Cn;image=image/max(max(image));image=image*128;
%     imagename=regexp(filename,'_output_analysis','split');imagename=imagename{1};
%     imagename(ismember(imagename,'_'))=[];
%     imagename=strcat('AVG_',imagename,'.tif');
%     image=double(imread(imagename));image=image/max(max(image));image=image*128;
%     image=uint8(image);
%     image2=zeros(size(image(:,:,1)));
%     image3=repmat(image,1,1,3);
%     %idx_rsquared=find(MatFiles(idx).rsquared>0.05);
%     idx_rsquared=find(MatFiles(idx).rsquared>0.01);
%     ROIs=MatFiles(idx).ROI;
%     ROIs=ROIs(:,idx_rsquared);
%     for k = 1 : size(ROIs,2)
%         image2=zeros(size(image(:,:,1)));
%         ROI=full(reshape(ROIs(:,k),size(image,1),size(image,2)));
%         image2=image2+ROI;image2=(image2/max(max(image2)));image2=uint8(image2);
%         for j=1:3
%             image3(:,:,j)=image3(:,:,j)+image2*colors(j)*MatFiles(idx).rsquared(idx_rsquared(k))*5;
%         end
%     end
%     %image3(:,:,3)=image;
%     %name=strcat('_Rsquared',filename(1:end-3),'.tif');
%     name=strcat('_Rsquared',imagename(4:end));
%     imwrite(image3,name,'tif');
% end
% %image3=uint8(image3);
% 
% Fishes=[133 134 136 137 138 139 140 141];
% counter=1;xplot=floor(sqrt(length(Fishes)));yplot=ceil(length(Fishes)/xplot);
% Fighandle=figure;
% set(Fighandle, 'Position', [100, 100, 1200, 900]);
% 
% for i=[133 134 136 137 138 139 140 141]
%     idx_fish=find(idx_Fish(idx_rsq2)==i);
%     subplot(xplot,yplot,counter);plot(mean(ZS(idx_rsq2(idx_fish),:),1));
%     counter=counter+1;
% end
% 
% rsquare_loom2_back=rsquare_loom2;


for i=1:length(MatFiles)
    if i==1        
        MatFiles(i).rsquared2=rsquare_loom2(1:MatFiles(i).GoodNumber);
    else        
        MatFiles(i).rsquared2=rsquare_loom2(MatFiles(i-1).GoodNumber+1:MatFiles(i).GoodNumber);
    end
end

% colors = [256 0 0];
% for idx=1:length(MatFiles)
%     filename=MatFiles(idx).name;    
%     %image=MatFiles(idx).Cn;image=image/max(max(image));image=image*128;
%     imagename=regexp(filename,'_output_analysis','split');imagename=imagename{1};
%     imagename(ismember(imagename,'_'))=[];
%     imagename=strcat('AVG_',imagename,'.tif');
%     image=double(imread(imagename));image=image/max(max(image));image=image*128;
%     image=uint8(image);
%     image2=zeros(size(image(:,:,1)));
%     image3=repmat(image,1,1,3);
%     %idx_rsquared=find(MatFiles(idx).rsquared>0.05);
%     idx_rsquared2=find(MatFiles(idx).rsquared2>0.2);
%     ROIs=MatFiles(idx).ROI;
%     ROIs=ROIs(:,idx_rsquared2);
%     for k = 1 : size(ROIs,2)
%         image2=zeros(size(image(:,:,1)));
%         ROI=full(reshape(ROIs(:,k),size(image,1),size(image,2)));
%         image2=image2+ROI;image2=(image2/max(max(image2)));image2=uint8(image2);
%         for j=1:3
%             image3(:,:,j)=image3(:,:,j)+image2*colors(j)*MatFiles(idx).rsquared2(idx_rsquared2(k))/2;
%         end
%     end
%     %image3(:,:,3)=image;
%     %name=strcat('_Rsquared',filename(1:end-3),'.tif');
%     name=strcat('_Rsquared2',imagename(4:end));
%     imwrite(image3,name,'tif');
% end
% %image3=uint8(image3);

coefficients={};
for idx=1:length(ModelResults2)
    coef=[ModelResults2(idx).coef];
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

%[coef_sorted idx_sort]=sortrows(coefficients,[1 2 3 4 5]);

coefficients_bool=coefficients(:,1:6)>0;
coefficients_bool=sum(coefficients_bool,2);
rsquare_loom2(coefficients_bool<4)=0;

coefficients_Nb=coefficients>0;
coefficients_Nb=sum(coefficients_Nb,2);

for i=1:length(MatFiles)
    if i==1        
        MatFiles(i).rsquared2=rsquare_loom2(1:MatFiles(i).GoodNumber);
        MatFiles(i).CoefNb=coefficients_Nb(1:MatFiles(i).GoodNumber);
    else        
        MatFiles(i).rsquared2=rsquare_loom2(MatFiles(i-1).GoodNumber+1:MatFiles(i).GoodNumber);
        MatFiles(i).CoefNb=coefficients_Nb(MatFiles(i-1).GoodNumber+1:MatFiles(i).GoodNumber);
    end
end

colors = [256 0 0];coloring=jet(20);
for idx=1:length(MatFiles)
    filename=MatFiles(idx).name;    
    %image=MatFiles(idx).Cn;image=image/max(max(image));image=image*128;
    imagename=regexp(filename,'_output_analysis','split');imagename=imagename{1};
    imagename(ismember(imagename,'_'))=[];
    imagename=strcat('AVG_',imagename,'.tif');
    image=double(imread(imagename));image=image/max(max(image));image=image*128;
    image=uint8(image);
    image2=zeros(size(image(:,:,1)));
    image3=repmat(image,1,1,3);
    %idx_rsquared=find(MatFiles(idx).rsquared>0.05);
    idx_rsquared2=find(MatFiles(idx).rsquared2>0.2);
    ROIs=MatFiles(idx).ROI;
    ROIs=ROIs(:,idx_rsquared2);
    Coefs=MatFiles(idx).CoefNb(idx_rsquared2');
    for k = 1 : size(ROIs,2)
        image2=zeros(size(image(:,:,1)));
        ROI=full(reshape(ROIs(:,k),size(image,1),size(image,2)));
        image2=image2+ROI;image2=(image2/max(max(image2)));image2=uint8(image2);
        for j=1:3
            image3(:,:,j)=image3(:,:,j)+image2*coloring(Coefs(k),j)*256;
        end
    end
    %image3(:,:,3)=image;
    %name=strcat('_Rsquared',filename(1:end-3),'.tif');
    name=strcat('_RsquaredCoef',imagename(4:end));
    imwrite(image3,name,'tif');
end
%image3=uint8(image3);


for idx=1:length(MatFiles)
    filename=MatFiles(idx).name;    
    %image=MatFiles(idx).Cn;image=image/max(max(image));image=image*128;
    imagename=regexp(filename,'_output_analysis','split');imagename=imagename{1};
    imagename(ismember(imagename,'_'))=[];
    imagename=strcat('AVG_',imagename,'.tif');
    image=double(imread(imagename));image=image/max(max(image));image=image*128;
    image=uint8(image);
    image2=zeros(size(image(:,:,1)));
    image3=repmat(image,1,1,3);
    %idx_rsquared=find(MatFiles(idx).rsquared>0.05);
    idx_rsquared=find(MatFiles(idx).rsquared>0.025);
    ROIs=MatFiles(idx).ROI;
    ROIs=ROIs(:,idx_rsquared);
    for k = 1 : size(ROIs,2)
        image2=zeros(size(image(:,:,1)));
        ROI=full(reshape(ROIs(:,k),size(image,1),size(image,2)));
        image2=image2+ROI;image2=(image2/max(max(image2)));image2=uint8(image2);
        for j=1:3
            image3(:,:,j)=image3(:,:,j)+image2*colors(j)*MatFiles(idx).rsquared(idx_rsquared(k))*5;
        end
    end
    %image3(:,:,3)=image;
    %name=strcat('_Rsquared',filename(1:end-3),'.tif');
    name=strcat('_Rsquared',imagename(4:end));
    imwrite(image3,name,'tif');
end
%image3=uint8(image3);