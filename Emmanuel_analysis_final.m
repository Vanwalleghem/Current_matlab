MatFiles=dir('*analysis_matlab.mat');
name=strcat(MatFiles(1).name);
Calcium=load(name, 'DenoisedTraces');
Calcium=Calcium.DenoisedTraces;
%MatFiles(1).number=size(Calcium,1);
%Spikes=load(name, 'Spikes');
%Spikes=Spikes.Spikes;
%Noise=load(name, 'Noise');
%Noise=Noise.Noise;
%DF=load(name, 'dFonF');
%DF=DF.dFonF;
Fitness=load(name, 'idx_components');
Fitness=Fitness.idx_components+1;
GoodCalcium=Calcium(Fitness,:);
%GoodSpikes=Spikes(Fitness,:);
%GoodNoise=Noise(Fitness,:);
%GoodDF=DF(Fitness,:);
MatFiles(1).GoodNumber=length(Fitness);
for i = 2:length(MatFiles)
name=strcat(MatFiles(i).name);
C=load(name, 'DenoisedTraces');
C=C.DenoisedTraces;
%     if i==3
%         C=[C(:,1) C(:,1) C(:,1:58)];
%     end
% S=load(name, 'Spikes');
% S=S.Spikes;
%N=load(name, 'Noise');
%N=N.Noise;
F=load(name, 'idx_components');
F=F.idx_components+1;
%D=load(name, 'dFonF');
%D=D.dFonF;
GC=C(F,:);
%GS=S(F,:);
%GD=D(F,:);
if size(GC,2)==660
    %Noise=vertcat(Noise,N(:,21:660));
    %GN=N(F,:);
    %Calcium=vertcat(Calcium,C(:,21:660));
    %DF=vertcat(DF,D);
    %Spikes=vertcat(Spikes,S);
    Fitness=horzcat(Fitness,F);
    GoodCalcium=vertcat(GoodCalcium,GC(:,21:660)); %The fish 20+ are longer
    %GoodNoise=vertcat(GoodNoise,GN(:,21:660));
    %GoodDF=vertcat(GoodDF,GD);
    %GoodSpikes=vertcat(GoodSpikes,GS);
else
    %Noise=vertcat(Noise,N);
    %GN=N(F,:);
    %Calcium=vertcat(Calcium,C);
    %DF=vertcat(DF,D);
    %Spikes=vertcat(Spikes,S);
    Fitness=horzcat(Fitness,F);
    GoodCalcium=vertcat(GoodCalcium,GC); %The fish 20+ are longer
    %GoodNoise=vertcat(GoodNoise,GN);
    %GoodDF=vertcat(GoodDF,GD);
    %GoodSpikes=vertcat(GoodSpikes,GS);
end
%MatFiles(i).number=size(Calcium,1);
MatFiles(i).GoodNumber=MatFiles(i-1).GoodNumber+length(F);
end
clearvars GC C S F N name i GS GN N;
ZS=zscore(GoodCalcium,1,2);
ZS=detrend(ZS')';

options = statset('UseParallel',1); [idxKmeans_ZS Cmap_ZS]=kmeans(ZS,50,'Options',options,'Distance','cityblock','Replicates',5,'MaxIter',1000,'Display','final');
[Model_ZS,GoodBetas_ZS]=Test_Regress(Cmap_ZS,Regressors,idxKmeans_ZS,0.1);
GoodBetas=GoodBetas_ZS([11 16 17 29]);

Numbers=[1 [MatFiles.GoodNumber]];
counter=1;
idx_Plane=nan(length(ZS),1);
idx_Fish=nan(length(ZS),1);
name=strcat(MatFiles(1).name);
for i=1:length(MatFiles)	
    name=strcat(MatFiles(i).name);
    [Plane,~]=regexp(name,'Slice(\d+)_','tokens','match');Plane=str2num(Plane{1}{1});    
    idx_Plane(Numbers(i):Numbers(i+1))=Plane;    
    [Fish,~]=regexp(name,'fish(\d+)_','tokens','match');Fish=str2num(Fish{1}{1});
    idx_Fish(Numbers(i):Numbers(i+1))=Fish;
end
clearvars i Fish Plane name counter

test=Numbers;
for i=2:length(test)
    test(i)=test(i)-Numbers(i-1);
end

Stimuli=zeros(10,size(ZS,2));
GCaMP6=[0,1.69644104899772,5.13796058542217,8.27886020152244,10.3756715204800,11.8173714529814,12.2425184714093,10.8571417354877,8.80831829681196,6.91339112244670,5.46959264663869,4.30868766622567,3.42533619066766,2.75378443486879,2.18017250852183,1.72816235135824,1.32732537295463,1.00684435500268,0.730210038304555,0.530242444093118,0.362253250339685,0.227668255288566,0.0869242416152502,0.000718266708050853,-0.0828334873368325]';
idxStart=53;
for i=1:10
    Stimuli(i,(idxStart+(i-1)*60):(idxStart+(i-1)*60)+size(GCaMP6,1)-1)=GCaMP6;
end

ModelResults=[];
parfor i=1:size(ZS,1)
    mdl=stepwiselm(Stimuli',ZS(i,:),'linear','Criterion','adjrsquared','Intercept',true,'Upper','interactions','Verbose',0);
    %mdl=stepwiselm(Stimuli',ZS(i,:),'linear','Criterion','adjrsquared','Upper','interactions','Verbose',0);
    ModelResults(i).coef=mdl.Coefficients;
    ModelResults(i).MSE=mdl.MSE;
    ModelResults(i).Fitted=mdl.Fitted;
    ModelResults(i).rsquared=mdl.Rsquared.Adjusted;
end
rsquare_loom=[ModelResults.rsquared];
idx_rsq=find(rsquare_loom>0.3 & rsquare_loom<1);
figure;
imagesc(ZS(idx_rsq,:), [-0.5 4]);colormap hot

coefficients={};
for idx=1:length(ModelResults)
    coef=[ModelResults(idx).coef];
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

coefficients_bool=coefficients>0;
coefficients_bool=sum(coefficients_bool,2);
rsquare_tone(coefficients_bool<4)=0;

coefficients_Nb=coefficients>0;
coefficients_Nb=sum(coefficients_Nb,2);
idx_coef=find(coefficients_Nb>2);
idx_coef_rsq=intersect(idx_rsq,idx_coef);


options = statset('UseParallel',1); [idxKmeans_ZS_rsq Cmap_ZS_rsq]=kmeans(ZS(idx_coef_rsq,:),10,'Options',options,'Distance','cityblock','Replicates',5,'MaxIter',1000,'Display','final');
[Model_ZS2,GoodBetas_ZS2]=Test_Regress(Cmap_ZS_rsq,Stimuli,idxKmeans_ZS_rsq,0.3);

options = statset('UseParallel',1); [idxKmeans_ZS_rsq2 Cmap_ZS_rsq2]=kmeans(ZS(idx_coef_rsq,:),5,'Options',options,'Distance','cityblock','Replicates',5,'MaxIter',1000,'Display','final');
[Model_ZS3,GoodBetas_ZS3]=Test_Regress(Cmap_ZS_rsq2,Stimuli,idxKmeans_ZS_rsq2,0.3);

GoodBetas=[1 3 4 5 6 10 13 16];
GoodBetas=GoodBetas_ZS2(GoodBetas);
rows=length(GoodBetas);
counter=1;
Fighandle=figure;
set(Fighandle, 'Position', [100, 100, 1300, 900]);
for i=GoodBetas
    idx_temp=find(idxKmeans_ZS_rsq==i);
    subplot(rows,4,counter);plot(mean(ZS(idx_coef_rsq(idx_temp),:),1));
    subplot(rows,4,counter+1);imagesc(ZS(idx_coef_rsq(idx_temp),:),[0 3]);
    subplot(rows,4,counter+2);histogram(idx_Plane(idx_coef_rsq(idx_temp)));    
    subplot(rows,4,counter+3);histogram(idx_Fish(idx_coef_rsq(idx_temp)));
    counter=counter+4;
end

GoodBetas=[2 4 6 7 8 10];
GoodBetas=GoodBetas_ZS3(GoodBetas);
rows=length(GoodBetas);
counter=1;
Fighandle=figure;
set(Fighandle, 'Position', [100, 100, 1300, 900]);
for i=GoodBetas
    idx_temp=find(idxKmeans_ZS_rsq2==i);
    subplot(rows,4,counter);plot(mean(ZS(idx_coef_rsq(idx_temp),:),1));
    subplot(rows,4,counter+1);imagesc(ZS(idx_coef_rsq(idx_temp),:),[0 3]);
    subplot(rows,4,counter+2);histogram(idx_Plane(idx_coef_rsq(idx_temp)));    
    subplot(rows,4,counter+3);histogram(idx_Fish(idx_coef_rsq(idx_temp)));
    counter=counter+4;
end

Test=[];i=1;
name=strcat(MatFiles(i).name);
Rs=load(name, 'ROIs');
Rs=Rs.ROIs;
F=load(name, 'idx_components');
F=F.idx_components+1;
Rs=Rs(:,F);
MatFiles(i).ROI=Rs;
Test(i)=size(Rs,2)==MatFiles(i).GoodNumber;
for i = 3:length(MatFiles)
    name=strcat(MatFiles(i).name);
    Rs=load(name, 'ROIs');
    Rs=Rs.ROIs;
    F=load(name, 'idx_components');
    F=F.idx_components+1;
    Rs=Rs(:,F);
    MatFiles(i).ROI=Rs;
    Test(i)=size(Rs,2)==(MatFiles(i).GoodNumber-MatFiles(i-1).GoodNumber);
end
clearvars GC C S F N name i;

idxKmeans_final=zeros(size(ZS,1),1);
idxKmeans_final(idx_coef_rsq)=idxKmeans_ZS_rsq2;

temp=[];
counter=1;
for i=GoodBetas
    idx_temp=find(idxKmeans_final==i);   
    temp{counter}=idx_temp;    
    counter=counter+1;    
end

Numbers(1)=0;
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
for idx=1:length(MatFiles)
    filename=MatFiles(idx).name;
    ROIsNb=[];ClusterNb=[];
    %for k = 1 : length(temp)
    for k = 1 : length(temp)
        tempROIsNb=find([temp{k}]<=Numbers(idx+1) & [temp{k}]>Numbers(idx));
        if tempROIsNb            
            ROIsNb=[ROIsNb; temp{k}(tempROIsNb)];
            %temp{k}(tempROIsNb)=[];
            ClusterNb=[ClusterNb; repmat(k,length(tempROIsNb),1)];
        end
    end
    if ROIsNb
        imagename=regexp(filename,'_output_analysis','split');
        %imagename=regexp(imagename,'_output_analysis_matlab2.mat','split');
        imagename=strcat(imagename{1},'_mean.tif');
        image=double(imread(imagename));image=image/max(max(image));image=image*64;
        image=uint8(image);
        image2=zeros(size(image(:,:,1)));
        image3=repmat(image,1,1,3);
        ROIs=MatFiles(idx).ROI;       
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
set(Fighandle, 'Position', [100, 100, 1400, 900]);
counter=1;counter2=1;xplot=floor(sqrt(length(GoodBetas)));yplot=ceil(length(GoodBetas)/xplot);
for i=GoodBetas
%     if counter==3
%         counter=counter+1;
%     end
    idx_temp=idx_coef_rsq(find(idxKmeans_ZS_rsq2==i));   
    %idx_temp=idx_rsq(idx_temp);            
    subplot(xplot,yplot,counter);plot(mean(ZS(idx_temp,:),1),'color',colors(counter2,:)/256);%axis([0 131 -1 4]);rectangle('FaceColor','r','Position',[11 -1 10 0.25]);rectangle('FaceColor','r','Position',[51 -1 10 0.25]);rectangle('FaceColor','r','Position',[91 -1 10 0.25]);rectangle('FaceColor','b','Position',[31 -1 10 0.25]);rectangle('FaceColor','b','Position',[71 -1 10 0.25]);rectangle('FaceColor','b','Position',[111 -1 10 0.25]);    
    counter=counter+1;
    counter2=counter2+1;
end

Fish_movements={};
for fish=1:size(Behavioural_responses,2)
    for i=1:4;
        temp=find(Behavioural_responses(:,fish)==i);
        temp=round((temp-112)/3);
        temp(temp<1)=[];
        temp(temp>640)=[];        
        Fish_movements{fish,i}=temp;
    end
    if fish==1
        temp=find(Behavioural_responses(:,fish)==5);
        temp=round((temp-112)/3);
        temp(temp<1)=[];
        temp(temp>640)=[];        
        Fish_movements{fish,1}=temp;
    end
end

Movements=zeros(size(Behavioural_responses,2)-1,4,size(ZS,2));
%GCaMP6=[5.13796058542217,10.3756715204800,11.8173714529814,12.2425184714093,10.8571417354877,8.80831829681196,6.91339112244670,5.46959264663869,4.30868766622567,3.42533619066766,2.75378443486879,2.18017250852183,1.72816235135824,1.32732537295463,1.00684435500268,0.730210038304555,0.530242444093118,0.362253250339685,0.227668255288566,0.0869242416152502]';
GCaMP6=[0,0.5,1.69644104899772,5.13796058542217,8.27886020152244,10.3756715204800,11.8173714529814,12.2425184714093,10.8571417354877,8.80831829681196,6.91339112244670,5.46959264663869,4.30868766622567,3.42533619066766,2.75378443486879,2.18017250852183,1.72816235135824,1.32732537295463,1.00684435500268,0.730210038304555,0.530242444093118,0.362253250339685,0.227668255288566,0.0869242416152502,0.000718266708050853,-0.0828334873368325]';
for fish=1:size(Behavioural_responses,2)-1
    for movement=1:4
        idx=Fish_movements{fish+1,movement};
        if idx
            for i=idx'
                Movements(fish,movement,i:i+size(GCaMP6,1)-1)=GCaMP6;
            end
        end
    end
end
Movements=Movements(:,:,1:size(ZS,2));

Correlation_movement={};
Fish_with_behav=[23 24 25 26 27 28 30 31 32 34 35 36 37 38 39 40 41];counter=1;
for fish=Fish_with_behav
    temp_ZS=ZS(find(idx_Fish==fish),:);Correlation=[];
    parfor idx=1:size(temp_ZS,1)
        for movement=1:4
            temp=corrcoef(squeeze(squeeze(Movements(counter,movement,:))),temp_ZS(idx,:));
            Correlation(movement,idx)=temp(1,2);            
        end
    end
    Correlation_movement{counter}=Correlation;
    counter=counter+1;
end

Threshold=0.2;
for fish=1:length(Fish_with_behav)
    figure;
    temp_ZS=ZS(find(idx_Fish==Fish_with_behav(fish)),:);
    for i=1:4
        Correlation=Correlation_movement{fish};
        subplot(2,2,i);plot(mean(temp_ZS(find(Correlation(i,:)>Threshold),:),1));
        mdl=stepwiselm(Stimuli',mean(temp_ZS(find(Correlation(i,:)>Threshold),:),1),'linear','Criterion','adjrsquared','Intercept',true,'Upper','interactions','Verbose',0);
        title(strcat('lin_reg to looms : ',num2str([mdl.Rsquared.Adjusted]),' nb of ROIs : ', num2str(length(find(Correlation(i,:)>Threshold)))));
    end
end

% loom_temp=find(sum(Stimuli,1));
% for fish=1:length(Fish_with_behav)
%     for movement = 1:4
%         mov_temp=find(Movements(fish,movement,:)>0);
%         Correlation_movements_looms{fish,movement}=intersect(loom_temp,mov_temp);
%     end
% end

Correlation_movements_loomsb={};
for fish=1:length(Fish_with_behav)    
    Mov_temp=Fish_movements{fish+1,1};
    for i=2:4
        Mov_temp=vertcat(Mov_temp,Fish_movements{fish+1,i});
    end
    Correlation_movements_loomsb{fish}=intersect([Fish_movements{1,1}],Mov_temp);
end

Movements_looms=zeros(size(Behavioural_responses,2)-1,size(ZS,2));
GCaMP6=[0,0.5,1.69644104899772,5.13796058542217,8.27886020152244,10.3756715204800,11.8173714529814,12.2425184714093,10.8571417354877,8.80831829681196,6.91339112244670,5.46959264663869,4.30868766622567,3.42533619066766,2.75378443486879,2.18017250852183,1.72816235135824,1.32732537295463,1.00684435500268,0.730210038304555,0.530242444093118,0.362253250339685,0.227668255288566,0.0869242416152502,0.000718266708050853,-0.0828334873368325]';
for fish=1:length(Fish_with_behav)
    if Correlation_movements_loomsb{fish}
        for i=Correlation_movements_loomsb{fish}'
            Movements_looms(fish,i:i+size(GCaMP6,1)-1)=GCaMP6;
        end
    end
end
Movements_looms=Movements_looms(:,1:size(ZS,2));

Correlation_movement_loom={};
Fish_with_behav=[23 24 25 26 27 28 30 31 32 34 35 36 37 38 39 40 41];counter=1;
for fish=Fish_with_behav
    temp_ZS=ZS(find(idx_Fish==fish),:);Correlation=[];
    parfor idx=1:size(temp_ZS,1)
        temp=corrcoef(squeeze(squeeze(Movements_looms(counter,:))),temp_ZS(idx,:));
        Correlation(idx)=temp(1,2);        
    end
    Correlation_movement_loom{counter}=Correlation;
    counter=counter+1;
end



figure;Threshold=0.2;
counter=1;counter2=1;xplot=floor(sqrt(length(Fish_with_behav)));yplot=ceil(length(Fish_with_behav)/xplot);
for fish=1:length(Fish_with_behav)    
    temp_ZS=ZS(find(idx_Fish==Fish_with_behav(fish)),:);    
    Corr_temp=Correlation_movement_loom{fish};
    subplot(xplot,yplot,fish);plot(mean(temp_ZS(find(Corr_temp>Threshold),:),1));
    for i=1:10
        rectangle('FaceColor','r','Position',[40+((i-1)*60) -1 10 0.25]);
    end
    %mdl=stepwiselm(Stimuli',mean(temp_ZS(find(Corr_temp>Threshold),:),1),'linear','Criterion','adjrsquared','Intercept',true,'Upper','interactions','Verbose',0);
    %title(num2str([mdl.Rsquared.Adjusted]));
    title(strcat(' nb of ROIs : ', num2str(length(find(Corr_temp>Threshold)))));
end

figure;
xplot=floor(sqrt(length(Fish_with_behav)));yplot=ceil(length(Fish_with_behav)/xplot);
for fish=1:length(Fish_with_behav)    
    temp_ZS=ZS(find(idx_Fish==Fish_with_behav(fish)),:);    
    Corr_temp=Correlation_movement_loom{fish};
    subplot(xplot,yplot,fish);imagesc(temp_ZS(find(Corr_temp>Threshold),:),[-0.5 3]);    
end

for i=1:4
    figure;
    xplot=floor(sqrt(length(Fish_with_behav)));yplot=ceil(length(Fish_with_behav)/xplot);
    for fish=1:length(Fish_with_behav)
        temp_ZS=ZS(find(idx_Fish==Fish_with_behav(fish)),:);        
        Correlation=Correlation_movement{fish};
        subplot(xplot,yplot,fish);plot(mean(temp_ZS(find(Correlation(i,:)>Threshold),:),1));
        %mdl=stepwiselm(Stimuli',mean(temp_ZS(find(Correlation(i,:)>Threshold),:),1),'linear','Criterion','adjrsquared','Intercept',true,'Upper','interactions','Verbose',0);
        %title(num2str([mdl.Rsquared.Adjusted]));
        title(strcat(' nb of ROIs : ', num2str(length(find(Correlation(i,:))))));
    end
end

for i=1:4
    figure;
    xplot=floor(sqrt(length(Fish_with_behav)));yplot=ceil(length(Fish_with_behav)/xplot);
    for fish=1:length(Fish_with_behav)
        temp_ZS=ZS(find(idx_Fish==Fish_with_behav(fish)),:);        
        Correlation=Correlation_movement{fish};
        subplot(xplot,yplot,fish);imagesc(temp_ZS(find(Correlation(i,:)>0.25),:),[-0.5 3]);       
    end
end


figure;Threshold=0.2;
xplot=floor(sqrt(length(Fish_with_behav)));yplot=ceil(length(Fish_with_behav)/xplot);
for fish=1:length(Fish_with_behav)    
    idx_fish=find(idx_Fish==Fish_with_behav(fish));
    temp_ZS=ZS(intersect(idx_coef_rsq,idx_fish),:);    
    Corr_temp=Correlation_movement_loom{fish};
    [~,idx_interecpt]=intersect(idx_fish,idx_coef_rsq);    
    subplot(xplot,yplot,fish);imagesc(temp_ZS(find(Corr_temp(idx_interecpt)>Threshold),:),[-0.5 3]);
    for j=1:10
        rectangle('FaceColor','r','Position',[40+((i-1)*60) -1 10 0.25]);
    end
end


for i=1:4
    figure;
    xplot=floor(sqrt(length(Fish_with_behav)));yplot=ceil(length(Fish_with_behav)/xplot);
    for fish=1:length(Fish_with_behav)
        idx_fish=find(idx_Fish==Fish_with_behav(fish));
        temp_ZS=ZS(find(idx_Fish==Fish_with_behav(fish)),:);        
        Correlation=Correlation_movement{fish};
        [~,idx_interecpt]=intersect(idx_fish,idx_coef_rsq);    
        subplot(xplot,yplot,fish);imagesc(temp_ZS(find(Correlation(i,idx_interecpt)>Threshold),:),[-0.5 3]);       
        for j=1:10
            rectangle('FaceColor','r','Position',[40+((j-1)*60) -1 10 0.25]);
        end
    end
end

