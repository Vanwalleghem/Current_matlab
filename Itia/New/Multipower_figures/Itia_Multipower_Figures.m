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
i=1;colors{i}=[0.14 1 0.14; 0.7 0.4 1];regionName=ItiaList{i};LinReg.(regionName).GoodBetas_AVG_final=LinReg.(regionName).GoodBetas_AVG([1 2]);
i=2;colors{i}=[0.14 1 0.14; 0.7 0.4 1];regionName=ItiaList{i};LinReg.(regionName).GoodBetas_AVG_final=LinReg.(regionName).GoodBetas_AVG([1 2]);
i=3;colors{i}=[0.14 1 0.14; 0.7 0.4 1];regionName=ItiaList{i};LinReg.(regionName).GoodBetas_AVG_final=LinReg.(regionName).GoodBetas_AVG([3 1]);
i=4;colors{i}=[0.14 1 0.14];regionName=ItiaList{i};LinReg.(regionName).GoodBetas_AVG_final=LinReg.(regionName).GoodBetas_AVG([1]);
i=5;colors{i}=[0.7 0.4 1;0 0.6 0.6];regionName=ItiaList{i};LinReg.(regionName).GoodBetas_AVG_final=LinReg.(regionName).GoodBetas_AVG([1 2]);
%i=6;colors{i}=[0.14 1 0.14; 0.7 0.4 1; 0 0.6 0.6];regionName=ItiaList{i};LinReg.(regionName).GoodBetas_AVG_final=LinReg.(regionName).GoodBetas_AVG([1 5 4]);
i=6;colors{i}=[0.14 1 0.14; 0.7 0.4 1; 0 0.6 0.6];regionName=ItiaList{i};LinReg.(regionName).GoodBetas_AVG_final=LinReg.(regionName).GoodBetas_AVG([4 5]);
i=7;colors{i}=[0.14 1 0.14; 0.7 0.4 1];regionName=ItiaList{i};LinReg.(regionName).GoodBetas_AVG_final=LinReg.(regionName).GoodBetas_AVG([1 2]);
i=8;colors{i}=[0.14 1 0.14];regionName=ItiaList{i};LinReg.(regionName).GoodBetas_AVG_final=LinReg.(regionName).GoodBetas_AVG([1]);
i=9;colors{i}=[0.14 1 0.14; 0.7 0.4 1;0 0.6 0.6];regionName=ItiaList{i};LinReg.(regionName).GoodBetas_AVG_final=LinReg.(regionName).GoodBetas_AVG([3 2 1]);
i=10;colors{i}=[0.14 1 0.14; 0.7 0.4 1];regionName=ItiaList{i};LinReg.(regionName).GoodBetas_AVG_final=LinReg.(regionName).GoodBetas_AVG([4 3]);

figure;
counter=1;xplot=length(ItiaList);yplot=3;counter2=1;
x = linspace(0.25,246/4,246);
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
        subplot(xplot,yplot,counter+(j-1));        
        for k=0:5
            meanToPlot=mean(LinReg.(regionName).ZS_AVG(idx_temp2,5+(k*40):40+(k*40)),1);
            %conf_95=confidence_intervals(LinReg.(regionName).ZS_AVG(idx_temp2,5+(k*40):40+(k*40)),95);
            std_95=std(LinReg.(regionName).ZS_AVG(idx_temp2,5+(k*40):40+(k*40)),1,1);
            hold on;%plot(x(5+(k*40):40+(k*40)),meanToPlot-mean(meanToPlot(1:4)),'color',colors{i}(j,:));title(num2str(length(idx_temp2)));
            H=shadedErrorBar(x(5+(k*40):40+(k*40)),meanToPlot-mean(meanToPlot(1:4)), std_95);axis([0 60 -3 3]);
            H.mainLine.Color=colors{i}(j,:);
            H.patch.FaceColor=colors{i}(j,:)/2;
            H.edge(1).Color=colors{i}(j,:)/2;
            H.edge(2).Color=colors{i}(j,:)/2;
        end
        hold off;
    end
    counter2=counter2+yplot;
end
clearvars i j k counter counter2 H idx_temp GoodBet_temp std_95 meanToPlot

Index_ELO=strfind(Fish_list, 'ELO');
ELO_fish=find(not(cellfun('isempty', Index_ELO)));
ERO_fish=find(cellfun('isempty', Index_ELO));


fname = 'Huc_H2B_RFP_8bit.tif';
info = imfinfo(fname);
num_images = numel(info);
Template=zeros(info(1).Height,info(1).Width,length(info),'uint8');
for k = 1:num_images
    image_temp = imread(fname, k, 'Info', info)';image_temp=double(image_temp);
    image_temp=image_temp/max(max(image_temp));image_temp=image_temp*128;
    Template(:,:,k) = image_temp';    
end
clearvars info fname num_images k i j 

imageSizeY = size(Template,1);
imageSizeX = size(Template,2);
[columnsInImage rowsInImage] = meshgrid(1:imageSizeX, 1:imageSizeY);radius =3;
for i=1:length(ItiaList)
    Template_ELO=zeros([size(Template),3]);
    Template_ERO=zeros([size(Template),3]);
    regionName=ItiaList{i};
    idx_temp=find(LinReg.(regionName).rsquared>0.1);
    ROIs_temp=ROIsPerBrain.(regionName).ROIs(idx_temp,:);
    Numbers=ROIsPerBrain.(regionName).Numbers;Numbers=[0 Numbers];
    for j=1:length(LinReg.(regionName).GoodBetas_AVG_final)
        idx_goodROIs=find(LinReg.(regionName).KmeansIdx_select_AVG==LinReg.(regionName).GoodBetas_AVG_final(j));
        idx_goodROIs_correct=idx_temp(idx_goodROIs);
        idx_ELO=[];
        idx_ERO=[];
        for fishNb=1:length(Fish_list)
            if Index_ELO{fishNb}                
                idx_ELOtemp=find(Numbers(fishNb) < idx_goodROIs_correct & idx_goodROIs_correct <= Numbers(fishNb+1));
                if isempty(idx_ELO)
                    idx_ELO=idx_ELOtemp;
                else
                    idx_ELO=horzcat(idx_ELO, idx_ELOtemp);
                end
            else
                idx_EROtemp=find(Numbers(fishNb) < idx_goodROIs_correct & idx_goodROIs_correct <= Numbers(fishNb+1));
                if isempty(idx_ERO)
                    idx_ERO=idx_EROtemp;
                else
                    idx_ERO=horzcat(idx_ERO, idx_EROtemp);
                end
            end
        end        
        goodROIs=ROIs_temp(idx_goodROIs,:);        
        for slice=1:size(Template_ELO,3)
            idx_slice=find(goodROIs(:,3)==slice);
            if idx_slice
                for roi_nb=1:length(idx_slice)
                    xcoord=goodROIs(idx_slice(roi_nb),1);
                    ycoord=goodROIs(idx_slice(roi_nb),2);
                    if xcoord>1 & ycoord>1
                        circlePixels = (rowsInImage - ycoord).^2 + (columnsInImage - xcoord).^2 <= radius.^2;
                        for col=1:3                           
                            if ismember(idx_slice(roi_nb),idx_ELO)
                                image_temp=squeeze(squeeze(Template_ELO(:,:,slice,col)));
                                image_temp(circlePixels)=colors{i}(j,col)*256;
                                Template_ELO(:,:,slice,col)=image_temp;
                            elseif ismember(idx_slice(roi_nb),idx_ERO)
                                image_temp=squeeze(squeeze(Template_ERO(:,:,slice,col)));
                                image_temp(circlePixels)=colors{i}(j,col)*256;
                                Template_ERO(:,:,slice,col)=image_temp;
                            else
                                error_roi=[idx_slice roi_nb];
                                break;
                            end
                        end
                    end
                end
            end
        end        
    end
    %InBrain_select_ELO{i}=uint8(Template_ELO);
    InBrain_select_ERO{i}=uint8(Template_ERO);
end
clearvars goodROIs ROIs_temp regionName idx_temp Template2 idx_slice roi_nb xcoord ycoord circlePixels image_temp columnsInImage rowsInImage col corr_temp Fighdandle framerate GoodBet_temp H idx_g Model_cereb

for i=1:length(ItiaList)
    regionName=ItiaList{i};
    OutputName=strcat('Zbrain_Multipower_ELO_HighThresh',regionName,'.tif');    
    delete(OutputName);
    Template2=InBrain_select_ELO{i};
    for slice=1:size(Template2,3)
        image_temp=uint8(squeeze(Template2(:, :,slice,:)));
        imwrite(image_temp, OutputName, 'WriteMode', 'append');
    end
end
clearvars OutputName radius slice Template Threshold xplot ZS_temp
clearvars goodROIs ROIs_temp regionName idx_temp Template2 idx_slice roi_nb xcoord ycoord circlePixels image_temp columnsInImage rowsInImage col corr_temp Fighdandle framerate GoodBet_temp H idx_g Model_cereb

for i=1:length(ItiaList)
    regionName=ItiaList{i};
    OutputName=strcat('Zbrain_Multipower_ERO_HighThresh',regionName,'.tif');    
    delete(OutputName);
    Template2=InBrain_select_ERO{i};
    for slice=1:size(Template2,3)
        image_temp=uint8(squeeze(Template2(:, :,slice,:)));
        imwrite(image_temp, OutputName, 'WriteMode', 'append');
    end
end
clearvars OutputName radius slice Template Threshold xplot ZS_temp
clearvars goodROIs ROIs_temp regionName idx_temp Template2 idx_slice roi_nb xcoord ycoord circlePixels image_temp columnsInImage rowsInImage col corr_temp Fighdandle framerate GoodBet_temp H idx_g Model_cereb
clearvars idx_ELO idx_ELOtemp idx_ERO idx_EROtemp idx_goodROIs idx_goodROIs_correct Template_ELO Template_ERO
clearvars InBrain_select_ELO InBrain_select_ERO

for i=1:length(ItiaList)
    regionName=ItiaList{i};
    idx_temp=find(LinReg.(regionName).rsquared>0.1);
    ROIs_temp=ROIsPerBrain.(regionName).ROIs(idx_temp,:);
    Numbers=ROIsPerBrain.(regionName).Numbers;Numbers=[0 Numbers];
    for j=1:length(LinReg.(regionName).GoodBetas_AVG_final)
        idx_goodROIs=find(LinReg.(regionName).KmeansIdx_select_AVG==LinReg.(regionName).GoodBetas_AVG_final(j));
        idx_goodROIs_correct=idx_temp(idx_goodROIs);
        idx_ELO=[];
        idx_ERO=[];
        for fishNb=1:length(Fish_list)
            if Index_ELO{fishNb}                
                idx_ELOtemp=find(Numbers(fishNb) < idx_goodROIs_correct & idx_goodROIs_correct <= Numbers(fishNb+1));
                if isempty(idx_ELO)
                    idx_ELO=idx_ELOtemp;
                else
                    idx_ELO=horzcat(idx_ELO, idx_ELOtemp);
                end
            else
                idx_EROtemp=find(Numbers(fishNb) < idx_goodROIs_correct & idx_goodROIs_correct <= Numbers(fishNb+1));
                if isempty(idx_ERO)
                    idx_ERO=idx_EROtemp;
                else
                    idx_ERO=horzcat(idx_ERO, idx_EROtemp);
                end
            end
        end
        goodROIs=ROIs_temp(idx_goodROIs,:); 
        ROI_coord_ELO{i,j}=goodROIs(idx_ELO,:);
        ROI_coord_ERO{i,j}=goodROIs(idx_ERO,:);    
    end
end
clearvars goodROIs ROIs_temp regionName idx_temp Template2 idx_slice roi_nb xcoord ycoord circlePixels image_temp columnsInImage rowsInImage col corr_temp Fighdandle framerate GoodBet_temp H idx_g Model_cereb

for i=1:length(ItiaList)
    temp=[];
    for j=1:3
        if ~isempty(ROI_coord_ERO{i,j})
            if isempty(temp)
                temp=ROI_coord_ERO{i,j};
            else
                temp=vertcat(temp,ROI_coord_ERO{i,j});
            end 
        end
    end        
	[pca_coord{i},score_coord{i},~,~,explained_coord{i},~] = pca(temp);
end
    

for i=1:length(ItiaList)
    regionName=ItiaList{i};
    for j=1:3
        if length(LinReg.(regionName).GoodBetas_AVG_final)==2
            [~,p_val_ks{i,j},~] = kstest2(score_coord{i}(1:length(ROI_coord_ERO{i,1}),j),score_coord{i}(length(ROI_coord_ERO{i,1}+1):end,j));
        elseif length(LinReg.(regionName).GoodBetas_AVG_final)==3
            [~,p_val_temp,~] = kstest2(score_coord{i}(1:length(ROI_coord_ERO{i,1}),j),score_coord{i}(length(ROI_coord_ERO{i,1})+1:length(ROI_coord_ERO{i,1})+1+length(ROI_coord_ERO{i,2}),j));
            [~,p_val_temp2,~] = kstest2(score_coord{i}(1:length(ROI_coord_ERO{i,1}),j),score_coord{i}(length(ROI_coord_ERO{i,2})+1:end,j));
            [~,p_val_temp3,~] = kstest2(score_coord{i}(length(ROI_coord_ERO{i,1})+1:length(ROI_coord_ERO{i,1})+1+length(ROI_coord_ERO{i,2}),j),score_coord{i}(length(ROI_coord_ERO{i,2})+1:end,j));
            p_val_ks{i,j}=[p_val_temp;p_val_temp2;p_val_temp3];
        end
    end
end

for i=1:length(ItiaList)
    regionName=ItiaList{i};
    for j=1:3
        if length(LinReg.(regionName).GoodBetas_AVG_final)==2
            [~,p_val_coord{i,j},~] = kstest2(ROI_coord_ERO{i,1}(:,j),ROI_coord_ERO{i,2}(:,j));
        elseif length(LinReg.(regionName).GoodBetas_AVG_final)==3
            [~,p_val_temp,~] = kstest2(ROI_coord_ERO{i,1}(:,j),ROI_coord_ERO{i,2}(:,j));
            [~,p_val_temp2,~] = kstest2(ROI_coord_ERO{i,1}(:,j),ROI_coord_ERO{i,3}(:,j));
            [~,p_val_temp3,~] = kstest2(ROI_coord_ERO{i,2}(:,j),ROI_coord_ERO{i,3}(:,j));
            p_val_coord{i,j}=[p_val_temp;p_val_temp2;p_val_temp3];
        end
    end
end

figure;
xplot=3;yplot=length(ItiaList);
for i=1:length(ItiaList)
    regionName=ItiaList{i};
    for j=1:3
        if length(LinReg.(regionName).GoodBetas_AVG_final)==2
            [~,p_val_ks{i,j},~] = kstest2(score_coord{i}(1:length(ROI_coord_ERO{i,1}),j),score_coord{i}(length(ROI_coord_ERO{i,1}+1):end,j));
        elseif length(LinReg.(regionName).GoodBetas_AVG_final)==3
            [~,p_val_temp,~] = kstest2(score_coord{i}(1:length(ROI_coord_ERO{i,1}),j),score_coord{i}(length(ROI_coord_ERO{i,1})+1:length(ROI_coord_ERO{i,1})+1+length(ROI_coord_ERO{i,2}),j));
            [~,p_val_temp2,~] = kstest2(score_coord{i}(1:length(ROI_coord_ERO{i,1}),j),score_coord{i}(length(ROI_coord_ERO{i,2})+1:end,j));
            [~,p_val_temp3,~] = kstest2(score_coord{i}(length(ROI_coord_ERO{i,1})+1:length(ROI_coord_ERO{i,1})+1+length(ROI_coord_ERO{i,2}),j),score_coord{i}(length(ROI_coord_ERO{i,2})+1:end,j));
            p_val_ks{i,j}=[p_val_temp;p_val_temp2;p_val_temp3];
        end
    end
end

Fighandle=figure;
set(Fighandle, 'Position', [100, 100, 960, 264]);
counter=1;xplot=3;yplot=1;counter2=1;
x = linspace(0.25,246/4,246);
i=1;
regionName=ItiaList{i};
idx_temp=LinReg.(regionName).KmeansIdx_select_AVG;
GoodBet_temp=LinReg.(regionName).GoodBetas_AVG_final;

ha = tight_subplot(2,1,[.01 .01],[.01 .01],[.01 .01]); 

counter=counter2;
if i==5
    counter=counter+1;
end
for j=1:length(GoodBet_temp)
    idx_temp2=find(idx_temp==GoodBet_temp(j));
    %subplot(xplot,yplot,counter+(j-1));
    axes(ha(j));
    for k=0:5
        meanToPlot=mean(LinReg.(regionName).ZS_AVG(idx_temp2,5+(k*40):40+(k*40)),1);
        %conf_95=confidence_intervals(LinReg.(regionName).ZS_AVG(idx_temp2,5+(k*40):40+(k*40)),95);
        std_95=std(LinReg.(regionName).ZS_AVG(idx_temp2,5+(k*40):40+(k*40)),1,1);
        hold on;%plot(x(5+(k*40):40+(k*40)),meanToPlot-mean(meanToPlot(1:4)),'color',colors{i}(j,:));title(num2str(length(idx_temp2)));
        H=shadedErrorBar(x(5+(k*40):40+(k*40)),meanToPlot-mean(meanToPlot(1:4)), std_95);axis([0 60 -3 4]);set(gca,'YTickLabel',[]);set(gca,'XTickLabel',[]);
        H.mainLine.Color=colors{i}(j,:);
        H.patch.FaceColor=colors{i}(j,:)/2;
        H.edge(1).Color=colors{i}(j,:)/2;
        H.edge(2).Color=colors{i}(j,:)/2;
    end
    hold off;
end
counter2=counter2+yplot;

clearvars i j k counter counter2 H idx_temp GoodBet_temp std_95 meanToPlot

%Laterality
PrismLaterality={};
Midline=306;imageSizeY = size(Template,1);
Index_ELO=strfind(Fish_list, 'ELO');
ELO_fish=find(not(cellfun('isempty', Index_ELO)));
ERO_fish=find(cellfun('isempty', Index_ELO));counter=1;
ELO_coord_mean=nan(length(ELO_fish),length(ItiaList),3);
ERO_coord_mean=nan(length(ERO_fish),length(ItiaList),3);
for i=1:length(ItiaList)    
    regionName=ItiaList{i};
    idx_temp=find(LinReg.(regionName).rsquared>0.1);
    ROIs_temp=ROIsPerBrain.(regionName).ROIs(idx_temp,:);
    Mask_limits(1,:)=min(ROIsPerBrain.(regionName).ROIs,[],1);
    Mask_limits(2,:)=max(ROIsPerBrain.(regionName).ROIs,[],1);
    Numbers=ROIsPerBrain.(regionName).Numbers;Numbers=[0 Numbers];
    counter=1;
    Prism_temp=nan(2,length(LinReg.(regionName).GoodBetas_AVG_final)*3);
    goodROIs_ELO_y=[];
    goodROIs_ERO_y=[];
    for j=1:length(LinReg.(regionName).GoodBetas_AVG_final)
        idx_goodROIs=find(LinReg.(regionName).KmeansIdx_select_AVG==LinReg.(regionName).GoodBetas_AVG_final(j));
        idx_goodROIs_correct=idx_temp(idx_goodROIs);
        idx_ELO=[];
        idx_ERO=[];
        counter_ELO=1;
        counter_ERO=1;
        for fishNb=1:length(Fish_list)
            if Index_ELO{fishNb}                
                idx_ELOtemp=find(Numbers(fishNb) < idx_goodROIs_correct & idx_goodROIs_correct <= Numbers(fishNb+1));
                if isempty(idx_ELO)
                    idx_ELO=idx_ELOtemp;
                else
                    idx_ELO=horzcat(idx_ELO, idx_ELOtemp);
                end
                ROIs_ELO=ROIsPerBrain.(regionName).ROIs(idx_goodROIs_correct(idx_ELOtemp),2);
                ELO_coord_mean(counter_ELO,i,j)=mean((ROIs_ELO-Midline)/(Mask_limits(2,2)-Midline));
                counter_ELO=counter_ELO+1;
            else
                idx_EROtemp=find(Numbers(fishNb) < idx_goodROIs_correct & idx_goodROIs_correct <= Numbers(fishNb+1));
                if isempty(idx_ERO)
                    idx_ERO=idx_EROtemp;
                else
                    idx_ERO=horzcat(idx_ERO, idx_EROtemp);
                end
                ROIs_ERO=ROIsPerBrain.(regionName).ROIs(idx_goodROIs_correct(idx_EROtemp),2);
                ERO_coord_mean(counter_ERO,i,j)=mean((ROIs_ERO-Midline)/(Mask_limits(2,2)-Midline));
                counter_ERO=counter_ERO+1;
            end
        end
        if length(idx_goodROIs)==(length(idx_ELO)+length(idx_ERO))
            goodROIs=ROIs_temp(idx_goodROIs,:);
            goodROIs_ELO_y{i,j}=(goodROIs(idx_ELO,2)-Midline)/(Mask_limits(2,2)-Midline);
            goodROIs_ERO_y{i,j}=(goodROIs(idx_ERO,2)-Midline)/(Mask_limits(2,2)-Midline);
            if ~isempty(goodROIs_ELO_y{i,j})
                Prism_temp(1,counter)=nanmean(goodROIs_ELO_y{i,j});
                Prism_temp(1,counter+1)=nanstd(goodROIs_ELO_y{i,j});
                Prism_temp(1,counter+2)=length(ELO_fish);
            end
            if ~isempty(goodROIs_ERO_y{i,j})
                Prism_temp(2,counter)=nanmean(goodROIs_ERO_y{i,j});
                Prism_temp(2,counter+1)=nanstd(goodROIs_ERO_y{i,j});
                Prism_temp(2,counter+2)=length(ERO_fish);
            end
            counter=counter+3;
        else
            break
        end
    end
    PrismLaterality{i}=Prism_temp;
end
clearvars goodROIs ROIs_temp regionName idx_temp Template2 idx_slice roi_nb xcoord ycoord circlePixels image_temp columnsInImage rowsInImage col corr_temp Fighdandle framerate GoodBet_temp H idx_g Model_cereb

Midline=306;imageSizeY = size(Template,1);
Index_ELO=strfind(Fish_list, 'ELO');
ELO_fish=find(not(cellfun('isempty', Index_ELO)));
ERO_fish=find(cellfun('isempty', Index_ELO));counter=1;
ELO_coord_count=nan(length(ELO_fish),length(ItiaList)*2,3);
ERO_coord_count=nan(length(ERO_fish),length(ItiaList)*2,3);
for i=1:length(ItiaList)
    regionName=ItiaList{i};
    idx_temp=find(LinReg.(regionName).rsquared>0.1);
    ROIs_temp=ROIsPerBrain.(regionName).ROIs(idx_temp,:);
    Mask_limits(1,:)=min(ROIsPerBrain.(regionName).ROIs,[],1);
    Mask_limits(2,:)=max(ROIsPerBrain.(regionName).ROIs,[],1);
    Numbers=ROIsPerBrain.(regionName).Numbers;Numbers=[0 Numbers];
    counter=1;    
    if i>10
        nb_clusters=length(LinReg.(regionName).GoodBeta_merge);        
    else
        nb_clusters=length(LinReg.(regionName).GoodBetas_AVG_final);        
    end    
    for j=1:nb_clusters
        if i>10            
            idx_goodROIs=find(LinReg.(regionName).KmeansIdx_merge==LinReg.(regionName).GoodBeta_merge(j));
            idx_goodROIs_correct=idx_temp(idx_goodROIs);
        else            
            idx_goodROIs=find(LinReg.(regionName).KmeansIdx_select_AVG==LinReg.(regionName).GoodBetas_AVG_final(j));
            idx_goodROIs_correct=idx_temp(idx_goodROIs);
        end
        counter_ELO=1;
        counter_ERO=1;
        for fishNb=1:length(Fish_list)
            if Index_ELO{fishNb}
                idx_ELOtemp=find(Numbers(fishNb) < idx_goodROIs_correct & idx_goodROIs_correct <= Numbers(fishNb+1));
                ROIs_ELO=ROIsPerBrain.(regionName).ROIs(idx_goodROIs_correct(idx_ELOtemp),2);
                ELO_coord_count(counter_ELO,(i*2)-1,j)=sum(ROIs_ELO>Midline);
                ELO_coord_count(counter_ELO,i*2,j)=sum(ROIs_ELO<Midline);
                counter_ELO=counter_ELO+1;
            else
                idx_EROtemp=find(Numbers(fishNb) < idx_goodROIs_correct & idx_goodROIs_correct <= Numbers(fishNb+1));
                ROIs_ERO=ROIsPerBrain.(regionName).ROIs(idx_goodROIs_correct(idx_EROtemp),2);
                ERO_coord_count(counter_ERO,(i*2)-1,j)=sum(ROIs_ERO>Midline);
                ERO_coord_count(counter_ERO,i*2,j)=sum(ROIs_ERO<Midline);
                counter_ERO=counter_ERO+1;
            end
        end
    end
end
clearvars goodROIs ROIs_temp regionName idx_temp Template2 idx_slice roi_nb xcoord ycoord circlePixels image_temp columnsInImage rowsInImage col corr_temp Fighdandle framerate GoodBet_temp H idx_g Model_cereb
clearvars counter_ERO ROIs_ERO ROIs_ELO counter_ELO Mask_limits Numbers idx_EROtemp idx_ELOtemp idx_goodROIs_correct idx_goodROIs counter counter2

test=squeeze(sum(ERO_coord_count,1));
test=squeeze(ERO_coord_count(:,:,1));

%Normalize count
Midline=306;imageSizeY = size(Template,1);
Index_ELO=strfind(Fish_list, 'ELO');
ELO_fish=find(not(cellfun('isempty', Index_ELO)));
ERO_fish=find(cellfun('isempty', Index_ELO));counter=1;
ELO_total_count=nan(length(ELO_fish),length(ItiaList)*2);
ERO_total_count=nan(length(ERO_fish),length(ItiaList)*2);
for i=1:length(ItiaList)
    regionName=ItiaList{i};
    idx_temp=find(LinReg.(regionName).rsquared>0.1);
    ROIs_temp=ROIsPerBrain.(regionName).ROIs;
    Numbers=ROIsPerBrain.(regionName).Numbers;Numbers=[0 Numbers];
    counter=1;
    goodROIs_ELO_y=[];
    goodROIs_ERO_y=[];    
    counter_ERO=1;
    counter_ELO=1;
    for fishNb=1:length(Fish_list)
        if Index_ELO{fishNb}
            ROIs_fish=ROIs_temp(Numbers(fishNb)+1:Numbers(fishNb+1),:);            
            ELO_total_count(counter_ELO,(i*2)-1)=sum(ROIs_fish(:,2)>Midline);
            ELO_total_count(counter_ELO,i*2)=sum(ROIs_fish(:,2)<Midline);
            counter_ELO=counter_ELO+1;
        else
            ROIs_fish=ROIs_temp(Numbers(fishNb)+1:Numbers(fishNb+1),:);                        
            ERO_total_count(counter_ERO,(i*2)-1)=sum(ROIs_fish(:,2)>Midline);
            ERO_total_count(counter_ERO,i*2)=sum(ROIs_fish(:,2)<Midline);
            counter_ERO=counter_ERO+1;
        end
    end
end

clearvars goodROIs ROIs_temp regionName idx_temp Template2 idx_slice roi_nb xcoord ycoord circlePixels image_temp columnsInImage rowsInImage col corr_temp Fighdandle framerate GoodBet_temp H idx_g Model_cereb
clearvars counter_ERO ROIs_ERO ROIs_ELO counter_ELO Mask_limits Numbers idx_EROtemp idx_ELOtemp idx_goodROIs_correct idx_goodROIs

test=squeeze(sum(ERO_total_count,1));


%Normalizing coordinates
 load('__MultiBrainRegions_final.mat', 'ROIsPerBrain');
 Norm_coord_ERO=ROI_coord_ERO;
for i=1:length(ItiaList)    
    regionName=ItiaList{i};
    Mask_limits=zeros(2,3);
    Mask=ROIsPerBrain.(regionName).ROIs;
    Mask_limits(1,:)=min(Mask,[],1);
    Mask_limits(2,:)=max(Mask,[],1);
    for j=1:3
        if ~isempty(ROI_coord_ERO{i,j});
            for axis=1:3
                Norm_coord_ERO{i,j}(:,axis)=(Norm_coord_ERO{i,j}(:,axis)-Mask_limits(1,axis))/(Mask_limits(2,axis)-Mask_limits(1,axis));
            end
        end
    end
end
    

for i=1:length(ItiaList)
    temp=[];
    for j=1:3
        if ~isempty(Norm_coord_ERO{i,j})
            if isempty(temp)
                temp=Norm_coord_ERO{i,j};
            else
                temp=vertcat(temp,Norm_coord_ERO{i,j});
            end 
        end
    end        
	[pca_coord_norm{i},score_coord_norm{i},~,~,explained_coord_norm{i},~] = pca(temp);
end

for i=1:length(ItiaList)
    regionName=ItiaList{i};
    for j=1:3
        if length(LinReg.(regionName).GoodBetas_AVG_final)==2
            [~,p_val_ks_norm{i,j},~] = kstest2(score_coord_norm{i}(1:length(Norm_coord_ERO{i,1}),j),score_coord_norm{i}(length(Norm_coord_ERO{i,1}+1):end,j));
        elseif length(LinReg.(regionName).GoodBetas_AVG_final)==3
            [~,p_val_temp,~] = kstest2(score_coord_norm{i}(1:length(Norm_coord_ERO{i,1}),j),score_coord_norm{i}(length(Norm_coord_ERO{i,1})+1:length(Norm_coord_ERO{i,1})+1+length(Norm_coord_ERO{i,2}),j));
            [~,p_val_temp2,~] = kstest2(score_coord_norm{i}(1:length(Norm_coord_ERO{i,1}),j),score_coord_norm{i}(length(Norm_coord_ERO{i,2})+1:end,j));
            [~,p_val_temp3,~] = kstest2(score_coord_norm{i}(length(Norm_coord_ERO{i,1})+1:length(Norm_coord_ERO{i,1})+1+length(Norm_coord_ERO{i,2}),j),score_coord_norm{i}(length(Norm_coord_ERO{i,2})+1:end,j));
            p_val_ks_norm{i,j}=[p_val_temp;p_val_temp2;p_val_temp3];
        end
    end
end

for i=1:length(ItiaList)
    regionName=ItiaList{i};
    for j=1:3
        if length(LinReg.(regionName).GoodBetas_AVG_final)==2
            [~,p_val_coord_norm{i,j},~] = kstest2(Norm_coord_ERO{i,1}(:,j),Norm_coord_ERO{i,2}(:,j));
        elseif length(LinReg.(regionName).GoodBetas_AVG_final)==3
            [~,p_val_temp,~] = kstest2(Norm_coord_ERO{i,1}(:,j),Norm_coord_ERO{i,2}(:,j));
            [~,p_val_temp2,~] = kstest2(Norm_coord_ERO{i,1}(:,j),Norm_coord_ERO{i,3}(:,j));
            [~,p_val_temp3,~] = kstest2(Norm_coord_ERO{i,2}(:,j),Norm_coord_ERO{i,3}(:,j));
            p_val_coord_norm{i,j}=[p_val_temp;p_val_temp2;p_val_temp3];
        end
    end
end

test=cell2mat(p_val_coord([1:8 10],:));
test1=cell2mat(p_val_coord(9,:));
test2=[test ; test1];
test2=test2(:);
FDR = mafdr(test2,'BHFDR', true);

Planes_KS_pvalues=[test ; test1];
[pVals_KS KS_idx]=sort(test2);
Multiple_comparison=length(Planes_KS_pvalues(:));
pval=0.05;
for i=1:length(pVals_KS)
    if pVals_KS(i) >= pval/(Multiple_comparison-(i-1))
        break
    end
end
pVals_KS(i:end)=nan;
KS_idx(isfinite(pVals_KS))=[];
Planes_KS_pvalues(KS_idx)=1;

Planes_KS_FDR=[test ; test1];
Planes_KS_FDR(FDR>=0.05)=1;

%Get Dorsal vs Ventral thalamus
regionName=ItiaList{1};
IndexC=strfind({Zbrain_Masks{:,2}}, regionName);
IndexC=find(not(cellfun('isempty', IndexC)));
dors_ROIs=Zbrain_Masks{IndexC(1),3};
ROI_temp=ROIsPerBrain.Thalamus.ROIs;
IsInDorsal=ismember(ROI_temp,dors_ROIs,'rows');
idx_dorsal=find(IsInDorsal==1);
clearvars ROI_temp IndexC dors_ROIs


Fighandle=figure;
set(Fighandle, 'Position', [100, 100, 960, 264]);
x = linspace(0.25,246/4,246);
regionName=ItiaList{1};
idx_temp=LinReg.(regionName).KmeansIdx_select_AVG;
GoodBet_temp=LinReg.(regionName).GoodBetas_AVG_final;
IsInDorsal_temp=IsInDorsal(find(LinReg.Thalamus.rsquared>0.1));
idx_dorsal=find(IsInDorsal_temp==1);
idx_ventral=find(IsInDorsal_temp==0);
for j=1:length(GoodBet_temp)
    counter=1+j-1;
    idx_temp2=find(idx_temp==GoodBet_temp(j));
    %idx_temp_d=find(idx_temp(idx_dorsal)==GoodBet_temp(j));
    idx_temp_db=ismember(idx_temp2,idx_dorsal);
    %idx_temp_v=find(idx_temp(idx_ventral)==GoodBet_temp(j));
    subplot(xplot,yplot,counter);
    for k=0:5
        meanToPlot=mean(LinReg.(regionName).ZS_AVG(idx_temp2(idx_temp_db),5+(k*40):40+(k*40)),1);        
        std_95=std(LinReg.(regionName).ZS_AVG(idx_temp2(idx_temp_db),5+(k*40):40+(k*40)),1,1);
        hold on;
        H=shadedErrorBar(x(5+(k*40):40+(k*40)),meanToPlot-mean(meanToPlot(1:4)), std_95);axis([0 60 -3 3]);
        H.mainLine.Color=colors{i}(j,:);
        H.patch.FaceColor=colors{i}(j,:)/2;
        H.edge(1).Color=colors{i}(j,:)/2;
        H.edge(2).Color=colors{i}(j,:)/2;
    end
    hold off;
    counter=4+j-1;
    subplot(xplot,yplot,counter);    
    for k=0:5
        meanToPlot=mean(LinReg.(regionName).ZS_AVG(idx_temp2(~idx_temp_db),5+(k*40):40+(k*40)),1);        
        std_95=std(LinReg.(regionName).ZS_AVG(idx_temp2(~idx_temp_db),5+(k*40):40+(k*40)),1,1);
        hold on;
        H=shadedErrorBar(x(5+(k*40):40+(k*40)),meanToPlot-mean(meanToPlot(1:4)), std_95);axis([0 60 -3 3]);
        H.mainLine.Color=colors{i}(j,:);
        H.patch.FaceColor=colors{i}(j,:)/2;
        H.edge(1).Color=colors{i}(j,:)/2;
        H.edge(2).Color=colors{i}(j,:)/2;
    end
    hold off;
end

idx_bet1=find(idx_temp==GoodBet_temp(1));
idx_bet2=find(idx_temp==GoodBet_temp(2));
idx_temp_db1=ismember(idx_bet1,idx_dorsal);
idx_temp_db2=ismember(idx_bet2,idx_dorsal);
subplot(xplot,yplot,3);bar([sum(idx_temp_db1) sum(idx_temp_db2)]);
subplot(xplot,yplot,6);bar([sum(~idx_temp_db1) sum(~idx_temp_db2)]);
idx_rsq_thal=find(LinReg.Thalamus.rsquared>0.1);
Numbers_fish=ROIsPerBrain.Thalamus.Numbers;
idx_fish_temp=zeros(size(idx_temp));
idx_fish_temp(find(idx_rsq_thal<=Numbers_fish(1)))=1;
for fish=2:length(Numbers_fish)    
    fish_idx=find(idx_rsq_thal<=Numbers_fish(fish) & idx_rsq_thal>Numbers_fish(fish-1));
    idx_fish_temp(fish_idx)=fish;    
end

distrib_thal=zeros(2,2,13);
for fish=1:13
    distrib_thal(1,1,fish)=length(find(idx_fish_temp(idx_bet1(idx_temp_db1))==fish));
    distrib_thal(2,1,fish)=length(find(idx_fish_temp(idx_bet2(idx_temp_db2))==fish));
    distrib_thal(1,2,fish)=length(find(idx_fish_temp(idx_bet1(~idx_temp_db1))==fish));
    distrib_thal(2,2,fish)=length(find(idx_fish_temp(idx_bet2(~idx_temp_db2))==fish));
end

distrib_thal_norm=zeros(2,2,13);
for fish=1:13
    distrib_thal_norm(1,1,fish)=length(find(idx_fish_temp(idx_bet1(idx_temp_db1))==fish))/length(find(idx_fish_temp(idx_bet1)==fish));
    distrib_thal_norm(2,1,fish)=length(find(idx_fish_temp(idx_bet2(idx_temp_db2))==fish))/length(find(idx_fish_temp(idx_bet2)==fish));
    distrib_thal_norm(1,2,fish)=length(find(idx_fish_temp(idx_bet1(~idx_temp_db1))==fish))/length(find(idx_fish_temp(idx_bet1)==fish));
    distrib_thal_norm(2,2,fish)=length(find(idx_fish_temp(idx_bet2(~idx_temp_db2))==fish))/length(find(idx_fish_temp(idx_bet2)==fish));
end

clearvars i j k counter counter2 H idx_temp GoodBet_temp std_95 meanToPlot

for region_nb=1:length(ItiaList)
    ForCSVExport=[];
    for cluster=1:3
        if ROI_coord_ERO{region_nb,cluster}
            if ForCSVExport
                CSV_temp=zeros(length(ROI_coord_ERO{region_nb,cluster}),4);
                CSV_temp(:,1:3)=ROI_coord_ERO{region_nb,cluster};
                CSV_temp(:,4)=cluster;
                ForCSVExport=[ForCSVExport;CSV_temp];
            else
                ForCSVExport=zeros(length(ROI_coord_ERO{region_nb,cluster}),4);
                ForCSVExport(:,1:3)=ROI_coord_ERO{region_nb,cluster};
                ForCSVExport(:,4)=cluster;
            end
        end
    end
    ForCSVExport(:,3)=ForCSVExport(:,3)*2;
    filename=strcat('ROIs_coord',ItiaList{region_nb},'.csv');    
    csvwrite(filename,ForCSVExport);
end

%Separate DT and VT
region_nb=1
ForCSVExport=[];
for cluster=1:3
    if ROI_coord_ERO{region_nb,cluster}
        if ForCSVExport
            CSV_temp=zeros(length(ROI_coord_ERO{region_nb,cluster}),4);
            CSV_temp(:,1:3)=ROI_coord_ERO{region_nb,cluster};
            CSV_temp(:,4)=cluster;
            ForCSVExport=[ForCSVExport;CSV_temp];
        else
            ForCSVExport=zeros(length(ROI_coord_ERO{region_nb,cluster}),4);
            ForCSVExport(:,1:3)=ROI_coord_ERO{region_nb,cluster};
            ForCSVExport(:,4)=cluster;
        end
    end
end
regionName=ItiaList{1};
IndexC=strfind({Zbrain_Masks{:,2}}, regionName);
IndexC=find(not(cellfun('isempty', IndexC)));
dors_ROIs=Zbrain_Masks{IndexC(1),3};
IsInDorsal=ismember(ForCSVExport(:,1:3),dors_ROIs,'rows');
idx_dorsal=find(IsInDorsal==1);
idx_ventral=find(IsInDorsal==0);
ForCSVExport(:,3)=ForCSVExport(:,3)*2;
filename=strcat('ROIs_coord_dorsal',ItiaList{region_nb},'.csv');
csvwrite(filename,ForCSVExport(idx_dorsal,:));
filename=strcat('ROIs_coord_ventral',ItiaList{region_nb},'.csv');
csvwrite(filename,ForCSVExport(idx_ventral,:));

%Graph figure with thick lines and no SD
counter=1;xplot=3;yplot=1;
x = linspace(0.25,246/4,246);
for region_nb=1:length(ItiaList);
    Fighandle=figure;
    set(Fighandle, 'Position', [100, 100, 600, 250]);
    counter=1;
    regionName=ItiaList{region_nb};
    idx_temp=LinReg.(regionName).KmeansIdx_select_AVG;
    GoodBet_temp=LinReg.(regionName).GoodBetas_AVG_final;
    if region_nb==5
        counter=counter+1;
    end
    for j=1:length(GoodBet_temp)
        idx_temp2=find(idx_temp==GoodBet_temp(j));
        for k=0:5
            meanToPlot=mean(LinReg.(regionName).ZS_AVG(idx_temp2,4+(k*41):40+(k*41)),1);
            hold on;
            plot(x(4+(k*40):40+(k*40)),meanToPlot-mean(meanToPlot(1:4)),'color',colors{region_nb}(j,:),'LineWidth',3);axis([0 60 -3 3]);            
            if (j==1 & k<5)
                rectangle('EdgeColor','none','FaceColor',[0.25, 0.25, 0.25, 0.4-k*0.06],'Position',[x(9+(k*40)) -3 1 7]);
            end
        end
        hold on;
    end
    hold off;
    
    print(Fighandle,regionName,'-dsvg','-r0');
    close all;
end

%Graph figure for PVL and Neuropil
Fighandle=figure;
set(Fighandle, 'Position', [100, 100, 600, 250]);
counter=1;xplot=3;yplot=1;
x = linspace(0.25,246/4,246);
region_nb=6;
counter=1;
regionName=ItiaList{region_nb};
ROI_temp=ROIsPerBrain.(regionName).ROIs;
idx_rsq=find(LinReg.(regionName).rsquared>0.1);
IndexC=strfind({Zbrain_Masks{:,2}}, regionName);
IndexC=find(not(cellfun('isempty', IndexC)));
PVL=Zbrain_Masks{IndexC(1),3};
IsInPVL=ismember(ROI_temp(idx_rsq,1:3),PVL,'rows');
idx_PVL=find(IsInPVL==1);
idx_neuropil=find(IsInPVL==0);
idx_temp=LinReg.(regionName).KmeansIdx_select_AVG;
GoodBet_temp=LinReg.(regionName).GoodBetas_AVG_final;
for j=1:length(GoodBet_temp)
    idx_temp_pvl=find(idx_temp(idx_PVL)==GoodBet_temp(j));
    idx_temp_neuropil=find(idx_temp(idx_neuropil)==GoodBet_temp(j));
    subplot(2,1,1);
    for k=0:5
        meanToPlot=mean(LinReg.(regionName).ZS_AVG(idx_PVL(idx_temp_pvl),4+(k*41):40+(k*41)),1);
        hold on;
        plot(x(4+(k*40):40+(k*40)),meanToPlot-mean(meanToPlot(1:4)),'color',colors{region_nb}(j,:),'LineWidth',3);axis([0 60 -3 3]);        
        rectangle('EdgeColor','none','FaceColor',[0.5, 0.5, 0.5, 0.2-k*0.04],'Position',[x(10+(k*40)) -3 1 7]);
    end
    hold off;
    subplot(2,1,2);
    for k=0:5
        meanToPlot=mean(LinReg.(regionName).ZS_AVG(idx_neuropil(idx_temp_neuropil),4+(k*41):40+(k*41)),1);
        plot(x(4+(k*40):40+(k*40)),meanToPlot-mean(meanToPlot(1:4)),'color',colors{region_nb}(j,:),'LineWidth',3);axis([0 60 -3 3]);hold on;
        rectangle('EdgeColor','none','FaceColor',[0.5, 0.5, 0.5, 0.2-k*0.04],'Position',[x(10+(k*40)) -3 1 7]);
    end
end


%Laterality Telencephalon inhibitory neurons
Telencephalon_ROIs=ROI_coord_ERO{5, 1};
idx_telen_caud=find(Telencephalon_ROIs(:,1)<1200);
prct_Telen=prctile(Telencephalon_ROIs(idx_telen_caud,2),[25 75]);
prct_Telen=[275 375];
idx_telen_25=find(Telencephalon_ROIs(idx_telen_caud,2) >= prct_Telen(1) & Telencephalon_ROIs(idx_telen_caud,2) < prct_Telen(2));
idx_telen_75=find(Telencephalon_ROIs(idx_telen_caud,2) < prct_Telen(1) | Telencephalon_ROIs(idx_telen_caud,2) >= prct_Telen(2));
Midline=306;
Telen_prct_count=zeros(2,2);
Telen_prct_count(1,1)=sum(Telencephalon_ROIs(idx_telen_caud(idx_telen_25),2)>Midline);
Telen_prct_count(1,2)=sum(Telencephalon_ROIs(idx_telen_caud(idx_telen_25),2)<Midline);
Telen_prct_count(2,1)=sum(Telencephalon_ROIs(idx_telen_caud(idx_telen_75),2)>Midline);
Telen_prct_count(2,2)=sum(Telencephalon_ROIs(idx_telen_caud(idx_telen_75),2)<Midline);

%Rhombomeres 5-7
IndexC=strfind({Zbrain_Masks{:,2}},'Rhombomere');
IndexC=find(not(cellfun('isempty', IndexC)));
IndexC=[222 223 224];
regionName='Rhombomeres_5-7';
Mask=[];
for j=IndexC
    if isempty(Mask)
        Mask=Zbrain_Masks{j,3};
    else
        Mask=vertcat(Mask,Zbrain_Masks{j,3});
    end
end

isInR57=[];
IsInBrainRegion=ismember(ROI_coord_ERO{10,1},Mask,'rows');
isInR57{1}=ROI_coord_ERO{10,1}(IsInBrainRegion,:);
IsInBrainRegion=ismember(ROI_coord_ERO{10,2},Mask,'rows');
isInR57{2}=ROI_coord_ERO{10,2}(IsInBrainRegion,:);
ROIs_temp=isInR57{1};
test=ismember(ROIs_temp,MON_Masks{1,3},'rows');
isInR57{1}=isInR57{1}(~test,:);
ROIs_temp=isInR57{2};
test=ismember(ROIs_temp,MON_Masks{1,3},'rows');
isInR57{2}=isInR57{2}(~test,:);

ForCSVExport=isInR57{1};
ForCSVExport(:,4)=1;
temp=isInR57{2};
temp(:,4)=2;
ForCSVExport=[ForCSVExport;temp];
ForCSVExport(:,3)=ForCSVExport(:,3)*2;
filename=strcat('__ROIs_coord_Rhomb57','.csv');
csvwrite(filename,ForCSVExport);

%Graph figure with thick lines and no SD
counter=1;xplot=3;yplot=1;
x = linspace(0.25,246/4,246);
for region_nb=10
    Fighandle=figure;
    set(Fighandle, 'Position', [100, 100, 600, 250]);
    counter=1;
    regionName=ItiaList{region_nb};
    idx_temp=LinReg.(regionName).KmeansIdx_select_AVG;
    GoodBet_temp=LinReg.(regionName).GoodBetas_AVG_final;
    if region_nb==5
        counter=counter+1;
    end
    for j=1:length(GoodBet_temp)
        idx_temp2=find(idx_temp==GoodBet_temp(j));
        
        for k=0:5
            meanToPlot=mean(LinReg.(regionName).ZS_AVG(idx_temp2,4+(k*41):40+(k*41)),1);
            hold on;
            plot(x(4+(k*40):40+(k*40)),meanToPlot-mean(meanToPlot(1:4)),'color',colors{region_nb}(j,:),'LineWidth',3);axis([0 60 -3 3]);            
            if (j==1 & k<5)
                rectangle('EdgeColor','none','FaceColor',[0.25, 0.25, 0.25, 0.4-k*0.06],'Position',[x(9+(k*40)) -3 1 7]);
            end
        end
        hold on;
    end
    hold off;
    
    print(Fighandle,regionName,'-dsvg','-r0');
    close all;
end

ItiaList={'Thalamus','Cerebellum','NucMLF','Semicircularis','Telencephalon','Tectum','Longitudinalis','Tegmentum','Habenula','Hindbrain','MON','Pretectum'};
i=12;
regionName=ItiaList{i};
progressbar;
for i=12
    progressbar(i/length(ItiaList),[]);
    regionName=ItiaList{i};
    Mask=[];
    IndexC=strfind({Zbrain_Masks{:,2}}, regionName);
    IndexC=find(not(cellfun('isempty', IndexC)));
    if i==5
        IndexC=294;
    end
    for j=IndexC
        if isempty(Mask)
            Mask=Zbrain_Masks{j,3};
        else
            Mask=vertcat(Mask,Zbrain_Masks{j,3});
        end
    end
    for fish_nb=1:length(Fish_list)
        progressbar([],fish_nb/length(Fish_list));
        if iscell(Fish_list)
            Fish_name=Fish_list{fish_nb};
        else
            Fish_name=num2str(Fish_list(fish_nb));
        end
        IndexC=strfind({MatFiles.name}, Fish_name);
        MatFiles_fish = find(not(cellfun('isempty', IndexC)));
        ROI_name=strsplit(Fish_name,'Fish2017');
        while iscell(ROI_name)
            ROI_name=ROI_name{1};
        end
        IndexC=strfind({ROIs.name},ROI_name);
        ROI_fish=find(not(cellfun('isempty', IndexC)));        
        ROI_fish=ROIs(ROI_fish).coord;ROI_fish(:,1:2)=round(ROI_fish(:,1:2));
        ROI_fish(:,3)=round(((ROI_fish(:,3)-1)*2)+24);
        if MatFiles_fish(1)==1            
            numbersForROIs=[1 [MatFiles(MatFiles_fish).number]];
        else            
            numbersForROIs=[MatFiles(MatFiles_fish(1)-1).number+1 [MatFiles(MatFiles_fish).number]];
        end
        IsInBrainRegion=[];        
        IsInBrainRegion=ismember(ROI_fish,Mask,'rows');
        PerBrainRegions(fish_nb).(regionName).ROIsCent=ROI_fish(IsInBrainRegion,:);        
        temp=ZS(numbersForROIs(1):numbersForROIs(end),:);IsInBrainRegion=temp(IsInBrainRegion,:);        
        PerBrainRegions(fish_nb).(regionName).ZS=IsInBrainRegion;        
    end
end


i=12;
regionName=ItiaList{i};
for j=1:length(PerBrainRegions)
    if j==1;
        ZS_Brain.(regionName)=PerBrainRegions(j).(regionName).ZS;
    else
        ZS_Brain.(regionName)=vertcat(ZS_Brain.(regionName),PerBrainRegions(j).(regionName).ZS);
    end
end

for i=12
    regionName=ItiaList{i};
    for j=1:length(PerBrainRegions)
        if j==1;
            ROIsPerBrain.(regionName).ROIs=PerBrainRegions(j).(regionName).ROIsCent;
            ROIsPerBrain.(regionName).Numbers=length(PerBrainRegions(j).(regionName).ROIsCent);
            temp=length(PerBrainRegions(j).(regionName).ROIsCent);
        else
            ROIsPerBrain.(regionName).ROIs=vertcat(ROIsPerBrain.(regionName).ROIs, PerBrainRegions(j).(regionName).ROIsCent); 
            temp=temp+length(PerBrainRegions(j).(regionName).ROIsCent);
            ROIsPerBrain.(regionName).Numbers=horzcat(ROIsPerBrain.(regionName).Numbers,temp);
        end
    end    
end
clearvars i j k ans Cmap coef coefficients counter Fish_name fish_nb GCaMP6 GoodBetas_nuc idxStart idx_Fish idx_temp idxKmeans IndexC image_temp IsInBrainRegion
clearvars IsInBrainRegion_good Mask MatFiles_fish mdl numbersForROIs_good numbersForROIs regionName ROI_fish ROI_fish_good ROI_name rsquared  temp



for j=12
    regionName=ItiaList{j};
    temp=ZS_Brain.(regionName);	
    coefficients=struct();
    rsquared=zeros(1,length(temp));
    parfor i=1:size(temp,1)
        mdl=fitlm(Stimuli',temp(i,:));%,'interactions');
        coefficients(i).coef=mdl.Coefficients;
        rsquared(i)=mdl.Rsquared.Adjusted;
    end    
    LinReg.(regionName).coef=coefficients;
    LinReg.(regionName).rsquared=rsquared;    
end

%Kmeans of all
options = statset('UseParallel',1); 
for j=12
    regionName=ItiaList{j};
    idx_temp=find(LinReg.(regionName).rsquared>0.1);
    LinReg.(regionName).ZS_rsq=ZS_Brain.(regionName)(idx_temp,:);
    [idxKmeans Cmap]=kmeans(LinReg.(regionName).ZS_rsq,10,'Options',options,'Distance','cityblock','Replicates',3,'MaxIter',1000,'Display','final');
    LinReg.(regionName).KmeansCenter=Cmap;
    LinReg.(regionName).KmeansIdx=idxKmeans;
end

%Get the good ones
for j=12
    regionName=ItiaList{j};
    [~,LinReg.(regionName).GoodBetas]=Test_Regress( LinReg.(regionName).KmeansCenter,Stimuli,LinReg.(regionName).KmeansIdx,0.5);
end

%Correlate back to Kmeans to remove crap
Threshold=0.4;
for i=12
    regionName=ItiaList{i};
    temp=LinReg.(regionName).ZS_rsq;
    idx_temp=LinReg.(regionName).KmeansIdx;
    GoodBet_temp=LinReg.(regionName).GoodBetas;
    for j=1:length(GoodBet_temp)
        idx_g=find(idx_temp==GoodBet_temp(j));
        temp_g=temp(idx_g,:);
        corr_temp=zeros(size(idx_g));
        parfor k=1:length(idx_g)        
            temp=corrcoef(LinReg.(regionName).KmeansCenter(GoodBet_temp(j),:), temp_g(k,:));
            corr_temp(k)=temp(1,2);
        end        
        idx_temp(idx_g(find(corr_temp<=Threshold)))=0;
    end
    LinReg.(regionName).KmeansIdx_select=idx_temp;    
end
clearvars i j k temp ans Cmap coeff coefficients explained i Hindbrain_Mask idx_temp idxKmeans ish_nb 

figure;
counter=1;xplot=length(ItiaList);yplot=4;counter2=1;
for i=12
    regionName=ItiaList{i};
    idx_temp=LinReg.(regionName).KmeansIdx_select;
    GoodBet_temp=LinReg.(regionName).GoodBetas;
    counter=counter2;
    for j=1:length(GoodBet_temp)
        idx_temp2=find(idx_temp==GoodBet_temp(j));
        subplot(xplot,yplot,counter+(j-1));plot(mean(LinReg.(regionName).ZS_rsq(idx_temp2,:),1));
    end
    counter2=counter2+yplot;
end

%Building the AVG across 3 presentations
for j=12
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

GoodBet_temp=LinReg.(regionName).GoodBetas;
KmeansIdx_merge=LinReg.(regionName).KmeansIdx_select;
idx_temp=ismember(KmeansIdx_merge,GoodBet_temp([2 3]));
KmeansIdx_merge(idx_temp)=GoodBet_temp(3);
LinReg.(regionName).GoodBeta_merge=GoodBet_temp([3 1]);
LinReg.(regionName).KmeansIdx_merge=KmeansIdx_merge;

i=12;colors{i}=[0.14 1 0.14; 0.7 0.4 1];regionName=ItiaList{i};
%Graph figure with thick lines and no SD
counter=1;xplot=3;yplot=1;
x = linspace(0.25,246/4,246);
region_nb=12;
Fighandle=figure;
set(Fighandle, 'Position', [100, 100, 600, 250]);
counter=1;
regionName=ItiaList{region_nb};
idx_temp=LinReg.(regionName).KmeansIdx_merge;
GoodBet_temp=LinReg.(regionName).GoodBeta_merge;
if region_nb==5
    counter=counter+1;
end
for j=1:length(GoodBet_temp)
    idx_temp2=find(idx_temp==GoodBet_temp(j));
    for k=0:5
        meanToPlot=mean(LinReg.(regionName).ZS_AVG(idx_temp2,4+(k*41):40+(k*41)),1);
        hold on;
        plot(x(4+(k*40):40+(k*40)),meanToPlot-mean(meanToPlot(1:4)),'color',colors{region_nb}(j,:),'LineWidth',3);axis([0 60 -3 3]);
        if (j==1 & k<5)
            rectangle('EdgeColor','none','FaceColor',[0.25, 0.25, 0.25, 0.4-k*0.06],'Position',[x(9+(k*40)) -3 1 7]);
        end
    end
    hold on;
end
hold off;

print(Fighandle,regionName,'-dsvg','-r0');
close all;


for i=12
    regionName=ItiaList{i};
    idx_temp=find(LinReg.(regionName).rsquared>0.1);
    ROIs_temp=ROIsPerBrain.(regionName).ROIs(idx_temp,:);
    Numbers=ROIsPerBrain.(regionName).Numbers;Numbers=[0 Numbers];
    for j=1:length(LinReg.(regionName).GoodBeta_merge)
        idx_goodROIs=find(LinReg.(regionName).KmeansIdx_merge==LinReg.(regionName).GoodBeta_merge(j));
        idx_goodROIs_correct=idx_temp(idx_goodROIs);
        idx_ELO=[];
        idx_ERO=[];
        for fishNb=1:length(Fish_list)
            if Index_ELO{fishNb}                
                idx_ELOtemp=find(Numbers(fishNb) < idx_goodROIs_correct & idx_goodROIs_correct <= Numbers(fishNb+1));
                if isempty(idx_ELO)
                    idx_ELO=idx_ELOtemp;
                else
                    idx_ELO=horzcat(idx_ELO, idx_ELOtemp);
                end
            else
                idx_EROtemp=find(Numbers(fishNb) < idx_goodROIs_correct & idx_goodROIs_correct <= Numbers(fishNb+1));
                if isempty(idx_ERO)
                    idx_ERO=idx_EROtemp;
                else
                    idx_ERO=horzcat(idx_ERO, idx_EROtemp);
                end
            end
        end
        goodROIs=ROIs_temp(idx_goodROIs,:); 
        ROI_coord_ELO{i,j}=goodROIs(idx_ELO,:);
        ROI_coord_ERO{i,j}=goodROIs(idx_ERO,:);    
    end
end
clearvars goodROIs ROIs_temp regionName idx_temp Template2 idx_slice roi_nb xcoord ycoord circlePixels image_temp columnsInImage rowsInImage col corr_temp Fighdandle framerate GoodBet_temp H idx_g Model_cereb


for region_nb=12
    ForCSVExport=[];
    for cluster=1:3
        if ROI_coord_ERO{region_nb,cluster}
            if ForCSVExport
                CSV_temp=zeros(length(ROI_coord_ERO{region_nb,cluster}),4);
                CSV_temp(:,1:3)=ROI_coord_ERO{region_nb,cluster};
                CSV_temp(:,4)=cluster;
                ForCSVExport=[ForCSVExport;CSV_temp];
            else
                ForCSVExport=zeros(length(ROI_coord_ERO{region_nb,cluster}),4);
                ForCSVExport(:,1:3)=ROI_coord_ERO{region_nb,cluster};
                ForCSVExport(:,4)=cluster;
            end
        end
    end
    ForCSVExport(:,3)=ForCSVExport(:,3)*2;
    filename=strcat('__ROIs_coord',ItiaList{region_nb},'.csv');    
    csvwrite(filename,ForCSVExport);
end


