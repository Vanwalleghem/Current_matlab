CSV_Files=dir('*.csv');
ROIs=struct();
for i=1:length(CSV_Files);
    temp=csvread(CSV_Files(i).name,0);          
    ROIs(i).coord=temp(:,1:4);
    ROIs(i).idx=temp(:,4);    
end
clearvars i temp CSV_Files Fishname

i=1;ROI_pool=ROIs(i).coord;
for i=2:length(ROIs);
    if i==8
        ROI_temp=ROIs(i).coord;
        ROI_temp(ROI_temp(:,4)==2,4)=3;
        ROI_temp(ROI_temp(:,4)==1,4)=2;
        ROI_pool=[ROI_pool; ROI_temp];
    else
        ROI_pool=[ROI_pool; ROIs(i).coord];
    end
end

for i=1:3
    ROI_temp=ROI_pool(find(ROI_pool(:,4)==i),:);
    for prop_ROI=20:20:100
        CSV_temp=datasample(ROI_temp,round(size(ROI_temp,1)*(prop_ROI/100)));
        filename=strcat('ROIsWB_clust',num2str(i),'_prop',num2str(prop_ROI),'.csv');
        csvwrite(filename,CSV_temp);
    end
end

i=1;prop_ROI=10;
ROI_temp=ROI_pool(find(ROI_pool(:,4)==i),:);
CSV_temp=datasample(ROI_temp,round(size(ROI_temp,1)*(prop_ROI/100)));
        filename=strcat('ROIsWB_clust',num2str(i),'_prop',num2str(prop_ROI),'.csv');
        csvwrite(filename,CSV_temp);
        
        
figure;histogram(