%% 
%Get all the CaImAn outputs
HDF5list1=dir('*4D2b_new.hdf5');


%Load the functional data
CaImAnData=struct();
for i=1:length(HDF5list1)
    CaImAnData(i).idx_components=h5read(HDF5list1(i).name,'/estimates/idx_components');
    temp=h5read(HDF5list1(i).name,'/estimates/C');
    %CaImAnData(i).CalciumTraces=h5read(HDF5list1(i).name,'/estimates/C');
    CaImAnData(i).DF=DeltaF2(temp(:,1+CaImAnData(i).idx_components)',70,21);
end

%CaImAnData=rmfield(CaImAnData,'CalciumTraces');


%% Get ROIs

% Ensure parallel pool is properly initialized
if isempty(gcp('nocreate'))
    parpool('local');
end

% Iterate through each file in HDF5list1
for i = 1:length(HDF5list1)
    % Load data
    data = h5read(HDF5list1(i).name, '/estimates/A/data');
    indptr = h5read(HDF5list1(i).name, '/estimates/A/indptr');
    indices = h5read(HDF5list1(i).name, '/estimates/A/indices');
    shape  = h5read(HDF5list1(i).name, '/estimates/A/shape');
    dims = h5read(HDF5list1(i).name, '/dims');  % Assuming 3D data

    % Convert dimensions to row vector
    dims = dims(:)';

    % Convert the csc matrix into a matlab sparse
    A_sparse=sparse(shape(1),shape(2));
    for k=1:length(indptr)-1
        ij=indptr(k)+1:indptr(k+1);
        y=indices(ij)+1;
        A_sparse(y,k)=data(ij);
    end

    % Initialize ROI centroid storage
    idx_temp = CaImAnData(i).idx_components+1;
    ROI_centroid = zeros(length(idx_temp), 3);  % For 3D data

    % Parallel processing each ROI
    Img_sum = zeros(dims);
    parfor ROI_nb = 1:length(idx_temp)        
        ROI_idx=idx_temp(ROI_nb);

        temp=A_sparse(:,ROI_idx);
        Img=reshape(full(temp),dims);       
        s = regionprops(Img > 0, Img, 'WeightedCentroid');
        if ~isempty(s)
            ROI_centroid(ROI_nb, :) = s.WeightedCentroid;
        else
            ROI_centroid(ROI_nb, :) = [NaN, NaN, NaN];
        end
        Img_sum=Img_sum+Img;
    end
    CaImAnData(i).ROI_centroid = ROI_centroid;
    CaImAnData(i).ImgROIs=Img_sum;
end

% Shutdown parallel pool after processing
delete(gcp('nocreate'));


%figure;imagesc(squeeze(max(Img,[],3)))
%% Visualisation of ROI
%Load the tiff file, can be max intensity or what not
TiffImg = tiffreadVolume("GV_20200729_fish2_ENSGF_5DPF_range110_step5_exposure21_power60_Max.tif");
TiffImg = permute(TiffImg,[3 1 2]);
%TiffImgMax = zeros(size(TiffImg,1),size(TiffImg,2),size(TiffImg,3)/10);
%for i = 1:size(TiffImgMax,3)
%    TiffImgMax(:,:,i) = max(TiffImg(:,:,i:size(TiffImgMax,3)/10:size(TiffImg,3)),[],3);
%end
StackSliderAndROI(TiffImg,ROI_centroid)

%% figure with slider, need a separate function .m file
fig = uifigure;
sld = uislider(fig);
StackSlider(TiffImg(:,:,1:29))

%% Compile the data and index it

MatFiles_names={HDF5list1.name};

%Design a regular expression that match your naming scheme (https://regexr.com/)
%fin = cellfun(@(x)regexp(x,'fish(\d+)_ENS(.*)_(\d)DPF.*_slice(\d+)','tokens'), MatFiles_names, 'UniformOutput', false);
fin = cellfun(@(x)regexp(x,'fish(\d+)_ENS(.*)_(\d)DPF.*\.hdf5','tokens'), MatFiles_names, 'UniformOutput', false);

names=zeros(3,length(fin));
for i=1:length(fin)
    names(1, i) = str2num(fin{i}{1}{1});  % Fish number
    names(2, i) = str2num(fin{i}{1}{3});  % Days post fertilization
    
    if strcmp(fin{i}{1}{2},'GF')
        names(3,i)=2;
    elseif strcmp(fin{i}{1}{2},'Fed')
        names(3,i)=3;
    else
        names(3,i)=1;
    end      
end
FishList=unique(names);
clearvars i Fish fin
%% Need to get limits from the template binary images
Tiff_Templates=dir('*.tif');
Limits=zeros(4,length(Tiff_Templates));
for ij=1:length(Tiff_Templates)
    info = imfinfo(Tiff_Templates(ij).name);
    numberOfPages = length(info);
    ik = 1;thisPage = imread(Tiff_Templates(ij).name, ik);thisPage=logical(thisPage);
    Limits_temp=regionprops(thisPage,'BoundingBox');Limits_temp=Limits_temp.BoundingBox;
    [~,idx_temp]=max(Limits_temp(3:4));Limits_temp(2+idx_temp)=Inf;Limits_temp(1:2)=Limits_temp(1:2)-20;
    %% 
    idx_temp=setdiff([1 2],idx_temp);Limits_temp(2+idx_temp)=Limits_temp(2+idx_temp)+40;
    Limits(:,ij)=Limits_temp;
end
%% Check that ROIs are within the gut

%This did not work at all
% %Need the MatFiles_ordered from the suite2p analysis to have the proper
% %order of the files
% 
% MatFiles_ordered_names={MatFiles_ordered.name};
% 
% fin = cellfun(@(x)regexp(x,'fish(\d+)_ENS(.*)_(\d)DPF','tokens'), MatFiles_ordered_names, 'UniformOutput', false);
% names_gutLim=zeros(3,length(fin));
% for i=1:length(fin)    
%     names_gutLim(1,i)=str2num(fin{i}{1}{1});names_gutLim(2,i)=str2num(fin{i}{1}{3});
%     if strcmp(fin{i}{1}{2},'GF')
%         names_gutLim(3,i)=2;
%     elseif strcmp(fin{i}{1}{2},'Fed')
%         names_gutLim(3,i)=3;
%     else
%         names_gutLim(3,i)=1;
%     end      
% end
% clearvars i Fish fin
% 
% %Get the limits from the MIPs
% fid = fopen('Gut limits');
% Limits=inf(2,length(MatFiles_ordered_names));
% counter=1;
% while ~feof(fid) % feof(fid) is true when the file ends
%       textLineEntry = fgetl(fid); % read one line
%       if regexp(textLineEntry,'x')
%         Limits(1,counter)=str2num(textLineEntry(3:end));
%       elseif regexp(textLineEntry,'y')
%         Limits(2,counter)=str2num(textLineEntry(3:end));
%       else
%           warning('error');
%       end
%       counter=counter+1;
% end
% fclose(fid); % close the file

Tiff_Templates_names={Tiff_Templates.name};
%Now to check all the slices against the MIPs
for ij=1:length(MatFiles_names)
    name_temp=MatFiles_names{ij};
    %name_temp=split(name_temp,'_4D_');name_temp=name_temp{1};
    %idx_temp=find(contains(MatFiles_ordered_names,name_temp));
    name_temp=split(name_temp,'_range');name_temp=name_temp{1};
    idx_temp=find(contains(Tiff_Templates_names,name_temp));
    ROI_temp=CaImAnData(ij).ROI_centroid;
    CaImAnData(ij).idx_Gut=(ROI_temp(:,1)<Limits(1,idx_temp)+Limits(3,idx_temp) & ROI_temp(:,1)>Limits(1,idx_temp)) & (ROI_temp(:,2)<Limits(2,idx_temp)+Limits(4,idx_temp) & ROI_temp(:,2)>Limits(2,idx_temp));    
end


%%
% Define template file patterns and initialize limits array
Tiff_Templates = dir('*.tif');
Limits = zeros(6, length(Tiff_Templates));  % Adjust for 3D (BoundingBox has 6 elements in 3D)

for ij = 1:length(Tiff_Templates)
    info = imfinfo(Tiff_Templates(ij).name);
    numberOfPages = length(info);

    % Process the 3D volume in chunks to reduce memory usage
    chunkSize = 100;  % Define the chunk size
    numChunks = ceil(numberOfPages / chunkSize);

    % Initialize overall bounding box limits
    overallLimits = [Inf, Inf, Inf, 0, 0, 0];

    for chunkIdx = 1:numChunks
        startSlice = (chunkIdx - 1) * chunkSize + 1;
        endSlice = min(chunkIdx * chunkSize, numberOfPages);

        % Initialize 3D binary image for the current chunk
        thisVolume = false(info(1).Height, info(1).Width, endSlice - startSlice + 1);

        % Read each slice of the current chunk of the 3D TIFF
        for ik = startSlice:endSlice
            thisVolume(:,:,ik - startSlice + 1) = logical(imread(Tiff_Templates(ij).name, ik));
        end

        % Compute the bounding box using regionprops3 for the current chunk
        Limits_temp = regionprops3(thisVolume, 'BoundingBox');
        Limits_temp = Limits_temp.BoundingBox;

        % Update overall bounding box limits
        overallLimits(1:3) = min(overallLimits(1:3), Limits_temp(1:3));
        overallLimits(4:6) = max(overallLimits(4:6), Limits_temp(4:6) + Limits_temp(1:3));
    end

    % Adjust bounding box
    [~, idx_temp] = max(overallLimits(4:6));  % For 3D, max dimension indices 4:6
    overallLimits(3+idx_temp) = Inf;
    overallLimits(1:3) = overallLimits(1:3) - 20;

    idx_temp = setdiff([1 2 3], idx_temp);  % Adjust the other dimensions
    overallLimits(3+idx_temp) = overallLimits(3+idx_temp) + 40;

    Limits(:, ij) = overallLimits(:);  % Store limits
end

% Check that ROIs are within the gut
Tiff_Templates_names = {Tiff_Templates.name};

for ij = 1:length(MatFiles_names)
    name_temp = MatFiles_names{ij};
    name_temp = split(name_temp, '_range'); name_temp = name_temp{1};
    idx_temp = find(contains(Tiff_Templates_names, name_temp));

    ROI_temp = CaImAnData(ij).ROI_centroid;

    % Check if ROIs are within the gut limits in 3D
    CaImAnData(ij).idx_Gut = ...
        (ROI_temp(:,1) < Limits(1, idx_temp) + Limits(4, idx_temp) & ROI_temp(:,1) > Limits(1, idx_temp)) & ...
        (ROI_temp(:,2) < Limits(2, idx_temp) + Limits(5, idx_temp) & ROI_temp(:,2) > Limits(2, idx_temp)) & ...
        (ROI_temp(:,3) < Limits(3, idx_temp) + Limits(6, idx_temp) & ROI_temp(:,3) > Limits(3, idx_temp));
end

%% Visualize
for ij = 1:length(MatFiles_names)
    figure;
    scatter3(CaImAnData(ij).ROI_centroid(:,1), CaImAnData(ij).ROI_centroid(:,2), CaImAnData(ij).ROI_centroid(:,3), 'filled');
    title(['ROI Centroids for file ', num2str(ij)]);
    xlabel('X');
    ylabel('Y');
    zlabel('Z');
end

%% Compile relevant data in cells
DF_guts={};
for treat=1:3
    for dpf=3:7        
        for fish_nb=1:7
            idx_temp=find(ismember(names(1:3,:)',[fish_nb,dpf,treat],'rows'));
            for data_nb=1:length(idx_temp)
                if data_nb==1
                    DF_guts{treat,dpf-2,fish_nb}=CaImAnData(idx_temp(data_nb)).DF(CaImAnData(idx_temp(data_nb)).idx_Gut,:);
                else
                    DF_guts{treat,dpf-2,fish_nb}=vertcat(DF_guts{treat,dpf-2,fish_nb},CaImAnData(idx_temp(data_nb)).DF(CaImAnData(idx_temp(data_nb)).idx_Gut,:));
                end
            end
        end
    end    
end
%%

DF_guts={};
for treat=1:3
    for dpf=3:7        
        for fish_nb=1:7
            idx_temp=find(ismember(names(1:3,:)',[fish_nb,dpf,treat],'rows'));
            for data_nb=1:length(idx_temp)
                if data_nb==1
                    DF_guts{treat,dpf-2,fish_nb}=CaImAnData(idx_temp(data_nb)).DF(CaImAnData(idx_temp(data_nb)).idx_Gut,:);
                else
                    DF_guts{treat,dpf-2,fish_nb}=vertcat(DF_guts{treat,dpf-2,fish_nb},CaImAnData(idx_temp(data_nb)).DF(CaImAnData(idx_temp(data_nb)).idx_Gut,:));
                end
            end
        end
    end    
end

%% See rastermaps
Fighandle=figure;
set(Fighandle, 'Position', [100, 100, 1200, 400]);
set(findall(Fighandle,'type','text'),'fontSize',12,'fontWeight','bold','FontName','Arial')
xplot=7;yplot=5;
ha = tight_subplot(yplot,xplot,[.01 .01]*5,[.01 .01]*5,[.01 .01]*5);
counter=1;
treat=3;
for dpf=3:7
    counter=1+(dpf-3)*7;    
    for fish_nb=1:7
        ZS_temp=DF_guts{treat,dpf-2,fish_nb};        
        axes(ha(counter));
        imagesc(ZS_temp,[0 1]);colormap hot
        set(gca,'XTickLabel',[]);set(gca,'YTickLabel',[]);
        title(strcat(num2str(size(DF_guts{treat,dpf-2,fish_nb},1)),' fish nb ',num2str(fish_nb),' DPF ',num2str(dpf)));
        counter=counter+1;
    end
end

%%
ops.nCall=[30 100];
ops.useGPU=0;
smooth_ZS_fish={};
results_ENS=struct();
results_ENS_raw=struct();
for treat=1:3
    for dpf=3:7        
        for fish_nb=1:7
            ZS_temp=DF_guts{treat,dpf-2,fish_nb};
            if ~isempty(ZS_temp)
                ops.iPC=1:min([size(ZS_temp,1),size(ZS_temp,2)]);
                [isort1, isort2, Sm] = mapTmap(ZS_temp,ops);
                smooth_ZS_fish{treat,dpf-2,fish_nb,1}=isort1;
                smooth_ZS_fish{treat,dpf-2,fish_nb,2}=isort2;
                smooth_ZS_fish{treat,dpf-2,fish_nb,3}=Sm;
                pks_smooth=cell(1,size(Sm,1));
                locs_smooth=cell(1,size(Sm,1));                
                Periodicity_smooth=nan(size(Sm,1),1);
                pks_raw=cell(1,size(Sm,1));
                locs_raw=cell(1,size(Sm,1));                
                Periodicity_raw=nan(size(Sm,1),1);
                for roi_nb=1:size(ZS_temp,1)
                    [pks_smooth{roi_nb},locs_smooth{roi_nb}] = findpeaks(Sm(roi_nb,:),'MinPeakProminence',0.1,'MinPeakDistance',5);
                    [pks_raw{roi_nb},locs_raw{roi_nb}] = findpeaks(ZS_temp(roi_nb,:),'MinPeakProminence',0.1,'MinPeakDistance',5);
                    if length(locs_smooth{roi_nb})>5
                        locs_temp=locs_smooth{roi_nb}(2:end)-locs_smooth{roi_nb}(1:end-1);
                        locs_temp=locs_temp/2;
                        Periodicity_smooth(roi_nb)=mean(locs_temp);
                    end
                    if length(locs_raw{roi_nb})>5
                        locs_temp=locs_raw{roi_nb}(2:end)-locs_raw{roi_nb}(1:end-1);
                        locs_temp=locs_temp/2;
                        Periodicity_raw(roi_nb)=mean(locs_temp);
                    end
                end
                results_ENS(treat,dpf-2,fish_nb).peaks=pks_smooth;
                results_ENS(treat,dpf-2,fish_nb).locs=locs_smooth;
                results_ENS(treat,dpf-2,fish_nb).Periodicity=Periodicity_smooth;
                results_ENS_raw(treat,dpf-2,fish_nb).peaks=pks_raw;
                results_ENS_raw(treat,dpf-2,fish_nb).locs=locs_raw;
                results_ENS_raw(treat,dpf-2,fish_nb).Periodicity=Periodicity_raw;
            end
        end
    end
end
%% Get Data ready for Prism - Periodicity

PrismTemp=nan(3,5*7);
for treat=1:3
    for dpf=3:7
        for fish_nb=1:7
            PrismTemp(treat,fish_nb+(dpf-3)*7)=nanmean(results_ENS(treat,dpf-2,fish_nb).Periodicity);
        end
    end
end
PrismTemp(PrismTemp==0)=nan;

%% Periodicity distribution
edge=[0:3:100];
treat=3;
h=nan(length(edge)-1,5*7);
for dpf=3:7
    for fish_nb=1:7
        temp=results_ENS_raw(treat,dpf-2,fish_nb).Periodicity;
        h(:,fish_nb+(dpf-3)*7)=histcounts(temp(~isnan(temp)),edge,'Normalization','probability')';
    end
end

%% Same for raw

% raw prism data
PrismTemp=nan(3,5*7);
for treat=1:3
    for dpf=3:7
        for fish_nb=1:7
            PrismTemp(treat,fish_nb+(dpf-3)*7)=nanmean(results_ENS_raw(treat,dpf-2,fish_nb).Periodicity);
        end
    end
end
PrismTemp(PrismTemp==0)=nan;

edge=[0:3:100];
treat=1;
h3=nan(length(edge)-1,5*7);
for dpf=3:7
    for fish_nb=1:7
        temp=results_ENS_raw(treat,dpf-2,fish_nb).Periodicity;
        h3(:,fish_nb+(dpf-3)*7)=histcounts(temp(~isnan(temp)),edge,'Normalization','probability')';
    end
end


%% peak analysis

edge=[0:2:60];
h1=nan(length(edge)-1,5*7);
treat=1;
for dpf=3:7
    for fish_nb=1:7
        peak_temp_nb=[];
        peak_temp=results_ENS_raw(treat,dpf-2,fish_nb).peaks;
        for i=1:length(peak_temp)
            peak_temp_nb(i)=length(peak_temp{i});
        end
        peak_temp_nb=peak_temp_nb(~isnan(peak_temp_nb));
        h1(:,fish_nb+(dpf-3)*7)=histcounts(peak_temp_nb,edge,'Normalization','probability')';
    end
end
h1=cumsum(h1,1);

%lots of inactive neurons
treat=1;

edge=[0:2:80];
h1=nan(length(edge)-1,5*7);
for dpf=3:7
    for fish_nb=1:7
        peak_temp_nb=[];
        peak_temp=results_ENS_raw(treat,dpf-2,fish_nb).peaks;
        for i=1:length(peak_temp)
            peak_temp_nb(i)=length(peak_temp{i});
        end
        peak_temp_nb=peak_temp_nb(~isnan(peak_temp_nb));
        peak_temp_nb(peak_temp_nb==0)=[];
        h1(:,fish_nb+(dpf-3)*7)=histcounts(peak_temp_nb,edge,'Normalization','probability')';
    end
end
h1=cumsum(h1,1);
h1=[h1;ones(1,size(h1,2))];

% raw prism data
PrismTemp=nan(3,5*7);
for treat=1:3
    for dpf=3:7
        for fish_nb=1:7
            peak_temp_nb=[];
            peak_temp=results_ENS_raw(treat,dpf-2,fish_nb).peaks;
            for i=1:length(peak_temp)
                peak_temp_nb(i)=length(peak_temp{i});
            end
            PrismTemp(treat,fish_nb+(dpf-3)*7)=nanmean(peak_temp_nb(peak_temp_nb>0));
        end
    end
end
PrismTemp(PrismTemp==0)=nan;


% raw prism data
PrismTemp=nan(3,5*7);
for treat=1:3
    for dpf=3:7
        for fish_nb=1:7
            peak_temp_nb=[];
            peak_temp=results_ENS(treat,dpf-2,fish_nb).peaks;
            for i=1:length(peak_temp)
                peak_temp_nb(i)=length(peak_temp{i});
            end
            PrismTemp(treat,fish_nb+(dpf-3)*7)=nanmean(peak_temp_nb(peak_temp_nb>0));
        end
    end
end
PrismTemp(PrismTemp==0)=nan;


%% Cluster quickly each fish

K_clusters_5=cell([size(DF_guts) 2]);
for treat=1:3
    for dpf=3:7
        for fish_nb=1:7
            ZS_temp=DF_guts{treat,dpf-2,fish_nb};
            if ~isempty(ZS_temp)                
                [K_clusters_5{treat,dpf-2,fish_nb,1},K_clusters_5{treat,dpf-2,fish_nb,2}]=kmeans(ZS_temp,5,"Distance","cityblock");
            end
        end
    end
end


% See rastermaps
Fighandle=figure;
set(Fighandle, 'Position', [100, 100, 1200, 400]);
set(findall(Fighandle,'type','text'),'fontSize',12,'fontWeight','bold','FontName','Arial')
xplot=7;yplot=5;
ha = tight_subplot(yplot,xplot,[.01 .01]*5,[.01 .01]*5,[.01 .01]*5);
counter=1;
treat=3;
for dpf=3:7
    counter=1+(dpf-3)*7;
    for fish_nb=1:7
        ZS_temp=DF_guts{treat,dpf-2,fish_nb};
        idx_temp=K_clusters_5{treat,dpf-2,fish_nb,1};
        ZS_mean=zeros(5,1200);
        for ij=1:5
            if idx_temp
                temp=nanmean(ZS_temp(idx_temp==ij,:),1);
                if temp
                    ZS_mean(ij,:)=temp;
                end
            end
        end
        axes(ha(counter));
        imagesc(ZS_mean,[0 1]);colormap hot
        set(gca,'XTickLabel',[]);set(gca,'YTickLabel',[]);
        title(strcat(num2str(size(DF_guts{treat,dpf-2,fish_nb},1)),' fish nb ',num2str(fish_nb),' DPF ',num2str(dpf)));
        counter=counter+1;
    end
end