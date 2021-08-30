Dir_list=dir('suite2p*');
Dir_list=Dir_list([Dir_list.isdir]==1);

Suite2p_traces={};
S2p_concTraces=[];
ROIs={};threshold=0.5;
progressbar;
for i=1:length(Dir_list)
    progressbar(i/length(Dir_list));
    load(strcat(Dir_list(i).name,'\plane0\Fall.mat'));
    idx_cell=find(iscell(:,2)>threshold);
    Suite2p_traces{i}=F(idx_cell,:);
    temp=Suite2p_traces{i};
    if size(temp,2)<1200
        temp=padarray(temp,[0 1200-size(temp,2)],'replicate','post');
    end
    S2p_concTraces=vertcat(S2p_concTraces,temp);
    %     Below takes too much memory and no spars ND array in matlab
    ROIs_Suite2p=zeros(ops.Ly,ops.Lx,length(idx_cell));
    ROIs_Suite2p=ndSparse(ROIs_Suite2p,size(ROIs_Suite2p));
    for ij=1:length(idx_cell)
        temp=[stat{idx_cell(ij)}.ypix;stat{idx_cell(ij)}.xpix];
        for ik=1:length(temp)
            ROIs_Suite2p(temp(1,ik)+1,temp(2,ik)+1,ij)=stat{idx_cell(ij)}.lam(ik);
        end
    end
    ROIs{i}=ROIs_Suite2p;
end
close all;
clearvars i ij ik ROIs_Suite2p temp idx_cell stat spks F Fneu iscell ops

i=1;ROIsCont=ROIs{i};
for i=2:length(ROIs)
    ROIsCont=cat(3,ROIsCont,ROIs{i});
end

DF=DeltaF2(S2p_concTraces,50,11);
DF(~isfinite(DF))=0;

%Design a regular expression that match your naming scheme (https://regexr.com/)
MatFiles_names={Dir_list.name};
fin = cellfun(@(x)regexp(x,'fish(\d+)_ENS(.*)_(\d)DPF','tokens'), MatFiles_names, 'UniformOutput', false);
names=zeros(3,length(fin));
for i=1:length(fin)
    %Fish=str2num(strcat(fin{i}{1}{1},fin{i}{1}{2}));%Concatenate all the matched groups from the regex
    names(1,i)=str2num(fin{i}{1}{1});names(2,i)=str2num(fin{i}{1}{3});
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

%You can also build an index of the Fish and the plane of each ROI
Numbers= [0 cumsum(cellfun('size',Suite2p_traces,1))];
idx_DPF=nan(size(DF,1),1);
idx_Fish=nan(size(DF,1),1);
idx_Treatment=nan(size(DF,1),1);
for i=1:length(Dir_list)
    idx_DPF(Numbers(i)+1:Numbers(i+1))=names(2,i);
    idx_Fish(Numbers(i)+1:Numbers(i+1))=names(1,i);
    idx_Treatment(Numbers(i)+1:Numbers(i+1))=names(3,i);
end
clearvars i Fish Plane name counter

ROIs_med=zeros(size(S2p_concTraces,1),2);
counter=1;
for i=1:length(Dir_list)
    load(strcat(Dir_list(i).name,'\plane0\Fall.mat'),'stat');
    load(strcat(Dir_list(i).name,'\plane0\Fall.mat'),'iscell');
    idx_cell=find(iscell(:,2)>0.5);
    for ij=1:length(idx_cell)
        ROIs_med(counter,:)=stat{idx_cell(ij)}.med;
        counter=counter+1;
    end
end
%% Get the limits of the gut
% img_mean={};
% Fighandle=figure;
% set(Fighandle, 'Position', [100, 100, 640, 480]);
% for i=1:length(Dir_list)
%     load(strcat(Dir_list(i).name,'\plane0\Fall.mat'),'ops');
%     imagesc(ops.meanImg);title(num2str(i));
%     img_mean{i}=ops.meanImg;
%     pause;
% end


%% Read the gut limits

fid = fopen('Gut limits_suite2p.txt');
Limits=inf(2,length(Dir_list));
Limits(2,:)=0;
counter=1;
while ~feof(fid) % feof(fid) is true when the file ends
    textLineEntry = fgetl(fid); % read one line
    if regexp(textLineEntry,'x')
        Limits(1,counter)=str2num(textLineEntry(3:end));
    elseif regexp(textLineEntry,'y')
        Limits(2,counter)=str2num(textLineEntry(3:end));
    else
        warning('error');
    end
    counter=counter+1;
end
fclose(fid); % close the file

idx_Gut=nan(size(DF,1),1);
for i=1:length(Dir_list)
    ROI_temp=ROIs_med(Numbers(i)+1:Numbers(i+1),:);
    %idx_ROI=(ROI_temp(:,1)<Limits(2,i) & ROI_temp(:,1)>10) & (ROI_temp(:,2)<Limits(1,i) & ROI_temp(:,2)>10);
    idx_ROI=(ROI_temp(:,1)>Limits(2,i) & ROI_temp(:,2)<Limits(1,i));
    idx_Gut(Numbers(i)+1:Numbers(i+1))=idx_ROI;
end
clearvars i Fish Plane name counter
idx_Gut=boolean(idx_Gut);


ZS_guts={};
for treat=1:3
    for dpf=3:7
        idx_temp=find(idx_DPF==dpf & idx_Treatment==treat);
        idx_temp=idx_temp(idx_Gut(idx_temp));
        ZS_guts{treat,dpf-2,1}=DF(idx_temp,:);
        ZS_guts{treat,dpf-2,2}=idx_temp;
    end
end

%% run rastermap to smooth the activity

ops.nCall=[5 200];
ops.useGPU=1;
smooth_ZS_fish={};
% results_ENS=struct();
for treat=1:3
    for dpf=3:7
        idx_temp=ZS_guts{treat,dpf-2,2};
        for fish_nb=1:max(unique(idx_Fish(idx_temp)))
            ZS_temp=DF(idx_temp(idx_Fish(idx_temp)==fish_nb),:);
            if ~isempty(ZS_temp)
                ops.iPC=1:size(ZS_temp,1);
                [isort1, isort2, Sm] = mapTmap(ZS_temp,ops);
                smooth_ZS_fish{treat,dpf-2,fish_nb,1}=isort1;
                smooth_ZS_fish{treat,dpf-2,fish_nb,2}=isort2;
                smooth_ZS_fish{treat,dpf-2,fish_nb,3}=Sm;
                %                 pks_smooth=cell(1,size(Sm,1));
                %                 locs_smooth=cell(1,size(Sm,1));
                %                 Periodicity_smooth=nan(size(Sm,1),1);
                %                 for roi_nb=1:size(ZS_temp,1)
                %                     [pks_smooth{roi_nb},locs_smooth{roi_nb}] = findpeaks(Sm(roi_nb,:),'MinPeakProminence',1.5,'MinPeakDistance',5);
                %                     if length(locs_smooth{roi_nb})>5
                %                         locs_temp=locs_smooth{roi_nb}(2:end)-locs_smooth{roi_nb}(1:end-1);
                %                         locs_temp=locs_temp/2;
                %                         Periodicity_smooth(roi_nb)=mean(locs_temp);
                %                     end
                %                 end
                %                 results_ENS(treat,dpf-2,fish_nb).peaks=pks_smooth;
                %                 results_ENS(treat,dpf-2,fish_nb).locs=locs_smooth;
                %                 results_ENS(treat,dpf-2,fish_nb).Periodicity=Periodicity_smooth;
            end
        end
    end
end


%% compute basic metrics

results_ENS=struct();
results_ENS_raw=struct();
for treat=1:3
    for dpf=3:7
        idx_temp=ZS_guts{treat,dpf-2,2};
        for fish_nb=1:max(unique(idx_Fish(idx_temp)))
            ZS_temp=DF(idx_temp(idx_Fish(idx_temp)==fish_nb),:);
            if ~isempty(ZS_temp)
                Sm=smooth_ZS_fish{treat,dpf-2,fish_nb,3};
                pks_smooth=cell(1,size(Sm,1));
                locs_smooth=cell(1,size(Sm,1));
                Periodicity_smooth=nan(size(Sm,1),1);
                pks_raw=cell(1,size(Sm,1));
                locs_raw=cell(1,size(Sm,1));
                Periodicity_raw=nan(size(Sm,1),1);
                for roi_nb=1:size(ZS_temp,1)
                    [pks_smooth{roi_nb},locs_smooth{roi_nb}] = findpeaks(Sm(roi_nb,:),'MinPeakProminence',0.2,'MinPeakDistance',5);
                    [pks_raw{roi_nb},locs_raw{roi_nb}] = findpeaks(ZS_temp(roi_nb,:),'MinPeakProminence',0.2,'MinPeakDistance',5);
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

%% get ready for Prism and plotting

PrismTemp=nan(3,5*8);
for treat=1:3
    for dpf=3:7
        for fish_nb=1:8
            PrismTemp(treat,fish_nb+(dpf-3)*8)=nanmean(results_ENS(treat,dpf-2,fish_nb).Periodicity);
        end
    end
end
PrismTemp(PrismTemp==0)=nan;

edge=[0:1:20];
treat=3;
h=nan(length(edge)-1,5*8);
for dpf=3:7
    for fish_nb=1:8
        temp=results_ENS(treat,dpf-2,fish_nb).Periodicity;
        temp(isnan(temp))=[];
        h(:,fish_nb+(dpf-3)*8)=histcounts(temp,edge,'Normalization','probability')';
    end
end


Fighandle=figure;
set(Fighandle, 'Position', [100, 100, 1200, 400]);
set(findall(Fighandle,'type','text'),'fontSize',12,'fontWeight','bold','FontName','Arial')
xplot=8;yplot=5;
ha = tight_subplot(yplot,xplot,[.01 .01]*5,[.01 .01]*5,[.01 .01]*5);
counter=1;
treat=1;
for dpf=3:7
    counter=1+(dpf-3)*8;
    idx_temp=ZS_guts{treat,dpf-2,2};
    for fish_nb=1:max(unique(idx_Fish(idx_temp)))
        axes(ha(counter));
        imagesc(smooth_ZS_fish{treat,dpf-2,fish_nb,3});
        set(gca,'XTickLabel',[]);set(gca,'YTickLabel',[]);
        title(num2str(sum(idx_Fish(idx_temp)==fish_nb)));
        counter=counter+1;
    end
end

Fighandle=figure;
set(Fighandle, 'Position', [100, 100, 1200, 400]);
set(findall(Fighandle,'type','text'),'fontSize',12,'fontWeight','bold','FontName','Arial')
xplot=8;yplot=5;
ha = tight_subplot(yplot,xplot,[.01 .01],[.01 .01],[.01 .01]);
counter=1;
treat=1;
for dpf=3:7
    counter=1+(dpf-3)*8;
    idx_temp=ZS_guts{treat,dpf-2,2};
    for fish_nb=1:max(unique(idx_Fish(idx_temp)))
        axes(ha(counter));
        ZS_temp=DF(idx_temp(idx_Fish(idx_temp)==fish_nb),:);
        imagesc(ZS_temp,[0 0.1]);
        set(gca,'XTickLabel',[]);set(gca,'YTickLabel',[]);
        counter=counter+1;
    end
end

Fighandle=figure;
set(Fighandle, 'Position', [100, 100, 1200, 400]);
set(findall(Fighandle,'type','text'),'fontSize',12,'fontWeight','bold','FontName','Arial')
xplot=8;yplot=5;
ha = tight_subplot(yplot,xplot,[.01 .01],[.01 .01],[.01 .01]);
counter=1;
treat=1;
for dpf=3:7
    counter=1+(dpf-3)*8;
    idx_temp=ZS_guts{treat,dpf-2,2};
    for fish_nb=1:max(unique(idx_Fish(idx_temp)))
        axes(ha(counter));
        ZS_temp=full(squeeze(sum(ROIsCont(:,:,idx_temp(idx_Fish(idx_temp)==fish_nb)),3)));
        imagesc(ZS_temp);
        set(gca,'XTickLabel',[]);set(gca,'YTickLabel',[]);
        counter=counter+1;
    end
end

%% raw prism data
PrismTemp=nan(3,5*8);
for treat=1:3
    for dpf=3:7
        for fish_nb=1:8
            PrismTemp(treat,fish_nb+(dpf-3)*8)=nanmean(results_ENS_raw(treat,dpf-2,fish_nb).Periodicity);
        end
    end
end
PrismTemp(PrismTemp==0)=nan;

edge=[0:2:60];
treat=3;
h=nan(length(edge)-1,5*8);
for dpf=3:7
    for fish_nb=1:8
        temp=results_ENS_raw(treat,dpf-2,fish_nb).Periodicity;
        temp(isnan(temp))=[];
        h(:,fish_nb+(dpf-3)*8)=histcounts(temp,edge,'Normalization','probability')';
    end
end

%% Now correlation/distance
Correlation_ENS=struct();
for treat=1:3
    for dpf=3:7
        idx_temp=ZS_guts{treat,dpf-2,2};
        for fish_nb=1:max(unique(idx_Fish(idx_temp)))
            ZS_temp=DF(idx_temp(idx_Fish(idx_temp)==fish_nb),:);
            ROIs_temp=ROIs_med(idx_temp(idx_Fish(idx_temp)==fish_nb));
            if ~isempty(ZS_temp)
                Correlation_ENS(treat,dpf,fish_nb).correlation=pdist(ZS_temp,'correlation');
                Correlation_ENS(treat,dpf,fish_nb).distance=pdist(ROIs_temp);
            end
        end
    end
end

edges={[0:0.05:2],[0:0.05:2]};
Fighandle=figure;
set(Fighandle, 'Position', [100, 100, 1200, 400]);
set(findall(Fighandle,'type','text'),'fontSize',12,'fontWeight','bold','FontName','Arial')
xplot=8;yplot=5;
ha = tight_subplot(yplot,xplot,[.01 .01],[.01 .01],[.01 .01]);
counter=1;
treat=1;
for dpf=3:7
    counter=1+(dpf-3)*8;
    idx_temp=ZS_guts{treat,dpf-2,2};
    for fish_nb=1:max(unique(idx_Fish(idx_temp)))
        if sum(idx_temp(idx_Fish(idx_temp)==fish_nb))
            Temp=hist3([Correlation_ENS(treat,dpf,fish_nb).correlation' (Correlation_ENS(treat,dpf,fish_nb).distance/max(Correlation_ENS(treat,dpf,fish_nb).distance/2))'],'Edges',edges);
            axes(ha(counter));
            imagesc(Temp);
            set(gca,'XTickLabel',[]);set(gca,'YTickLabel',[]);
        end
        counter=counter+1;
    end
end

%Mean of corr vs dist across fish
edges={[0:0.05:2],[0:0.05:2]};
Fighandle=figure;
set(Fighandle, 'Position', [100, 100, 1200, 400]);
set(findall(Fighandle,'type','text'),'fontSize',12,'fontWeight','bold','FontName','Arial')
xplot=5;yplot=3;
ha = tight_subplot(yplot,xplot,[.01 .01],[.01 .01],[.01 .01]);
counter=1;
for treat=1:3
    for dpf=3:7
        idx_temp=ZS_guts{treat,dpf-2,2};
        TempCat=[];
        for fish_nb=1:max(unique(idx_Fish(idx_temp)))
            if sum(idx_temp(idx_Fish(idx_temp)==fish_nb))
                Temp=hist3([Correlation_ENS(treat,dpf,fish_nb).correlation' (Correlation_ENS(treat,dpf,fish_nb).distance/max(Correlation_ENS(treat,dpf,fish_nb).distance/2))'],'Edges',edges);
            end
            TempCat=cat(3,TempCat,Temp);
        end
        axes(ha(counter));
        imagesc(squeeze(mean(TempCat,3)));colormap magma;
        set(gca,'XTickLabel',[]);set(gca,'YTickLabel',[]);
        counter=counter+1;
    end
end

PrismTemp=nan(3,5*8);
for treat=1:3
    for dpf=3:7
        for fish_nb=1:max(unique(idx_Fish(idx_temp)))
            Temp=Correlation_ENS(treat,dpf,fish_nb).correlation;
            PrismTemp(treat,fish_nb+(dpf-3)*8)=nanmean(1-Temp);
        end
    end
end

%Distrib
edge=[-1:0.1:1];
h=nan(length(edge)-1,5*8);
treat=3;
for dpf=3:7
    idx_temp=ZS_guts{treat,dpf-2,2};
    for fish_nb=1:max(unique(idx_Fish(idx_temp)))
        ZS_temp=DF(idx_temp(idx_Fish(idx_temp)==fish_nb),:);
        if ~isempty(ZS_temp)
            h(:,fish_nb+(dpf-3)*8)=histcounts(1-Correlation_ENS(treat,dpf,fish_nb).correlation,edge,'Normalization','probability')';
        end
    end
end


%% network analysis
Graph_ENS=struct();Threshold=0.2;%1-correlation
for treat=1:3
    for dpf=3:7
        for fish_nb=1:max(unique(idx_Fish))
            temp=Correlation_ENS(treat,dpf,fish_nb).correlation;temp(isnan(temp))=1;
            if ~isempty(temp)
                Graph_ENS(treat,dpf,fish_nb).density=density_und(squareform(temp<Threshold));
                Graph_ENS(treat,dpf,fish_nb).degree=degrees_und(squareform(temp<Threshold));
            end
        end
    end
end

PrismTemp=nan(3,5*8);
for treat=1:3
    for dpf=3:7
        for fish_nb=1:max(unique(idx_Fish))
            Temp=Graph_ENS(treat,dpf,fish_nb).density;
            PrismTemp(treat,fish_nb+(dpf-3)*8)=nanmean(Temp);
        end
    end
end

PrismTemp=nan(3,5*8);
for treat=1:3
    for dpf=3:7
        for fish_nb=1:max(unique(idx_Fish))
            Temp=Graph_ENS(treat,dpf,fish_nb).degree;
            PrismTemp(treat,fish_nb+(dpf-3)*8)=nanmean(Temp);
        end
    end
end

%% peak analysis

edge=[0:1:20];
h=nan(length(edge)-1,5*8);
treat=1;
for dpf=3:7
    idx_temp=ZS_guts{treat,dpf-2,2};
    for fish_nb=1:max(unique(idx_Fish(idx_temp)))
        ZS_temp=DF(idx_temp(idx_Fish(idx_temp)==fish_nb),:);
        if ~isempty(ZS_temp)
            peak_temp_nb=[];
            peak_temp=results_ENS_raw(treat,dpf-2,fish_nb).peaks;
            for i=1:length(peak_temp)
                peak_temp_nb(i)=length(peak_temp{i});
            end
            h(:,fish_nb+(dpf-3)*8)=histcounts(peak_temp_nb,edge,'Normalization','probability')';
        end
    end
end

%lots of inactive neurons
edge=[0:2:60];
h=nan(length(edge)-1,5*8);
treat=3;
for dpf=3:7
    idx_temp=ZS_guts{treat,dpf-2,2};
    for fish_nb=1:max(unique(idx_Fish(idx_temp)))
        ZS_temp=DF(idx_temp(idx_Fish(idx_temp)==fish_nb),:);
        if ~isempty(ZS_temp)
            peak_temp_nb=[];
            peak_temp=results_ENS_raw(treat,dpf-2,fish_nb).peaks;
            for i=1:length(peak_temp)
                peak_temp_nb(i)=length(peak_temp{i});
            end
            peak_temp_nb(peak_temp_nb==0)=[];
            h(:,fish_nb+(dpf-3)*8)=histcounts(peak_temp_nb,edge,'Normalization','probability')';
        end
    end
end


PrismTemp=nan(3,5*8);
for treat=1:3
    for dpf=3:7
        for fish_nb=1:max(unique(idx_Fish(idx_temp)))
            if ~isempty(ZS_temp)
                peak_temp_nb=[];
                peak_temp=results_ENS_raw(treat,dpf-2,fish_nb).peaks;
                for i=1:length(peak_temp)
                    peak_temp_nb(i)=length(peak_temp{i});
                end
            end
            PrismTemp(treat,fish_nb+(dpf-3)*8)=nanmean(peak_temp_nb(peak_temp_nb>0));
        end
    end
end




%% Cross-correlation matrices
Cross_corr_temp=nan(size(ZS_temp,1));
Cross_corr_lag=nan(size(ZS_temp,1));
for rows=1:size(ZS_temp,1)
    for columns=1:size(ZS_temp,1)
        if columns>=rows
            [Cross_corr_temp(rows,columns) Cross_corr_lag(rows,columns)]=max(xcorr(ZS_temp(rows,:)',ZS_temp(columns,:)',10,'coeff'));
        end
    end
end
Cross_corr_temp(isnan(Cross_corr_temp))=0;
Cross_corr_temp(1:1+size(Cross_corr_temp,1):end) = 1;
temp=Cross_corr_temp+tril(Cross_corr_temp',-1);


CrossCorrelation_ENS=struct();
for treat=1:3
    for dpf=3:7
        idx_temp=ZS_guts{treat,dpf-2,2};
        for fish_nb=1:max(unique(idx_Fish(idx_temp)))
            ZS_temp=DF(idx_temp(idx_Fish(idx_temp)==fish_nb),:);
            ROIs_temp=ROIs_med(idx_temp(idx_Fish(idx_temp)==fish_nb));
            if ~isempty(ZS_temp)
                Cross_corr_temp=nan(size(ZS_temp,1));
                Cross_corr_lag=nan(size(ZS_temp,1));
                for rows=1:size(ZS_temp,1)
                    parfor columns=1:size(ZS_temp,1)
                        if columns>=rows
                            [Cross_corr_temp(rows,columns) Cross_corr_lag(rows,columns)]=max(xcorr(ZS_temp(rows,:)',ZS_temp(columns,:)',10,'coeff'));
                        end
                    end
                end
                CrossCorrelation_ENS(treat,dpf,fish_nb).correlation=Cross_corr_temp;
                CrossCorrelation_ENS(treat,dpf,fish_nb).lag=Cross_corr_lag;
            end
        end
    end
end
