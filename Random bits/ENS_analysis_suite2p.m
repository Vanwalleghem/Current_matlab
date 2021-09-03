ZS=zscore(F(iscell(:,1)==1,:),1,2);
framerate=2;
x = linspace(0.2,size(ZS,2)/framerate,size(ZS,2));
y = linspace(1,size(ZS,1)/100,size(ZS,1));

CorrMatrix_dist=pdist(ZS,'correlation');
CorrMatrix=squareform(CorrMatrix_dist);

[RdBu]=cbrewer('div','RdBu',101);
[BuGn]=cbrewer('seq','BuGn',51);

Fighandle=figure;
set(Fighandle, 'Position', [100, 100, 1000, 1000]);set(gca,'xtick',[]);set(gca,'ytick',[]);set(gca,'xticklabel',[]);set(gca,'yticklabel',[]);
imagesc(1-CorrMatrix,[-1 1]);colormap(RdBu);

Fighandle=figure;
set(Fighandle, 'Position', [100, 100, 1000, 500]);
imagesc(x,y,ZS,[-0.5 3]);colormap(hot);%colormap(flipud(BuGn));

ROI_centroid=zeros(sum(iscell(:,1)),2);
idx_temp=find(iscell(:,1)==1);
for ij=1:sum(iscell(:,1))
   ROI_centroid(ij,:)=stat{idx_temp(ij)}.med;
end
clearvars ij idx_temp
figure;scatter(ROI_centroid(:,2),ROI_centroid(:,1));axis([0 size(ops.meanImg,1) 0 size(ops.meanImg,2)]);

Correlation_ROIs=pdist(ZS,'correlation');

edges=[0:0.05:2];
Fighandle=figure;
%set(Fighandle, 'Position', [10,10, 800, 200]);
histogram(Correlation_ROIs,edges,'normalization','probability','EdgeAlpha',0.5,'EdgeColor','k','FaceColor','g');hold on;
ax = gca;
outerpos = ax.OuterPosition;
ti = ax.TightInset; 
left = outerpos(1) + ti(1);
bottom = outerpos(2) + ti(2);
ax_width = outerpos(3) - ti(1) - ti(3);
ax_height = outerpos(4) - ti(2) - ti(4);
ax.Position = [left bottom ax_width ax_height];

Corr_temp=(1-squareform(Correlation_ROIs));
Corr_temp=weight_conversion(Corr_temp,'autofix');

Bu=cbrewer('div','RdBu',500);
Fighandle=figure;
set(Fighandle, 'Position', [10,10, 500, 500]);
imagesc(Corr_temp,[-1 1]);colormap(Bu)
ax = gca;
outerpos = ax.OuterPosition;
ti = ax.TightInset; 
left = outerpos(1) + ti(1);
bottom = outerpos(2) + ti(2);
ax_width = outerpos(3) - ti(1) - ti(3);
ax_height = outerpos(4) - ti(2) - ti(4);
ax.Position = [left bottom ax_width ax_height];

Correlation_matrix=threshold_absolute(Corr_temp,0.75);
Fighandle=figure;
graph_temp=graph(Correlation_matrix,'omitselfloops');
set(Fighandle, 'Position', [10,10, 1000, 1000]);
imagesc(ops.meanImg);hold on;
%scatter(ROI_centroid(~test,2),ROI_centroid(~test,1),20,'g','filled');hold on;
LWidths = 1*graph_temp.Edges.Weight/max(graph_temp.Edges.Weight);
LWidths(~isfinite(LWidths))=0.01;
LWidths(LWidths==0)=0.01;
LColors = ones(length(LWidths),3);LColors=LColors-LWidths;
plot(graph_temp,'EdgeColor','k','LineStyle','-','NodeColor','m','LineWidth',LWidths,'MarkerSize',5,'Marker','o','MarkerSize',2,'XData',ROI_centroid(:,2),'YData',ROI_centroid(:,1),'NodeLabel',char.empty(100,0));hold off;
ax = gca;
outerpos = ax.OuterPosition;
ti = ax.TightInset; 
left = outerpos(1) + ti(1);
bottom = outerpos(2) + ti(2);
ax_width = outerpos(3) - ti(1) - ti(3);
ax_height = outerpos(4) - ti(2) - ti(4);
ax.Position = [left bottom ax_width ax_height];
set(gca,'Visible','off')

X=F(iscell(:,1)==1,:);
L=40;
tic
Ws = {};
Hs = {};
numfits = 5; %number of fits to compare
for k = 1:20
    display(sprintf('running seqNMF with K = %i',k))
    for ii = 1:numfits
        [Ws{ii,k},Hs{ii,k}] = seqNMF(X,'K',k, 'L', L,'lambda', 0,'maxiter',30,'showplot',0); 
        % note that max iter set low (30iter) for speed in demo (not recommended in practice)
    end
    inds = nchoosek(1:numfits,2);
    for i = 1:size(inds,1) % consider using parfor for larger numfits
            Diss(i,k) = helper.DISSX(Hs{inds(i,1),k},Ws{inds(i,1),k},Hs{inds(i,2),k},Ws{inds(i,2),k});
    end
    
end

%% Plot Diss and choose K with the minimum average diss.
figure;
plot(1:20,Diss,'ko'), hold on
h1 = plot(1:20,median(Diss,1),'k-','linewidth',2);
%legend([h1 h2], {'median Diss','true K'})
xlabel('K')
ylabel('Diss')

%% Procedure for choosing lambda
nLambdas = 50; % increase if you're patient
K = 4; 
lambdas = sort([logspace(-1,-5,nLambdas)], 'ascend'); 
loadings = [];
regularization = []; 
cost = []; 
for li = 1:length(lambdas)
    [N,T] = size(X);
    [W, H, ~,loadings(li,:),power]= seqNMF(X,'K',K,'L',L,...
        'lambdaL1W', .1, 'lambda', lambdas(li), 'maxiter', 100, 'showPlot', 0); 
    [cost(li),regularization(li),~] = helper.get_seqNMF_cost(X,W,H);
    display(['Testing lambda ' num2str(li) '/' num2str(length(lambdas))])
end
%% plot costs as a function of lambda
figure;
windowSize = 3; 
b = (1/windowSize)*ones(1,windowSize);
a = 1;
Rs = filtfilt(b,a,regularization); 
minRs = prctile(regularization,10); maxRs= prctile(regularization,90);
Rs = (Rs-minRs)/(maxRs-minRs); 
R = (regularization-minRs)/(maxRs-minRs); 
Cs = filtfilt(b,a,cost); 
minCs =  prctile(cost,10); maxCs =  prctile(cost,90); 
Cs = (Cs -minCs)/(maxCs-minCs); 
C = (cost -minCs)/(maxCs-minCs); 

clf; hold on
plot(lambdas,Rs, 'b')
plot(lambdas,Cs,'r')
scatter(lambdas, R, 'b', 'markerfacecolor', 'flat');
scatter(lambdas, C, 'r', 'markerfacecolor', 'flat');
xlabel('Lambda'); ylabel('Cost (au)')
set(legend('Correlation cost', 'Reconstruction cost'), 'Box', 'on')
set(gca, 'xscale', 'log', 'ytick', [], 'color', 'none')
set(gca,'color','none','tickdir','out','ticklength', [0.025, 0.025])

%%
K=4;
L=10;
lambda=0.0001;
% run seqNMF with lambdaOrthoH -> events based
lambdaOrthoH = 0.2; % favor events-based (these can take any value, don't need to be zero and one)
lambdaOrthoW = 0;
display('Running seqNMF on simulated data, lambdaOrthoH -> events based')

X2=X-repmat(min(X,[],2),1,size(X,2));
%figure(3); 
[W,H] = seqNMF(X2,'K',K, 'L', L,'lambda', lambda, ...
    'lambdaOrthoH', lambdaOrthoH, 'lambdaOrthoW', lambdaOrthoW,'maxiter', 20, 'showPlot', 0);

% sort neurons and plot
[max_factor, L_sort, max_sort, hybrid] = helper.ClusterByFactor(W(:,:,:),1);
indSort = hybrid(:,3);
tstart = 1; % plot data starting at this timebin
figure;
WHPlot(W(indSort,:,:),H(:,tstart:end), X(indSort,tstart:end), 1); title('lambdaOrthoH -> events based')

% run seqNMF with lambdaOrthoW -> parts based
%figure(4); 
lambdaOrthoH = 0;  
lambdaOrthoW = 1; % favor parts-based (these can take any value, don't need to be zero and one)
display('Running seqNMF on simulated data, lambdaOrthoW -> parts based')
[W_parts,H_parts] = seqNMF(X2,'K',K, 'L', L,'lambda', lambda, ...
    'lambdaOrthoH', lambdaOrthoH, 'lambdaOrthoW', lambdaOrthoW, 'showPlot', 0);

% sort neurons and plot
[max_factor_parts, L_sort_parts, max_sort_parts, hybrid_parts] = helper.ClusterByFactor(W_parts(:,:,:),1);
indSort_parts = hybrid_parts(:,3);
figure;
WHPlot(W_parts(indSort_parts,:,:),H_parts(:,:), X(indSort_parts,:), 1); title('lambdaOrthoW -> parts based')
