Baseline=cell2mat(Baseline);
Baseline=Baseline(idx_components+1);
Traces=DenoisedTraces(idx_components+1,:);
Traces=Traces(Good_rois,:)-repmat(Baseline(Good_rois)',1,size(Traces,2));
X=Traces;
L=40;
tic
Ws = {};
Hs = {};
numfits = 5; %number of fits to compare
for k = 4:20
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
h2 = plot([3,3],[0,0.5],'r--');
legend([h1 h2], {'median Diss','true K'})
xlabel('K')
ylabel('Diss')

%% Procedure for choosing lambda
nLambdas = 30; % increase if you're patient
K = 2; 
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


%% choose lambda=.005; run multiple times, see number of sig factors
loadings = [];
pvals = []; 
is_significant = []; 
nIter = 20; % increase if patient
display('Running seqNMF multiple times for lambda=0.005')

[TrainData,idx] = datasample(Traces,ceil(size(Traces,1)*0.7),1);
idx_all=[1:1:size(Traces,1)];
testNEURAL = Traces(idx_all(~ismember(idx_all,idx)),:);
lambda = 0.001;
for iteri = 1:nIter
    [W, H, ~,loadings(iteri,:),power]= seqNMF(X,'K',K,'L',L,...
            'lambdaL1W', .1, 'lambda', lambda, 'maxiter', 100, 'showPlot', 0); 
    %p = .05;
    %[pvals(iteri,:),is_significant(iteri,:)] = test_significance(testNEURAL,W,p);
    %W = W(:,is_significant(iteri,:)==1,:); 
    %H = H(is_significant(iteri,:)==1,:); 
    [max_factor, L_sort, max_sort, hybrid] = helper.ClusterByFactor(W(:,:,:),1);
    indSort = hybrid(:,3);
    tstart = 300; 
    clf; WHPlot(W(indSort,1,:),H(1,tstart:end), X(indSort,tstart:end), 0)
    display(['seqNMF run ' num2str(iteri) '/' num2str(nIter)])
end
figure; hold on
h = histogram(sum(is_significant,2), 'edgecolor', 'w', 'facecolor', .7*[1 1 1]); 
%h.BinCounts = h.BinCounts/sum(h.BinCounts)*100; 
xlim([0 10]); 
xlabel('# significant factors')
ylabel('% seqNMF runs')

K=3;
lambda=1.0000e-04;

L=8;
% run seqNMF with lambdaOrthoH -> events based
lambdaOrthoH = .2; % favor events-based (these can take any value, don't need to be zero and one)
lambdaOrthoW = 0;
display('Running seqNMF on simulated data, lambdaOrthoH -> events based')
%figure(3); 
[W,H] = seqNMF(X,'K',K, 'L', L,'lambda', lambda, ...
    'lambdaOrthoH', lambdaOrthoH, 'lambdaOrthoW', lambdaOrthoW, 'showPlot', 0);

% sort neurons and plot
[max_factor, L_sort, max_sort, hybrid] = helper.ClusterByFactor(W(:,:,:),1);
indSort = hybrid(:,3);
tstart = 1; % plot data starting at this timebin
figure(1);
WHPlot(W(indSort,:,:),H(:,tstart:end), X(indSort,tstart:end), 1); title('lambdaOrthoH -> events based')

% run seqNMF with lambdaOrthoW -> parts based
%figure(4); 
lambdaOrthoH = 0;  
lambdaOrthoW = 1; % favor parts-based (these can take any value, don't need to be zero and one)
display('Running seqNMF on simulated data, lambdaOrthoW -> parts based')
[W_parts,H_parts] = seqNMF(X,'K',K, 'L', L,'lambda', lambda, ...
    'lambdaOrthoH', lambdaOrthoH, 'lambdaOrthoW', lambdaOrthoW, 'showPlot', 0);

% sort neurons and plot
[max_factor_parts, L_sort_parts, max_sort_parts, hybrid_parts] = helper.ClusterByFactor(W_parts(:,:,:),1);
indSort_parts = hybrid_parts(:,3);
figure(2);
WHPlot(W_parts(indSort_parts,:,:),H_parts(:,:), X(indSort_parts,:), 1); title('lambdaOrthoW -> parts based')