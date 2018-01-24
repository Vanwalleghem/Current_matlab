GMModels = {};
options = statset('Display','final','MaxIter',500,'TolFun',1e-6);
for k = 1:20
    GMModels{k} = fitgmdist(ZS2_rsq,k,'Options',options,'CovarianceType','diagonal','Options',options, 'Regularize', 1e-5);
    BIC(k)= GMModels{k}.BIC;
end
figure;plot(BIC);
AIC=[];
for i=2:length(GMModels)
    AIC(i)=GMModels{i}.AIC;
end
hold on;plot(AIC);

figure;
counter=1;xplot=floor(sqrt(18));yplot=ceil(18/xplot);
for i=1:18
    subplot(xplot,yplot,i);plot(mean(ZS2_rsq(find(idxKmeans_GM==i),:),1));
end

idxKmeans_GM=cluster(GMModels{18},ZS2_rsq);

BIC=[];% Bayesian Information Criterion
n=3100; % number of datapoints
temp=1;
for k=1:50  % number of clusters
    RSS=0;  % residual sum of squares
    [idx,C]=kmeans(y,k);  % Matlab command for k-mean clustering
    for i=1:3100
        RSS=RSS+sqrt((y(i,1)-C(idx(i),1))^2+(y(i,2)-C(idx(i),2))^2);
    end
    BIC(temp)=n*log(RSS/n)+(k*3)*log(n);
    temp=temp+1;
end
[p,l]=min(BIC);
plot(BIC)

figure;k=10;
idxKmeans_GM=cluster(GMModels{k},ZS2_rsq);
counter=1;xplot=floor(sqrt(k));yplot=ceil(k/xplot);
for i=1:k
    subplot(xplot,yplot,i);plot(mean(ZS2_rsq(find(idxKmeans_GM==i),:),1));
end


GMModels_AVG = {};
options = statset('Display','final','MaxIter',500,'TolFun',1e-6);
for k = 21:30
    GMModels_AVG{k} = fitgmdist(ZS_AVG2(idx_rsq,:),k,'Options',options,'CovarianceType','diagonal','Options',options, 'Regularize', 1e-5);
    BIC(k)= GMModels_AVG{k}.BIC;
end
figure;plot(BIC);
AIC=[];
for i=2:length(GMModels_AVG)
    AIC(i)=GMModels_AVG{i}.AIC;
end
hold on;plot(AIC);

figure;k=5;
ZS_AVG2_rsq=ZS_AVG2(idx_rsq,:);
idxKmeans_GMAVG=cluster(GMModels_AVG{k},ZS_AVG2_rsq);
counter=1;xplot=floor(sqrt(k));yplot=ceil(k/xplot);
for i=1:k
subplot(xplot,yplot,i);plot(mean(ZS_AVG2_rsq(find(idxKmeans_GMAVG==i),:),1));
end
[Model_ZS_GMMAVG,GoodBetas_ZS_GMMAVG]=Test_Regress(GMModels_AVG{k}.mu,Stimuli_AVG,idxKmeans_GMAVG,0.6);


D = pdist(ZS2_rsq,'cityblock');
Z=linkage(D,'single');
T = cluster(Z,'maxclust',30);

Z2 = linkage(ZS2_rsq,'single','correlation') 
Y2 = inconsistent(Z2);figure;histogram(Y2);
T2 = cluster(Z2,'maxclust',30);


clusters=[];
for i=1:max(T)
    clusters(i,:)=mean(ZS2_rsq(find(T==i),:),1);
end
figure;imagesc(clusters);

clusters2=[];
for i=1:max(T2)
    clusters2(i,:)=mean(ZS2_rsq(find(T2==i),:),1);
end
figure;imagesc(clusters2);
