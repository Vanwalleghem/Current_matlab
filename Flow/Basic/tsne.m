figure;

Y=tsne(ZS_rsq,[],2,300,100);
figure;
scatter(Y(:,1),Y(:,2),20,colorKmeans/255,'filled');

idxKmeans_ZS_goodmembers_rsq=idxKmeans_ZS_goodmembers(idx_rsq);
colorKmeans=zeros(length(idxKmeans_ZS_goodmembers_rsq),3);

for i=1:8
    for j=1:3
        colorKmeans(idxKmeans_ZS==GoodBetas_select(i),j)=colors(i,j);
    end
end

Yb=fast_tsne(ZS_rsq,2,400,50,0.1);
figure;
scatter(Yb(:,1),Yb(:,2),20,colorKmeans/255,'filled');

