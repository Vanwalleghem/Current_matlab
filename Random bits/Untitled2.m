
for fish_nb=1:Fish_list
    Clustering(fish_nb).name=Fish_list(fish_nb);
    idx_temp=find(idx_Fish==Fish_list(fish_nb));
    idxKmeans_final=zeros(size(ZS,1),1);
    idxKmeans_final(idx_temp)=Clustering(fish_nb).idxKmeans_ZS;
    Clustering(fish_nb).idxKmeans_final=idxKmeans_final;
    Numbers=[0 [ROIs_idx]];
    temp=[];
    counter=1;
    for i=1:length(Clustering(fish_nb).GoodBetas)
        temp{counter}=find(idxKmeans_final==Clustering(fish_nb).GoodBetas(i));
        counter=counter+1;
    end
end
