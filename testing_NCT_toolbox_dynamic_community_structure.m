
%%%% now I will be testing some of the parameters of the community
%%%% structure of the multilayer network from the script
%%%% testing_GenLouvain_v2_2_toolbox.m. 

%%% i got them from: http://commdetect.weebly.com/ 
%%% made by Danielle S. Bassett and her team. 

%%% I wil be testing it for the f20 weighted graphs of habituation

%%


%%%% just checking how many nodes I have per brain region and cluster

%%% the clustID tags are not in the right order... I will change it
oldclusttag=[4 2 5 6 1 7];
ClustID=zeros(size(Nodes.Nod_clustID(keep)));
for i=1:length(unique(Nodes.Nod_clustID(keep)))
    
    idx_temp=find(Nodes.Nod_clustID(keep)==oldclusttag(i));
    ClustID(idx_temp)=i;
   
end
 
clustnames={'fasthab1','fasthab2','fasthab3','modhab','weakhab','inhib'};

figure;
subplot(1,2,1);histogram(Nodes.Nod_brainID(keep)); xticks([1:9]); xticklabels(RegionList);xtickangle(45);
subplot(1,2,2);histogram(ClustID); xticks([1:6]); xticklabels(clustnames);xtickangle(45);

%% doing it for a multilayer matrix

%%% first for the averaged matrix for f20 in time 

%%% making the cell array
A={};

for loom=1:31
   temp=Data_corrMat2.f20.Mean_corrMat{1,loom}(keep,keep);
   temp(isnan(temp))=0;
    A{loom}=temp;  
end
clear temp

%%

%%% you need to set values for gamma and omega

gamma = 1;
omega = 0.5; %%% this makes a big influence... 

N=length(A{1});
T=length(A);
B=spalloc(N*T,N*T,N*N*T+2*N*T);
twomu=0;
for s=1:T
    k=sum(A{s});
    twom=sum(k);
    twomu=twomu+twom;
    indx=[1:N]+(s-1)*N;
    B(indx,indx)=A{s}-gamma*k'*k/twom;
end
twomu=twomu+2*omega*N*(T-1);
B = B + omega*spdiags(ones(N*T,2),[-N,N],N*T,N*T);
%[S,Q] = genlouvain(B);
[S,Q,nb_it] = iterated_genlouvain(B);
Q = Q/twomu
S = reshape(S,N,T);

%% checking results

figure;imagesc(S);colormap('jet');


%% testing flexibility

flexF20=flexibility(S','temp'); %%% i need to mind the orientation of the matrix to do it properly

figure;
plot(Zbrain_brainMask2D(:,1),Zbrain_brainMask2D(:,2),'k');xlim([400 1350])
%  hold on;
%  gscatter(Nodes.Nod_coor(keep,1),Nodes.Nod_coor(keep,2),Nodes.Nod_clustID(keep,:),'rgggbm',[],21); 
hold on;
scatter(Nodes.Nod_coor(keep,1),Nodes.Nod_coor(keep,2),20,flexF20,'filled'); colormap('cool'); caxis([0 1]);colorbar;
view(-90,90);
title('flexibility');
hold off;


clustnames={'strhab1','strhab2','strhab3','modhab','weakhab','inhib'};

    
figure;
subplot(1,2,1);boxplot(flexF20,Nodes.Nod_brainID(keep)),xticklabels(RegionList);xtickangle(45);ylim([0 1]);
subplot(1,2,2);boxplot(flexF20,ClustID),xticklabels(clustnames);xtickangle(45);ylim([0 1]);


%% cohesion strength and related

options.figureFlag	= 0;
options.colormap	= 'jet';

[Cij,node_cohesion,node_disjoint,node_flexibility,strength_cohesion,commChanges,commCohesion,commDisjoint,commIndex] = calc_node_cohesion(S,options);


figure;
plot(Zbrain_brainMask2D(:,1),Zbrain_brainMask2D(:,2),'k');xlim([400 1350])
% hold on;
% gscatter(Nodes.Nod_coor(keep,1),Nodes.Nod_coor(keep,2),Nodes.Nod_clustID(keep,:),'gggbrm',[],21); 
hold on;
scatter(Nodes.Nod_coor(keep,1),Nodes.Nod_coor(keep,2),20,node_cohesion,'filled'); colormap('cool'); caxis([0 1]);colorbar;
view(-90,90);
title('cohesion');
hold off;


figure;
subplot(1,2,1);boxplot(node_cohesion,Nodes.Nod_brainID(keep)),xticklabels(RegionList);xtickangle(45);ylim([0 1]);
subplot(1,2,2);boxplot(node_cohesion,ClustID),xticklabels(clustnames);xtickangle(45);ylim([0 1]);

%%

P = promiscuity(S');  %%% i need to mind the orientation of the matrix to do it properly


figure;
plot(Zbrain_brainMask2D(:,1),Zbrain_brainMask2D(:,2),'k');xlim([400 1350])
% hold on;
% gscatter(Nodes.Nod_coor(keep,1),Nodes.Nod_coor(keep,2),Nodes.Nod_clustID(keep,:),'gggbrm',[],21); 
hold on;
scatter(Nodes.Nod_coor(keep,1),Nodes.Nod_coor(keep,2),20,P,'filled'); colormap('cool'); caxis([0 1]);colorbar;
view(-90,90);
title('promiscuity');
hold off;
  
figure;
subplot(1,2,1);boxplot(P,Nodes.Nod_brainID(keep)),xticklabels(RegionList);xtickangle(45);ylim([0 1]);
subplot(1,2,2);boxplot(P,ClustID),xticklabels(clustnames);xtickangle(45);ylim([0 1]);

%%

%%% trying to make a module allegiance matrix needed for the 'integration'
%%% and 'recruitment' functions...  representing the probability that area 
%%% i and area j were assigned to the same functional community by 
%%% time-resolved clustering methods

MA=[];
for i=1:size(S,1)
for j=1:size(S,1)
   
    for loom=1:size(S,2) 
    temp_prob(loom)=S(i,loom)==S(j,loom);   
    end
    MA(i,j)=sum(temp_prob)/31;
end   
end

figure;imagesc(MA);colormap('parula');


%% Recruitment

%%% not working... I am not sure what kind of structure the systemByNode
%%% needs to be... 

% R = recruitment(MA,S);


%% persistence

%%%% is not working... I think I need a different kind of structure S... 

% pers = persistence(S);  %%% i need to mind the orientation of the matrix to do it properly
% 
% 
% figure;
% plot(Zbrain_brainMask2D(:,1),Zbrain_brainMask2D(:,2),'k');xlim([400 1350])
% % hold on;
% % gscatter(Nodes.Nod_coor(keep,1),Nodes.Nod_coor(keep,2),Nodes.Nod_clustID(keep,:),'gggbrm',[],21); 
% hold on;
% scatter(Nodes.Nod_coor(keep,1),Nodes.Nod_coor(keep,2),20,pers,'filled'); colormap('cool'); caxis([0 1]);colorbar;
% view(-90,90);
% title('persistence');
% hold off;
% 
% clustnames={'fasthab1','fasthab2','fasthab3','modhab','weakhab','inhib'};
%     
% figure;
% subplot(1,2,1);boxplot(pers,Nodes.Nod_brainID(keep)),xticklabels(RegionList);xtickangle(45);ylim([0 1]);
% subplot(1,2,2);boxplot(pers,ClustID),xticklabels(clustnames);xtickangle(45);ylim([0 1]);

%% testing some more things

%%% checking relationship between dynamic variables

figure;
subplot(1,3,1);scatter(flexF20,node_cohesion,'filled');xlim([0 1]);ylim([0 1]);title('flexVScohe');
subplot(1,3,2);scatter(flexF20,P,'filled');xlim([0 1]);ylim([0 1]);title('flexVSprom');
subplot(1,3,3);scatter(node_cohesion,P,'filled');xlim([0 1]);ylim([0 1]);title('coheVSprom');


%%%%% looking for the top 

threshold=mean(node_cohesion)+1*std(node_cohesion);

idx_top=find(node_cohesion>threshold);

figure;
plot(Zbrain_brainMask2D(:,1),Zbrain_brainMask2D(:,2),'k');xlim([400 1350])
%  hold on;
%  gscatter(Nodes.Nod_coor(keep,1),Nodes.Nod_coor(keep,2),Nodes.Nod_brainID(keep,:),[],[],21); 
  hold on;
  gscatter(Nodes.Nod_coor(keep,1),Nodes.Nod_coor(keep,2),ClustID,'gggbrm',[],21); 
hold on;
scatter(Nodes.Nod_coor(keep(idx_top),1),Nodes.Nod_coor(keep(idx_top),2),20,'k','filled');
view(-90,90);
title('top 1SD cohesion');
hold off;


figure;
subplot(1,2,1);histogram(Nodes.Nod_brainID(keep)); xticks([1:9]); xticklabels(RegionList);xtickangle(45);hold on;histogram(Nodes.Nod_brainID(keep(idx_top)));
subplot(1,2,2);histogram(ClustID); xticks([1:6]); xticklabels(clustnames);xtickangle(45);hold on;histogram(ClustID(idx_top)); 



%%%%% looking for the top 

threshold=mean(P)+1*std(P);

idx_top2=find(P>threshold);

figure;
plot(Zbrain_brainMask2D(:,1),Zbrain_brainMask2D(:,2),'k');xlim([400 1350])
%  hold on;
%  gscatter(Nodes.Nod_coor(keep,1),Nodes.Nod_coor(keep,2),Nodes.Nod_brainID(keep,:),[],[],21); 
  hold on;
  gscatter(Nodes.Nod_coor(keep,1),Nodes.Nod_coor(keep,2),ClustID,'gggbrm',[],21); 
hold on;
scatter(Nodes.Nod_coor(keep(idx_top2),1),Nodes.Nod_coor(keep(idx_top2),2),20,'k','filled');
view(-90,90);
title('top 1SD prom');
hold off;


figure;
subplot(1,2,1);histogram(Nodes.Nod_brainID(keep)); xticks([1:9]); xticklabels(RegionList);xtickangle(45);hold on;histogram(Nodes.Nod_brainID(keep(idx_top2)));
subplot(1,2,2);histogram(ClustID); xticks([1:6]); xticklabels(clustnames);xtickangle(45);hold on;histogram(ClustID(idx_top2)); 

%%%% merging two


idx_top3=intersect(idx_top,idx_top2);

figure;
plot(Zbrain_brainMask2D(:,1),Zbrain_brainMask2D(:,2),'k');xlim([400 1350])
%  hold on;
%  gscatter(Nodes.Nod_coor(keep,1),Nodes.Nod_coor(keep,2),Nodes.Nod_brainID(keep,:),[],[],21); 
  hold on;
  gscatter(Nodes.Nod_coor(keep,1),Nodes.Nod_coor(keep,2),ClustID,'gggbrm',[],21); 
hold on;
scatter(Nodes.Nod_coor(keep(idx_top3),1),Nodes.Nod_coor(keep(idx_top3),2),20,'k','filled');
view(-90,90);
title('top 0SD promXcohe');
hold off;


figure;
subplot(1,2,1);histogram(Nodes.Nod_brainID(keep)); xticks([1:9]); xticklabels(RegionList);xtickangle(45);hold on;histogram(Nodes.Nod_brainID(keep(idx_top3)));
subplot(1,2,2);histogram(ClustID); xticks([1:6]); xticklabels(clustnames);xtickangle(45);hold on;histogram(ClustID(idx_top3)); 

