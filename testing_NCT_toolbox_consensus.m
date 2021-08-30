
%%%% this script is to test the concensus functions from Danielle Bassett and her group. 
%%%% the idea is running multiple times (like 100 times) the multilayer network analysis and
%%%% then using these functions to find a representative partition that we
%%%% can then analyse furhter. 

for test=1:5
    
   
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

%figure;imagesc(S);colormap('jet');

%%
%     
%  
% flexF20=flexibility(S','temp'); %%% i need to mind the orientation of the matrix to do it properly
%  
% options.figureFlag	= 0;
% options.colormap	= 'jet';
% 
% [Cij,node_cohesion,node_disjoint,node_flexibility,strength_cohesion,commChanges,commCohesion,commDisjoint,commIndex] = calc_node_cohesion(S,options);
% 
% P = promiscuity(S');  %%% i need to mind the orientation of the matrix to do it properly
% 
% %%
% 
% 
% figure;
% subplot(1,3,1);scatter(flexF20,node_cohesion,'filled');xlim([0 1]);ylim([0 1]);title('flexVScohe');
% subplot(1,3,2);scatter(flexF20,P,'filled');xlim([0 1]);ylim([0 1]);title('flexVSprom');
% subplot(1,3,3);scatter(node_cohesion,P,'filled');xlim([0 1]);ylim([0 1]);title('coheVSprom');
% 
% 
% %%%%% looking for the top 
% 
% threshold=mean(node_cohesion)+1*std(node_cohesion);
% 
% idx_top=find(node_cohesion>threshold);
% 
% 
% figure;sgtitle('1SD cohesion');
% subplot(1,2,1);histogram(Nodes.Nod_brainID(keep)); xticks([1:9]); xticklabels(RegionList);xtickangle(45);hold on;histogram(Nodes.Nod_brainID(keep(idx_top)));
% subplot(1,2,2);histogram(ClustID); xticks([1:6]); xticklabels(clustnames);xtickangle(45);hold on;histogram(ClustID(idx_top)); 
% 
% 
% threshold=mean(P)+1*std(P);
% 
% idx_top2=find(P>threshold);
% 
% 
% figure;sgtitle('1SD prom');
% subplot(1,2,1);histogram(Nodes.Nod_brainID(keep)); xticks([1:9]); xticklabels(RegionList);xtickangle(45);hold on;histogram(Nodes.Nod_brainID(keep(idx_top2)));
% subplot(1,2,2);histogram(ClustID); xticks([1:6]); xticklabels(clustnames);xtickangle(45);hold on;histogram(ClustID(idx_top2)); 

time=[1 2 6 11 12];

low=min(min(S));
high=max(max(S))+1;


  figure;
  
  subplot(1,6,1);imagesc(S);colormap('gray');caxis([low high]);
  temp=[];
  for i=1:5
  subplot(1,6,i+1);
  plot(Zbrain_brainMask2D(:,1),Zbrain_brainMask2D(:,2),'k');xlim([400 1350])
  hold on;
%   gscatter(Nodes.Nod_coor(keep,1),Nodes.Nod_coor(keep,2),ClustID,'gggbrm',[],21,'off'); 
%     hold on;
    scatter(Nodes.Nod_coor(keep,1),Nodes.Nod_coor(keep,2),20,S(:,time(i)),'filled');colormap('gray');caxis([low high]);
    view(-90,90);
    %title('o=0.8/test=',num2str(test),'/time=',num2str(time(i)));
    hold off;  

    
    
  end
  
  
end


%% testing the converging partions 

%%%% I will try the consensus iterative from Danielle Bassett
%%% found at: http://commdetect.weebly.com/

%%%% I will first do it for the timepoints I was testing... later for the
%%%% whole 30 looms. 

%%% getting a sample of optimizations for the pre, 1st, 5th, 10th and 11th



%%%% making a multidimensional structure where to store things
S_test=[];

for test=1:100
     
gamma = 1;
omega = 1.3; %%% this makes a big influence... 

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
[S,Q] = genlouvain(B,10000,0);
%[S,Q,nb_it] = iterated_genlouvain(B);
Q = Q/twomu;
S = reshape(S,N,T);

S_test=cat(3,S_test,S);
  
  
end


%% trying consensus_iterative 


S_good=[];
for i=1:31
   C=S_test(:,i,:); 
   C=squeeze(C);
   [S2, Q2, X_new3, qpc] = consensus_iterative(C');
    S_good(:,i)=S2(i,:);
end



low=min(min(S_good));
high=max(max(S_good))+1;

time=[1 2 6 11 12];

  figure;
  
  subplot(1,7,1);imagesc(S_good);colormap('jet');caxis([low high]);
  temp=[];
  for i=1:5
  subplot(1,7,i+1);
  plot(Zbrain_brainMask2D(:,1),Zbrain_brainMask2D(:,2),'k');xlim([400 1350])
  hold on;
%    gscatter(Nodes.Nod_coor(keep,1),Nodes.Nod_coor(keep,2),ClustID,'gggbrm',[],21,'off'); 
%      hold on;
    scatter(Nodes.Nod_coor(keep,1),Nodes.Nod_coor(keep,2),20,S_good(:,time(i)),'filled');colormap('jet');caxis([low high]);
    view(-90,90);
    %title('o=0.8/test=',num2str(test),'/time=',num2str(time(i)));
    hold off;  
     
  end
  subplot(1,7,7);
  plot(Zbrain_brainMask2D(:,1),Zbrain_brainMask2D(:,2),'k');xlim([400 1350])
  hold on;
   gscatter(Nodes.Nod_coor(keep,1),Nodes.Nod_coor(keep,2),ClustID,'gggbrm',[],20,'off');     
     view(-90,90);
    %title('o=0.8/test=',num2str(test),'/time=',num2str(time(i)));
    hold off;  
    
    
    %%%% to check some of the green nodes recluted to other communites responses
    idx_temp=find(S_good(:,25)==3);
    figure;
    for i=1:length(idx_temp)
        plot(avgWTmat(idx_temp(i),:));
        pause(2);
        
    end

    
    
    %% trying consensus_similarity 
    
    
S_good_sim=[];
for i=1:31
   C=S_test(:,i,:); 
   C=squeeze(C);
   [consensus, consensus_simm, pairwise_simm] = consensus_similarity(C');
    S_good_sim(:,i)=S2(i,:);
end
    


low=min(min(S_good_sim));
high=max(max(S_good_sim))+1;

time=[1 2 6 11 12];

  figure;
  
  subplot(1,7,1);imagesc(S_good_sim);colormap('jet');caxis([low high]);
  temp=[];
  for i=1:5
  subplot(1,7,i+1);
  plot(Zbrain_brainMask2D(:,1),Zbrain_brainMask2D(:,2),'k');xlim([400 1350])
  hold on;
%    gscatter(Nodes.Nod_coor(keep,1),Nodes.Nod_coor(keep,2),ClustID,'gggbrm',[],21,'off'); 
%      hold on;
    scatter(Nodes.Nod_coor(keep,1),Nodes.Nod_coor(keep,2),20,S_good_sim(:,time(i)),'filled');colormap('jet');caxis([low high]);
    view(-90,90);
    %title('o=0.8/test=',num2str(test),'/time=',num2str(time(i)));
    hold off;  
     
  end
  subplot(1,7,7);
  plot(Zbrain_brainMask2D(:,1),Zbrain_brainMask2D(:,2),'k');xlim([400 1350])
  hold on;
   gscatter(Nodes.Nod_coor(keep,1),Nodes.Nod_coor(keep,2),ClustID,'gggbrm',[],20,'off');     
     view(-90,90);
    %title('o=0.8/test=',num2str(test),'/time=',num2str(time(i)));
    hold off;  
    
      
    %%%% to check some of the green nodes recluted to other communites responses
    idx_temp=find(S_good_sim(:,25)==3);
    
    figure;
    for i=1:length(idx_temp)
        plot(avgWTmat(idx_temp(i),:));ylim([-2 12]);
        pause(2);
        
    end

    
    