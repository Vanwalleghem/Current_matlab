CorrelationMatrix=squeeze(nanmean(CorrelationMatrices_perFish,3));
CorrelationMatrix=1-CorrelationMatrix;
CorrelationMatrix=weight_conversion(CorrelationMatrix,'autofix');
CorrelationMatrix=weight_conversion(CorrelationMatrix,'normalize');
figure;imagesc(CorrelationMatrix,[0 1]);
Correlation_bin_prop=threshold_proportional(CorrelationMatrix,0.25);
figure;imagesc(Correlation_bin_prop,[0 1]);
Correlation_bin=threshold_absolute(CorrelationMatrix,0.75);
figure;imagesc(Correlation_bin,[0 1]);


% Iterative community finetuning.
% W is the input connection matrix.
n  = size(Correlation_bin_prop,1);             % number of nodes
M  = 1:n;                   % initial community affiliations
Q0 = -1; Q1 = 0;            % initialize modularity values
while Q1-Q0>1e-5;           % while modularity increases
    Q0 = Q1;                % perform community detection
    [M, Q1] = community_louvain(Correlation_bin_prop, [], M);
end

color_Louvain=zeros(size(Node_all));
for i=1:8
    idx_temp=find(Node_ID==i);
    for node_nb=1:length(idx_temp)
        color_Louvain(idx_temp(node_nb),:)=colors(i,:)/256;
    end
end

Fighandle=figure;
set(Fighandle, 'Position', [1000,1000, 600, 1400]);
plot(Zbrain_brainMask2D(:,2),1400-Zbrain_brainMask2D(:,1),'k');hold on;
scatter(Node_all(:,2),Node_all(:,1),8,'k','filled');hold on;
scatter(Node_all(:,2),Node_all(:,1),5,M,'filled');hold on;

figure;
    Rw=rich_club_wu(Correlation_bin_prop);
    plot(Rw);


r=assortativity_wei(Correlation_bin_prop,0);
P=participation_coef(Correlation_bin_prop,Node_ID);
P_thr=participation_coef(Correlation_bin,Node_ID);

Node_all(find(round(Node_all(:,1))==109),2)=270;
Node_all(find(round(Node_all(:,1))==98),2)=275;
Node_all(find(round(Node_all(:,1))==100 & round(Node_all(:,2))==233),1)=147;


find(round(Node_all(:,1))==100 & round(Node_all(:,2))==233)

Fighandle=figure;
set(Fighandle, 'Position', [1000,1000, 600, 1400]);
plot(Zbrain_brainMask2D(:,2),1400-Zbrain_brainMask2D(:,1),'k');hold on;
scatter(Node_all(:,2),Node_all(:,1),55,'k','filled');hold on;
%scatter(Node_all(:,2),Node_all(:,1),50,color_temp,'filled');hold on;
%scatter(Node_all(:,2),Node_all(:,1),25,'k','filled');hold on;
scatter(Node_all(:,2),Node_all(:,1),45,P,'filled');colormap hot;hold off;

Fighandle=figure;
set(Fighandle, 'Position', [1000,1000, 600, 1400]);
plot(Zbrain_brainMask2D(:,2),1400-Zbrain_brainMask2D(:,1),'k');hold on;
scatter(Node_all(:,2),Node_all(:,1),30,'k','filled');hold on;
scatter(Node_all(:,2),Node_all(:,1),20,color_temp,'filled');hold off;

Fighandle=figure;
graph_temp=graph(Correlation_bin_prop,'upper');
set(Fighandle, 'Position', [1000,1000, 600, 1400]);
plot(Zbrain_brainMask2D(:,2),1400-Zbrain_brainMask2D(:,1),'k');hold on;
scatter(Node_all(:,2),Node_all(:,1),(degree(graph_temp))+5,'k','filled');hold on;
%plot(Nodes_graph,'NodeColor','k','LineStyle','None','Marker','o','MarkerSize',log(Nodes.density)+3,'XData',Node_all(:,2),'YData',Node_all(:,1),'NodeLabel',char.empty(100,0));hold on;
LWidths = 1*graph_temp.Edges.Weight/max(graph_temp.Edges.Weight);
LWidths(~isfinite(LWidths))=0.01;
LColors = ones(length(LWidths),3);LColors=LColors-LWidths;
plot(graph_temp,'EdgeColor',LColors,'NodeColor',color_temp,'LineStyle','-','LineWidth',LWidths,'Marker','o','MarkerSize',degree(graph_temp),'XData',Node_all(:,2),'YData',Node_all(:,1),'NodeLabel',char.empty(100,0));hold off;

[N,E] = rentian_scaling_3d(Correlation_bin_prop>0,Node_all,10000,1e-6);

[b,stats] = robustfit(log10(N),log10(E));[b(2) stats.se(2)]
figure; loglog(N,E,'.');axis([1 9 0 300]);


CorrelationMatrices_perFish_perRepeat=nan(size(Mean_allNodes_perFish,1),size(Mean_allNodes_perFish,1),size(Mean_allNodes_perFish,2),6);
for fish_nb=1:size(Mean_allNodes_perFish,2)
    for time_nb=1:6
        timing=50+(time_nb-1)*100;
        temp=squareform(pdist(squeeze(Mean_allNodes_perFish(:,fish_nb,timing:timing+70)),'correlation'));
        %temp=abs(1-temp);
        %temp=weight_conversion(temp,'autofix');
        %temp=weight_conversion(temp,'normalize');
        CorrelationMatrices_perFish_perRepeat(:,:,fish_nb,time_nb)=temp;
    end
end
CorrelationMatrix_PerRepeat=squeeze(nanmean(CorrelationMatrices_perFish_perRepeat,3));

Adjacency_matrix=CorrelationMatrix_PerRepeat;
for time_nb=1:6
    temp=squeeze(CorrelationMatrix_PerRepeat(:,:,time_nb));
    temp=abs(1-temp);
    temp=weight_conversion(temp,'autofix');
    temp=weight_conversion(temp,'normalize');
    Adjacency_matrix(:,:,time_nb)=temp;    
end


CorrelationMatrix_PerRepeat(isnan(CorrelationMatrix_PerRepeat))=1;

%%
%%%% testing some of the genlouvain functions

A=CorrelationMatrix_PerRepeat(:,:,1); %%% recovery loom

A(isnan(A))=0; %%%% I had to do this cause otherwise I just get NaN as an answer

gamma = 1;
k = full(sum(A)); %%% i had to change sum to nan sum because of the nan values
twom = sum(k);
B = full(A - gamma*k'*k/twom);
[S,Q] = genlouvain(B); %%%% you can play with some extra parameters here... 
Q = Q/twom

%% doing it for a multilayer matrix

gamma = 1.1;
omega = 0.95; %%% this makes a big influence...

N=(size(CorrelationMatrix_PerRepeat,1));
T=(size(CorrelationMatrix_PerRepeat,3));
B=spalloc(N*T,N*T,N*N*T+2*N*T);
twomu=0;
for s=1:T
    k=sum(CorrelationMatrix_PerRepeat(:,:,s));
    twom=sum(k);
    twomu=twomu+twom;
    indx=[1:N]+(s-1)*N;
    B(indx,indx)=CorrelationMatrix_PerRepeat(:,:,s)-gamma*k'*k/twom;
end
twomu=twomu+2*omega*N*(T-1);
B = B + omega*spdiags(ones(N*T,2),[-N,N],N*T,N*T);
[S,Q,nb_it] = iterated_genlouvain(B);
Q = Q/twomu
S = reshape(S,N,T);

figure;imagesc(S);
%%% to check community changes
for time_nb=1:6
    CommN(time_nb)=length(unique(S(:,time_nb)));
end
figure;plot(CommN);


CommCheck={};
figure;
N=(size(CorrelationMatrix_PerRepeat,1));
T=(size(CorrelationMatrix_PerRepeat,3));

gamma = 1;
omega = [0.05:0.05:1]; %%% this makes a big influence...
for o=1:length(omega)
    B=spalloc(N*T,N*T,N*N*T+2*N*T);
    twomu=0;
    for s=1:T
        k=sum(CorrelationMatrix_PerRepeat(:,:,s));
        twom=sum(k);
        twomu=twomu+twom;
        indx=[1:N]+(s-1)*N;
        B(indx,indx)=CorrelationMatrix_PerRepeat(:,:,s)-gamma*k'*k/twom;
    end
    
    
    
    twomu=twomu+2*omega(o)*N*(T-1);
    B = B + omega(o)*spdiags(ones(N*T,2),[-N,N],N*T,N*T);
    %[S,Q] = genlouvain(B);
    [S,Q,nb_it] = iterated_genlouvain(B);
    Q = Q/twomu;
    S = reshape(S,N,T);
    
    %%
    %%%% checking results
    
    subplot(4,5,o);imagesc(S);title(num2str(omega(o)));
    
    %%% to check community changes
    for time_nb=1:6
        CommN(time_nb)=length(unique(S(:,time_nb)));
        CommMax(time_nb)=max(unique(S(:,time_nb)));
    end
    %figure;plot(CommN);
    Omtemp(o)=omega(o);
    Qtemp(o)=Q;
    CommNtemp(o,:)=CommN;
    CommMaxtemp(o,:)=CommMax;
end

% CommCheck{1}=Omtemp;
% CommCheck{2}=Qtemp;
% CommCheck{3}=CommNtemp;
% CommCheck{4}=CommMaxtemp;

figure;
subplot(1,3,1);plot(Omtemp,Qtemp);title('Q'); 
subplot(1,3,2);plot(Omtemp,mean(CommNtemp,2));ylim([0 10]);title('average #Comm'); 
subplot(1,3,3);plot(Omtemp,max(CommMaxtemp,[],2));title('Max Comm');

Node_all_temp=round(Node_all);
Node_all_temp(:,3)=Node_all_temp(:,3)/2;
% Create rotation matrix
theta = 180; % to rotate 90 counterclockwise
R = [cosd(theta) -sind(theta); sind(theta) cosd(theta)];
Node_all_temp(:,[1 2])=Node_all_temp(:,[1 2])*R;
Node_all_temp(:,1)=1406+Node_all_temp(:,1);
Node_all_temp(:,2)=621+Node_all_temp(:,2);

idx_rand=randsample(length(Zbrain_AllMask),10000); %need to select a few members only
Fighandle=figure;
set(Fighandle, 'Position', [10,10, 600, 1400]);
scatter(Zbrain_AllMask(idx_rand,2),Zbrain_AllMask(idx_rand,1),'.')
hold on;
scatter(Node_all_temp(:,2),Node_all_temp(:,1),50,'r','filled'); hold off;

%NEED TO ROTATE 180 !!!

Node_region={};
for i=1:length(Node_ID)
    List_region=[];
    for j=1:length(Zbrain_Masks)        
        IsInEyes_temp=ismember(Zbrain_Masks{j,3},round(Node_all_temp(i,:)),'rows');IsInEyes_temp=find(IsInEyes_temp==1);
        if IsInEyes_temp
            List_region=[List_region j];
        end
    end
    Node_region{i,1}=List_region;
    if List_region
        Node_region{i,2}=strcat(Zbrain_Masks{List_region(1),1},' ',Zbrain_Masks{List_region(1),2});
    end
end

for i=1:length(Node_ID)    
    List_region=Node_region{i,1};
    if List_region
        for ij=1:length(List_region)
            Node_region{i,1+ij}=strcat(Zbrain_Masks{List_region(ij),1},'_',Zbrain_Masks{List_region(ij),2});
        end
    end
end


for i=1:length(Node_ID)    
    List_region=Node_region{i,1};
    if isempty(List_region)
        Min_dist=[inf,0];
        for j=1:length(Zbrain_Masks)
            dist_temp=min(pdist2(round(Node_all_temp(i,:)),Zbrain_Masks{j,3}));
            if dist_temp<Min_dist(1)
                Min_dist(1)=dist_temp;
                Min_dist(2)=j;
            end
        end
        Node_region{i,1}=Min_dist(2);
        Node_region{i,2}=strcat(Zbrain_Masks{Min_dist(2),1},'_',Zbrain_Masks{Min_dist(2),2});
    end
end


