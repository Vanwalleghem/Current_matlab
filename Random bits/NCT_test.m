Adjacency_matrix=abs(1-CorrelationMatrix_PerRepeat);
Adjacency_matrix(isnan(Adjacency_matrix))=0;
Qopt2=[];

gamma = 0.05:0.1:3;
omega = 0.05:0.1:3;

for test=1:100
    Qmat2=[];
    for i=1:length(gamma)
        parfor j=1:length(omega)           
            N=(size(Adjacency_matrix,1));
            T=(size(Adjacency_matrix,3));
            B=spalloc(N*T,N*T,N*N*T+2*N*T);
            twomu=0;
            for s=1:T
                k=sum(Adjacency_matrix(:,:,s));
                twom=sum(k);
                twomu=twomu+twom;
                indx=[1:N]+(s-1)*N;
                B(indx,indx)=Adjacency_matrix(:,:,s)-gamma(i)*k'*k/twom;
            end
            twomu=twomu+2*omega(j)*N*(T-1);
            B = B + omega(j)*spdiags(ones(N*T,2),[-N,N],N*T,N*T);
            [S,Q] = genlouvain(B,10000,0);
            %[S,Q,nb_it] = iterated_genlouvain(B);
            Q = Q/twomu;
            %S = reshape(S,N,T);
            
            Qmat2(i,j)=Q;
            
        end
    end
    
    Qopt2=cat(3,Qopt2,Qmat2);
    
end

meanQopt2=nanmean(Qopt2,3);
varQopt2=nanvar(Qopt2,0,3);

Fighandle=figure;
set(Fighandle, 'Position', [10,10, 1200, 300]);
subplot(1,2,1);imagesc(meanQopt2); colorbar;
subplot(1,2,2);imagesc(varQopt2); colorbar;

%Nodal null model
An_null={};

for s=1:T    
    randOrder = randperm(size(Adjacency_matrix,1),size(Adjacency_matrix,1));    
    temp=Adjacency_matrix(randOrder,randOrder,s);
    temp(isnan(temp))=0;
    An_null{s}=temp;
end
clear temp


Qopt2_null=[];
for test=1:100
    Qmat2_null=[];
    for i=1:length(gamma)
        parfor j=1:length(omega)            
            N=length(An_null{1});
            T=length(An_null);
            B=spalloc(N*T,N*T,N*N*T+2*N*T);
            twomu=0;
            for s=1:T
                k=sum(An_null{s});
                twom=sum(k);
                twomu=twomu+twom;
                indx=[1:N]+(s-1)*N;
                B(indx,indx)=An_null{s}-gamma(i)*k'*k/twom;
            end
            twomu=twomu+2*omega(j)*N*(T-1);
            B = B + omega(j)*spdiags(ones(N*T,2),[-N,N],N*T,N*T);
            [S,Q] = genlouvain(B,10000,0);
            %[S,Q,nb_it] = iterated_genlouvain(B);
            Q = Q/twomu;
            %S = reshape(S,N,T);
            
            Qmat2_null(i,j)=Q;
            
        end
    end
    
    Qopt2_null=cat(3,Qopt2_null,Qmat2_null);
    
end

meanQopt2_null=mean(Qopt2_null,3);
varQopt2_null=var(Qopt2_null,0,3);

Fighandle=figure;
set(Fighandle, 'Position', [10,10, 1200, 300]);
subplot(1,2,1);imagesc(meanQopt2_null); colorbar;
subplot(1,2,2);imagesc(varQopt2_null); colorbar;

Fighandle=figure;
set(Fighandle, 'Position', [10,10, 1200, 300]);
subplot(1,2,1);imagesc(meanQopt2-meanQopt2_null); 
subplot(1,2,2);imagesc(varQopt2); 

Fighandle=figure;
set(Fighandle, 'Position', [10,10, 400, 400]);
imagesc((meanQopt2-meanQopt2_null).*(-(varQopt2-max(varQopt2(:))))); 
opt_temp=(meanQopt2-meanQopt2_null).*(-(varQopt2-max(varQopt2(:))));
[M,row]=max(opt_temp,[],1);
[~,col]=max(M(1:20));
opt_temp(row(col),col)
gamma_opt=gamma(row(col))
omega_opt=omega(col)

gamma_opt=1.1
omega_opt=1.35

meanQopt2_dif=meanQopt2-meanQopt2_null;
low=min(min(meanQopt2_dif));
high=max(max(meanQopt2_dif));

meanQopt2_dif_norm=meanQopt2_dif-low;
meanQopt2_dif_norm=meanQopt2_dif_norm/(high-low);


varQopt2_norm=varQopt2;
low=min(min(varQopt2_norm));
high=max(max(varQopt2_norm));

varQopt2_norm=varQopt2_norm-low;
varQopt2_norm=varQopt2_norm/(high-low);

figure;imagesc(meanQopt2_dif_norm-varQopt2_norm);


%%%% making a multidimensional structure where to store things
S_test=[];


for test=1:100
    
    N=(size(Adjacency_matrix,1));
    T=(size(Adjacency_matrix,3));
    B=spalloc(N*T,N*T,N*N*T+2*N*T);
    twomu=0;
    for s=1:T
        k=sum(Adjacency_matrix(:,:,s));
        twom=sum(k);
        twomu=twomu+twom;
        indx=[1:N]+(s-1)*N;
        B(indx,indx)=Adjacency_matrix(:,:,s)-gamma_opt*k'*k/twom;
    end
    twomu=twomu+2*omega_opt*N*(T-1);
    B = B + omega_opt*spdiags(ones(N*T,2),[-N,N],N*T,N*T);
    [S,Q] = genlouvain(B,10000,0);
    %[S,Q,nb_it] = iterated_genlouvain(B);
    Q = Q/twomu;
    S = reshape(S,N,T);
    
    S_test=cat(3,S_test,S);
    
    
end


S_good=[];
for i=1:T
   C=S_test(:,i,:); 
   C=squeeze(C);
   [S2, Q2, X_new3, qpc] = consensus_iterative(C');
    S_good(:,i)=S2(i,:);
end



low=min(min(S_good));
high=max(max(S_good))+1;

figure;imagesc(S_good)


   
S_good_sim=[];
for i=1:T
   C=S_test(:,i,:); 
   C=squeeze(C);
   [consensus, consensus_simm, pairwise_simm] = consensus_similarity(C');
    S_good_sim(:,i)=S2(i,:);
end
    


low=min(min(S_good_sim));
high=max(max(S_good_sim))+1;
figure;imagesc(S_good_sim)



Fighandle=figure;
set(Fighandle, 'Position', [5,5, 2500, 700]);
subplot(1,8,1);imagesc(S_good_sim);colormap('jet');caxis([low high]);
temp=[];
for i=1:6
    subplot(1,8,i+1);
    
    plot(Zbrain_brainMask2D(:,2),1400-Zbrain_brainMask2D(:,1),'k');hold on;
    hold on;    
    scatter(Node_all(:,2),Node_all(:,1),20,S_good_sim(:,i),'filled');colormap(Qual_color);caxis([low high]);
    %view(-90,90);    
    hold off;    
end
subplot(1,8,8);
plot(Zbrain_brainMask2D(:,2),1400-Zbrain_brainMask2D(:,1),'k');hold on;
hold on;
scatter(Node_all(:,2),Node_all(:,1),20,color_temp,'filled');
%view(-90,90);
%title('o=0.8/test=',num2str(test),'/time=',num2str(time(i)));
hold off;

Flexibility_opt=flexibility(S_good_sim');

Fighandle=figure;
set(Fighandle, 'Position', [1000,1000, 600, 1400]);
plot(Zbrain_brainMask2D(:,2),1400-Zbrain_brainMask2D(:,1),'k');hold on;
scatter(Node_all(:,2),Node_all(:,1),65,'k','filled');hold on;
scatter(Node_all(:,2),Node_all(:,1),45,Flexibility_opt,'filled');colormap hot;hold off;

[Cij,node_cohesion,node_disjoint,node_flexibility,strength_cohesion,commChanges,commCohesion,commDisjoint,commIndex] = calc_node_cohesion(S_good_sim);

node_promiscuity=promiscuity(S_good_sim');

Fighandle=figure;
set(Fighandle, 'Position', [5,5, 2200, 600]);
subplot(1,4,1);
plot(Zbrain_brainMask2D(:,2),1400-Zbrain_brainMask2D(:,1),'k');hold on;
scatter(Node_all(:,2),Node_all(:,1),65,'k','filled');hold on;
scatter(Node_all(:,2),Node_all(:,1),45,node_flexibility,'filled');colormap hot;hold off;title('flexibility');
subplot(1,4,2);
plot(Zbrain_brainMask2D(:,2),1400-Zbrain_brainMask2D(:,1),'k');hold on;
scatter(Node_all(:,2),Node_all(:,1),65,'k','filled');hold on;
scatter(Node_all(:,2),Node_all(:,1),45,node_cohesion,'filled');colormap hot;hold off;title('cohesion');
subplot(1,4,3);
plot(Zbrain_brainMask2D(:,2),1400-Zbrain_brainMask2D(:,1),'k');hold on;
scatter(Node_all(:,2),Node_all(:,1),65,'k','filled');hold on;
scatter(Node_all(:,2),Node_all(:,1),45,node_disjoint,'filled');colormap hot;hold off;title('disjoint');
subplot(1,4,4);
plot(Zbrain_brainMask2D(:,2),1400-Zbrain_brainMask2D(:,1),'k');hold on;
scatter(Node_all(:,2),Node_all(:,1),65,'k','filled');hold on;
scatter(Node_all(:,2),Node_all(:,1),45,node_promiscuity,'filled');colormap hot;hold off;title('promiscuity');


