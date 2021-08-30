
%%%%% this script is to try to find the optimal gamma and omega values


Qopt=[];

gamma = 0.1:0.4:40;
omega = 0.1:0.4:40;

for test=1:100
Qmat=[];
parfor i=1:100   
for j=1:100

N=length(A{1});
T=length(A);
B=spalloc(N*T,N*T,N*N*T+2*N*T);
twomu=0;
for s=1:T
    k=sum(A{s});
    twom=sum(k);
    twomu=twomu+twom;
    indx=[1:N]+(s-1)*N;
    B(indx,indx)=A{s}-gamma(i)*k'*k/twom;
end
twomu=twomu+2*omega(j)*N*(T-1);
B = B + omega(j)*spdiags(ones(N*T,2),[-N,N],N*T,N*T);
[S,Q] = genlouvain(B,10000,0);
%[S,Q,nb_it] = iterated_genlouvain(B);
Q = Q/twomu;
%S = reshape(S,N,T);

Qmat(i,j)=Q;
   
end
end

Qopt=cat(3,Qopt,Qmat);

end

meanQopt=mean(Qopt,3);
varQopt=var(Qopt,0,3);

figure;
subplot(1,2,1);imagesc(meanQopt); colorbar;
subplot(1,2,2);imagesc(varQopt); colorbar;

%% close up



Qopt2=[];

gamma = 0.1:0.1:3.7;
omega = 0.1:0.1:3.7;

for test=1:100
Qmat2=[];
for i=1:length(gamma)   
parfor j=1:length(omega)

N=length(A{1});
T=length(A);
B=spalloc(N*T,N*T,N*N*T+2*N*T);
twomu=0;
for s=1:T
    k=sum(A{s});
    twom=sum(k);
    twomu=twomu+twom;
    indx=[1:N]+(s-1)*N;
    B(indx,indx)=A{s}-gamma(i)*k'*k/twom;
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

meanQopt2=mean(Qopt2,3);
varQopt2=var(Qopt2,0,3);

figure;
subplot(1,2,1);imagesc(meanQopt2); colorbar;
subplot(1,2,2);imagesc(varQopt2); colorbar;




%%
%save('QoptF20.mat','Qmat','Qopt','meanQopt','varQopt','Qmat2','Qopt2','meanQopt2','varQopt2');


%%


%%%% making a null nodal model of A

randOrder = randperm(99,99);

%unique(randOrder)


An_null={};

for loom=1:31
    
    randOrder = randperm(99,99);
    
   temp=Data_corrMat2.f20.Mean_corrMat{1,loom}(keep(randOrder),keep(randOrder));
   temp(isnan(temp))=0;
    An_null{loom}=temp;  
end
clear temp



Qopt2_null=[];

gamma = 0.1:0.1:3.7;
omega = 0.1:0.1:3.7;

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

figure;
subplot(1,2,1);imagesc(meanQopt2_null); colorbar;
subplot(1,2,2);imagesc(varQopt2_null); colorbar;


figure;imagesc(meanQopt2-meanQopt2_null);colorbar;

save('QoptF20.mat','A','Qmat','Qopt','meanQopt','varQopt','Qmat2','Qopt2','meanQopt2','varQopt2','An_null','randOrder','Qmat2_null','Qopt2_null','meanQopt2_null','varQopt2_null');


%%

%%% trying to find the minimum variance of Q and maximum difference in Q
%%% with the null model:

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

meanQopt2_dif_norm-varQopt2_norm;

%% trying also with a temporal null model



%%%%% this is changing the order of the layers


randOrder2 = randperm(31,31);

%unique(randOrder2)


At_null={};

for loom=1:31
    
   temp=Data_corrMat2.f20.Mean_corrMat{1,randOrder2(loom)}(keep,keep);
   temp(isnan(temp))=0;
    At_null{loom}=temp;  
end
clear temp



Qopt2_null2=[];

gamma = 0.1:0.1:3.7;
omega = 0.1:0.1:3.7;

for test=1:100
Qmat2_null2=[];
for i=1:length(gamma)   
parfor j=1:length(omega)

N=length(At_null{1});
T=length(At_null);
B=spalloc(N*T,N*T,N*N*T+2*N*T);
twomu=0;
for s=1:T
    k=sum(At_null{s});
    twom=sum(k);
    twomu=twomu+twom;
    indx=[1:N]+(s-1)*N;
    B(indx,indx)=At_null{s}-gamma(i)*k'*k/twom;
end
twomu=twomu+2*omega(j)*N*(T-1);
B = B + omega(j)*spdiags(ones(N*T,2),[-N,N],N*T,N*T);
[S,Q] = genlouvain(B,10000,0);
%[S,Q,nb_it] = iterated_genlouvain(B);
Q = Q/twomu;
%S = reshape(S,N,T);

Qmat2_null2(i,j)=Q;
   
end
end

Qopt2_null2=cat(3,Qopt2_null2,Qmat2_null2);

end

meanQopt2_null2=mean(Qopt2_null2,3);
varQopt2_null2=var(Qopt2_null2,0,3);

figure;
subplot(1,2,1);imagesc(meanQopt2_null2); colorbar;
subplot(1,2,2);imagesc(varQopt2_null2); colorbar;


figure;imagesc(meanQopt2-meanQopt2_null2);colorbar;

%%
save('QoptF20.mat','A','Qmat','Qopt','meanQopt','varQopt','Qmat2','Qopt2','meanQopt2','varQopt2','An_null','randOrder','Qmat2_null','Qopt2_null','meanQopt2_null','varQopt2_null','At_null','randOrder2','Qmat2_null2','Qopt2_null2','meanQopt2_null2','varQopt2_null2');

%%


%%% trying to find the minimum variance of Q and maximum difference in Q
%%% with the null2 model:

meanQopt2_dif2=meanQopt2-meanQopt2_null2;
low=min(min(meanQopt2_dif2));
high=max(max(meanQopt2_dif2));

meanQopt2_dif_norm2=meanQopt2_dif2-low;
meanQopt2_dif_norm2=meanQopt2_dif_norm2/(high-low);


varQopt2_norm2=varQopt2;
low=min(min(varQopt2_norm2));
high=max(max(varQopt2_norm2));

varQopt2_norm2=varQopt2_norm2-low;
varQopt2_norm2=varQopt2_norm2/(high-low);

meanQopt2_dif_norm2-varQopt2_norm2;

