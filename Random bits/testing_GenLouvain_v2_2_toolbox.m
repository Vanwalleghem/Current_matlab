%%%%% this script is to play with the GenLouvain toolbox following Dani's
%%%%% advice:

% 
% Hi Gilles and Emmanuel, 
% 
% Just following up on our discussion last week about dynamic community detection,
% so that we can assess quantitative markers of the modularity and its change over the looms (and fish). 
% The main code bank you want is this one:
% 
% http://netwiki.amath.unc.edu/GenLouvain/GenLouvain
% 
% You'll want to do a combination of Examples 2 and 3. Example 2 will show you how to link up looms (time) 
% in an ordered way. Example 3 will show you how to link up fish in a categorical way. 
% You'll need to choose 3 parameters: gamma, omega for the ordered links between times, 
% and omega for the categorical links between fish. We usually set gamma = 1 on the first pass unless 
% we have a reason not to (which in this case we don't); a large gamma will give you a few small communities, 
% and a small gamma will give you a few large communities. Next, you want to choose the two omegas, 
% which will determine the allowed variance in solutions across time and across fish, respectively. 
% Do you expect the variance in communities across time to be larger than across fish? 
% If so, then you want to choose the omega for fish (say 2?) to be larger than the omega for time (say 0.5?).
% The S matrix is a perfect way to visualize how the community partitions are different across fish and time. 
% Can you email us some of those when you get them? That will help us determine if the parameter values are reasonable.
% 
% Assuming they are, then the next step would be to calculate some quantitative statistics, which you can find here:
% http://commdetect.weebly.com/
% 
% under SectionIII: Dynamic Community Structure
% 
% Also, the papers I cited in my initial email about this would be relevant to read.
% 
% All the Best,
% 
% Dani

%%
%%%% I will first try it with the f20 dataset

load('Nodes_WTnFmr1Exp.mat')
load('Nodes_FnS_all_samples.mat')
load('graphs_fmr1Exp.mat')
load('graphs_FnS_all.mat')


%%
%%% this is how some of the mooments look. 
%%% note: the "keep" must be the adecuate one for this dataset (length=99)
figure;

counter=1;
subplot(1,8,counter);imagesc(Data_corrMat2.f20.Mean_corrMat{1,1}(keep,keep)); caxis([-1 1])%% for pre loom
subplot(1,8,counter+1);imagesc(Data_corrMat2.f20.Mean_corrMat{1,2}(keep,keep));caxis([-1 1]) %% for 1st loom
subplot(1,8,counter+2);imagesc(Data_corrMat2.f20.Mean_corrMat{1,3}(keep,keep));caxis([-1 1])
subplot(1,8,counter+3);imagesc(Data_corrMat2.f20.Mean_corrMat{1,4}(keep,keep));caxis([-1 1])
subplot(1,8,counter+4);imagesc(Data_corrMat2.f20.Mean_corrMat{1,5}(keep,keep));caxis([-1 1])
subplot(1,8,counter+5);imagesc(Data_corrMat2.f20.Mean_corrMat{1,6}(keep,keep));caxis([-1 1])
subplot(1,8,counter+6);imagesc(Data_corrMat2.f20.Mean_corrMat{1,11}(keep,keep));caxis([-1 1])%% for 10th loom
subplot(1,8,counter+7);imagesc(Data_corrMat2.f20.Mean_corrMat{1,12}(keep,keep));caxis([-1 1]) %% for 11th loom
title('f20');

%%
%%%% testing some of the genlouvain functions

A=Data_corrMat2.f20.Mean_corrMat{1,12}(keep,keep); %%% recovery loom

A(isnan(A))=0; %%%% I had to do this cause otherwise I just get NaN as an answer

gamma = 1;
k = full(sum(A)); %%% i had to change sum to nan sum because of the nan values
twom = sum(k);
B = full(A - gamma*k'*k/twom);
[S,Q] = genlouvain(B); %%%% you can play with some extra parameters here... 
Q = Q/twom

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

%% 
%%%% checking results

figure;imagesc(S);

%%% to check community changes
for loom=1:31
   CommN(loom)=length(unique(S(:,loom))); 
end
%figure;plot(CommN);


%% per fish?

%%% now I will try per fish... but not sure how to do it. I will try at a
%%% particular timepoint... loom 11? or moment 1 to 13?

for loom=[1 2 3 11 12]
%%% making the cell array

fish=fieldnames(Data_corrMat2.f20);
fish=fish(1:end-1);

Af={};

for f=1:length(fish)
   temp=Data_corrMat2.f20.(fish{f}).loomsR{1,loom}(keep,keep);
   temp(isnan(temp))=0;
    Af{f}=temp;  
end
clear temp



%%

%%% you need to set values for gamma and omega

gamma = 1;
omega = 0.5; %%% this makes a big influence... 

N=length(Af{1});
T=length(Af);
Bf=spalloc(N*T,N*T,(N+T)*N*T);
twomu=0;
for s=1:T
    k=sum(Af{s});
    twom=sum(k);
    twomu=twomu+twom;
    indx=[1:N]+(s-1)*N;
    Bf(indx,indx)=Af{s}-gamma*k'*k/twom;
end
twomu=twomu+T*omega*N*(T-1);
all2all = N*[(-T+1):-1,1:(T-1)];
Bf = Bf + omega*spdiags(ones(N*T,2*T-2),all2all,N*T,N*T);
[Sf,Qf,nb_it_f] = iterated_genlouvain(Bf);
Qf = Qf/twomu
Sf = reshape(Sf,N,T);


%% 
%%%% checking results

figure;imagesc(Sf);

%%% to check community changes
for f=1:11
   CommNf(f)=length(unique(Sf(:,f))); 
end
%figure;plot(CommNf);
end



%% playing with omega values

%%% I will check the number of communities and how many new ones appear
%%% when I change omega... 

%%% first for the f20 habituation


gamma = 1;
omega = [0.1:0.1:2]; %%% this makes a big influence... 

CommCheck={};
figure;
for o=1:length(omega)
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
for loom=1:31
   CommN(loom)=length(unique(S(:,loom))); 
   CommMax(loom)=max(unique(S(:,loom)));
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

%% now with thresholded graphs


At={};

for loom=1:31
   temp=Data_corrMat2.f20.Mean_corrMat{1,loom}(keep,keep);
   temp(isnan(temp))=0;
   temp = threshold_absolute(abs(temp),0.75);  
   At{loom}=temp;  
end
clear temp


%%% this is how some of the mooments look. 
%%% note: the "keep" must be the adecuate one for this dataset (length=99)
figure;

counter=1;
subplot(1,8,counter);imagesc(At{1,1}); caxis([0 1])%% for pre loom
subplot(1,8,counter+1);imagesc(At{1,2});caxis([0 1]) %% for 1st loom
subplot(1,8,counter+2);imagesc(At{1,3});caxis([0 1])
subplot(1,8,counter+3);imagesc(At{1,4});caxis([0 1])
subplot(1,8,counter+4);imagesc(At{1,5});caxis([0 1])
subplot(1,8,counter+5);imagesc(At{1,6});caxis([0 1])
subplot(1,8,counter+6);imagesc(At{1,11});caxis([0 1])%% for 10th loom
subplot(1,8,counter+7);imagesc(At{1,12});caxis([0 1]) %% for 11th loom
title('f20 thresholded');


%%% I will check the number of communities and how many new ones appear
%%% when I change omega... 

%%% first for the f20 habituation


gamma = 1;
omega = [0.1:0.1:2]; %%% this makes a big influence... 

CommCheck={};
figure;
for o=1:length(omega)
N=length(At{1});
T=length(At);
B=spalloc(N*T,N*T,N*N*T+2*N*T);
twomu=0;
for s=1:T
    k=sum(At{s});
    twom=sum(k);
    twomu=twomu+twom;
    indx=[1:N]+(s-1)*N;
    B(indx,indx)=At{s}-gamma*k'*k/twom;
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
for loom=1:31
   CommN(loom)=length(unique(S(:,loom))); 
   CommMax(loom)=max(unique(S(:,loom)));
end
%figure;plot(CommN);
OmtempT(o)=omega(o);
QtempT(o)=Q;
CommNtempT(o,:)=CommN;
CommMaxtempT(o,:)=CommMax;
end

% CommCheck{1}=Omtemp;
% CommCheck{2}=Qtemp;
% CommCheck{3}=CommNtemp;
% CommCheck{4}=CommMaxtemp;

figure;
subplot(1,3,1);plot(Omtemp,Qtemp);ylim([0 1]); hold on; plot(OmtempT,QtempT);title('Q'); hold off;
subplot(1,3,2);plot(Omtemp,mean(CommNtemp,2));ylim([0 10]);hold on; plot(OmtempT,mean(CommNtempT,2));title('average #Comm'); hold off;
subplot(1,3,3);plot(Omtemp,max(CommMaxtemp,[],2));ylim([0 16]);hold on;plot(OmtempT,max(CommMaxtempT,[],2));title('Max Comm'); hold off;

%% now for fish, at loom 11


loom=12;
Af={};

for f=1:length(fish)
   temp=Data_corrMat2.f20.(fish{f}).loomsR{1,loom}(keep,keep);
   temp(isnan(temp))=0;
    Af{f}=temp;  
end
clear temp

gamma = 1;
omega = [0.1:0.1:2]; %%% this makes a big influence... 


figure;
for o=1:length(omega)

N=length(Af{1});
T=length(Af);
Bf=spalloc(N*T,N*T,(N+T)*N*T);
twomu=0;
for s=1:T
    k=sum(Af{s});
    twom=sum(k);
    twomu=twomu+twom;
    indx=[1:N]+(s-1)*N;
    Bf(indx,indx)=Af{s}-gamma*k'*k/twom;
end
twomu=twomu+T*omega(o)*N*(T-1);
all2all = N*[(-T+1):-1,1:(T-1)];
Bf = Bf + omega(o)*spdiags(ones(N*T,2*T-2),all2all,N*T,N*T);
[Sf,Qf,nb_it_f] = iterated_genlouvain(Bf);
Qf = Qf/twomu;
Sf = reshape(Sf,N,T);

%% 
%%%% checking results

subplot(4,5,o);imagesc(Sf);title(num2str(omega(o)));

%%% to check community changes
for f=1:11
   CommNf(f)=length(unique(Sf(:,f))); 
   CommMaxf(f)=max(unique(Sf(:,f)));
end
%figure;plot(CommN);
OmtempF(o)=omega(o);
QtempF(o)=Qf;
CommNtempF(o,:)=CommNf;
CommMaxtempF(o,:)=CommMaxf;
end


figure;
subplot(1,3,1);plot(OmtempF,QtempF);title('Q'); 
subplot(1,3,2);plot(OmtempF,mean(CommNtempF,2));ylim([0 10]);title('average #Comm'); 
subplot(1,3,3);plot(OmtempF,max(CommMaxtempF,[],2));title('Max Comm'); 


%%% with the thresholded data
%loom=[1 2 3 11 12];
loom=12;
Aft={};

for f=1:length(fish)
   temp=Data_corrMat2.f20.(fish{f}).loomsR{1,loom}(keep,keep);
   temp(isnan(temp))=0;
   temp = threshold_absolute(abs(temp),0.75);  
   Aft{f}=temp;  
end
clear temp


gamma = 1;
omega = [0.1:0.1:2]; %%% this makes a big influence... 

figure;
for o=1:length(omega)

N=length(Af{1});
T=length(Af);
Bf=spalloc(N*T,N*T,(N+T)*N*T);
twomu=0;
for s=1:T
    k=sum(Af{s});
    twom=sum(k);
    twomu=twomu+twom;
    indx=[1:N]+(s-1)*N;
    Bf(indx,indx)=Af{s}-gamma*k'*k/twom;
end
twomu=twomu+T*omega(o)*N*(T-1);
all2all = N*[(-T+1):-1,1:(T-1)];
Bf = Bf + omega(o)*spdiags(ones(N*T,2*T-2),all2all,N*T,N*T);
[Sf,Qf,nb_it_f] = iterated_genlouvain(Bf);
Qf = Qf/twomu;
Sf = reshape(Sf,N,T);

%% 
%%%% checking results

subplot(4,5,o);imagesc(Sf);title(num2str(omega(o)));

%%% to check community changes
for f=1:11
   CommNf(f)=length(unique(Sf(:,f))); 
   CommMaxf(f)=max(unique(Sf(:,f)));
end
%figure;plot(CommN);
OmtempFT(o)=omega(o);
QtempFT(o)=Qf;
CommNtempFT(o,:)=CommNf;
CommMaxtempFT(o,:)=CommMaxf;
end

figure;
subplot(1,3,1);plot(OmtempF,QtempF);ylim([0 1]); hold on; plot(OmtempFT,QtempFT);title('Q'); hold off;
subplot(1,3,2);plot(OmtempF,mean(CommNtempF,2));ylim([0 10]);hold on; plot(OmtempFT,mean(CommNtempFT,2));title('average #Comm'); hold off;
subplot(1,3,3);plot(OmtempF,max(CommMaxtempF,[],2));ylim([0 16]);hold on;plot(OmtempFT,max(CommMaxtempFT,[],2));title('Max Comm'); hold off;
