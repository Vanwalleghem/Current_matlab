NewFlow=zeros(4,size(ZS,2));
fwd=[0 40 100 190 210 282 374 404]; %Infuse
fwd_off=[10 60 130 200 230 292 384 424];
back=[20 70 140 160 240 252 294 314 394 434]; %Withdraw
back_off=[30 90 150 180 250 262 304 344 404 444];
framerate=2.1;start=50;
back= round(start*framerate+ back*framerate); %Infuse
back_off=round(start*framerate+ back_off*framerate);
fwd=  round(start*framerate+ fwd*framerate); %Withdraw
fwd_off= round(start*framerate+ fwd_off*framerate);
for i=1:length(back)
NewFlow(1,back(i):back(i)+size(spike,1)-1)=spike';
NewFlow(2,back(i):back_off(i))=1;
end
for i=1:length(fwd)
NewFlow(3,fwd(i):fwd(i)+size(spike,1)-1)=spike';
NewFlow(4,fwd(i):fwd_off(i))=1;
end
clearvars GCaMP6 back back_off fwd fwd_off;


figure;plot(x,Cmap_ZS(26,:));hold on;plot(x,sum(NewFlow([1 3],:),1));

NewFlow(5,:)=[0:size(ZS,2):1];
figure;plot(NewFlow([1 3],:)');ylim([0 2])


figure;plot(Cmap_ZS(3,:))

figure;plot(Cmap_ZS(13,:))
hold on;plot(NewFlow(1:4,:)');

[Int_model,Int_GoodBetas]=Test_Regress(Cmap_ZS(:,1:400),NewFlow(:,1:400),idxKmeans_ZS,0.1);

fwd=[0 40 100 190 210 282 374 404]; %Infuse
fwd_off=[10 60 130 200 230 292 384 424];
back=[20 70 140 160 240 252 294 314 394 434]; %Withdraw
back_off=[30 90 150 180 250 262 304 344 404 444];

StimLength=1200;
x = linspace(0,StimLength/framerate,StimLength);
figure;plot(x,Cmap_ZS(26,:));hold on;plot(x,sum(NewFlow([1 3],:),1));

figure;plot(x,Cmap_ZS(4,:));hold on;plot(x,sum(NewFlow([1 3],:),1));