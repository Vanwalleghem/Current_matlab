%% Dantrolene Pref 1
%Composite2_MCcrop
%First ROI : makeOval(137, 83, 46, 46);
window_baseline=21;
window_smooth=7;

Pref1=[];
Pref1b=[];
Pref1c=[];

DF=DeltaF2(Pref1(:,2:3)',window_baseline,window_smooth);
DF_norm=DF;
for i=1:size(DF,1)
    DF_norm(i,:)=DF(i,:)/max(abs(DF(i,:)));
end
%figure;plot(DF_norm');hold on;plot(Pref1(:,3));

%Second ROI : makeOval(172, 246, 46, 46);
DFb=DeltaF2(Pref1b(:,2:3)',window_baseline,window_smooth);
DF_normb=DFb;
for i=1:size(DFb,1)
    DF_normb(i,:)=DFb(i,:)/max(abs(DFb(i,:)));
end
%figure;plot(DF_normb');hold on;plot(Pref1b(:,4));

%Third ROI : makeOval(248, 307, 62, 65);
DFc=DeltaF2(Pref1c(:,2:3)',window_baseline,window_smooth);
DF_normc=DFc;
for i=1:size(DFc,1)
    DF_normc(i,:)=DFc(i,:)/max(abs(DFc(i,:)));
end
%figure;plot(DF_normc');hold on;plot(Pref1c(:,4));

Fighandle=figure;
set(Fighandle, 'Position', [10, 10, 1000, 500]);
ha=tight_subplot(3,1);
axes(ha(1));
plot(DF_norm');hold on;plot(Pref1(:,4));
axes(ha(2));
plot(DF_normb');hold on;plot(Pref1b(:,4));
axes(ha(3));
plot(DF_normc');hold on;plot(Pref1c(:,4));

ProcessedData=struct();
ProcessedData(1).Raw{1}=Pref1;
ProcessedData(1).Raw{2}=Pref1b;
ProcessedData(1).Raw{3}=Pref1c;
ProcessedData(1).DF{1}=DF;
ProcessedData(1).DF{2}=DFb;
ProcessedData(1).DF{3}=DFc;

%% Dantrolene Pref 5
%Composite2_MCcrop
%First ROI : makeOval(127, 256, 55, 64);
window_baseline=21;
window_smooth=7;

Pref1=[];
Pref1b=[];
Pref1c=[];

DF=DeltaF2(Pref1(:,2:3)',window_baseline,window_smooth);
DF_norm=DF;
for i=1:size(DF,1)
    DF_norm(i,:)=DF(i,:)/max(abs(DF(i,:)));
end
%figure;plot(DF_norm');hold on;plot(Pref1(:,3));

%Second ROI : makeOval(36, 222, 95, 52);
DFb=DeltaF2(Pref1b(:,2:3)',window_baseline,window_smooth);
DF_normb=DFb;
for i=1:size(DFb,1)
    DF_normb(i,:)=DFb(i,:)/max(abs(DFb(i,:)));
end
%figure;plot(DF_normb');hold on;plot(Pref1b(:,4));

%Third ROI : makeOval(209, 43, 55, 53);
DFc=DeltaF2(Pref1c(:,2:3)',window_baseline,window_smooth);
DF_normc=DFc;
for i=1:size(DFc,1)
    DF_normc(i,:)=DFc(i,:)/max(abs(DFc(i,:)));
end
%figure;plot(DF_normc');hold on;plot(Pref1c(:,4));

Fighandle=figure;
set(Fighandle, 'Position', [10, 10, 1000, 500]);
ha=tight_subplot(3,1);
axes(ha(1));
plot(DF_norm');hold on;plot(Pref1(:,4));
axes(ha(2));
plot(DF_normb');hold on;plot(Pref1b(:,4));
axes(ha(3));
plot(DF_normc');hold on;plot(Pref1c(:,4));

ProcessedData=struct();
ProcessedData(2).Raw{1}=Pref1;
ProcessedData(2).Raw{2}=Pref1b;
ProcessedData(2).Raw{3}=Pref1c;
ProcessedData(2).DF{1}=DF;
ProcessedData(2).DF{2}=DFb;
ProcessedData(2).DF{3}=DFc;

