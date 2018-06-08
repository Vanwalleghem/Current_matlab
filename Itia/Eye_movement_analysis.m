load('D:\temp\allfisheyes.mat')
idx_ELO=find(sum(ismember(namefish,'ELO'),2)==3);
ZS_eyes=zscore(matrixeyes,1,2);

x = linspace(0.2,size(ZS_eyes,2)/5,size(ZS_eyes,2));
Fighandle=figure;
set(Fighandle, 'Position', [0, 0, 1280, 1024]);
set(findall(Fighandle,'type','text'),'fontSize',12,'fontWeight','bold','FontName','Arial')
xplot=2;yplot=2;
temp_ELO=ZS_eyes(idx_ELO,:);
temp_ERO=ZS_eyes(find(sum(ismember(namefish,'ERO'),2)==3),:);
subplot(xplot,yplot,1);
imagesc(temp_ELO,[-5 5]);
subplot(xplot,yplot,2);
temp=mean(temp_ELO,1);std_temp=std(temp_ELO,1,1);
shadedErrorBar(x, temp, std_temp);axis([0 150 -6 6]);
subplot(xplot,yplot,3);
imagesc(temp_ERO,[-5 5]);
subplot(xplot,yplot,4);
temp=mean(temp_ERO,1);std_temp=std(temp_ERO,1,1);
shadedErrorBar(x, temp, std_temp);axis([0 150 -6 6]);

Fighandle=figure;
set(Fighandle, 'Position', [0, 0, 1280, 1024]);
set(findall(Fighandle,'type','text'),'fontSize',12,'fontWeight','bold','FontName','Arial')
xplot=2;yplot=2;
temp_ELO=matrixeyes(idx_ELO([1 3 4 5 6]),:);
temp_ERO=matrixeyes(find(sum(ismember(namefish,'ERO'),2)==3),:);
subplot(xplot,yplot,1);
imagesc(temp_ELO,[-5 5]);
subplot(xplot,yplot,2);
temp=mean(temp_ELO,1);std_temp=std(temp_ELO,1,1);
shadedErrorBar(x, temp, std_temp);axis([0 150 -6 6]);
subplot(xplot,yplot,3);
imagesc(temp_ERO,[-5 5]);
subplot(xplot,yplot,4);
temp=mean(temp_ERO,1);std_temp=std(temp_ERO,1,1);
shadedErrorBar(x, temp, std_temp);axis([0 150 -6 6]);

%Sup figure for variability
x = linspace(0.2,size(ZS_eyes,2)/5,size(ZS_eyes,2));
Fighandle=figure;
set(Fighandle, 'Position', [100,100, 1200, 200]);
set(findall(Fighandle,'type','text'),'fontSize',12,'fontWeight','bold','FontName','Arial')
xplot=1;yplot=2;ha = tight_subplot(xplot,yplot,[.01 .01],[.01 .01],[.01 .01]);
temp_ELO=ZS_eyes(idx_ELO,:);
temp_ERO=ZS_eyes(find(sum(ismember(namefish,'ERO'),2)==3),:);
axes(ha(1));
temp=mean(temp_ELO,1);std_temp=std(temp_ELO,1,1);
shadedErrorBar(x, temp, std_temp);axis([0 150 -6 6]);set(gca,'XTickLabel',[]);set(gca,'YTickLabel',[]);
axes(ha(2));
temp=mean(temp_ERO,1);std_temp=std(temp_ERO,1,1);
shadedErrorBar(x, temp, std_temp);axis([0 150 -6 6]);set(gca,'XTickLabel',[]);set(gca,'YTickLabel',[]);

figure;
counter=1;
x = linspace(0.25,760/4,760);
order=[5 9 1 12 6 8 3 2 10 11];xplot=length(order);yplot=3;counter2=1;
ha = tight_subplot(xplot,yplot,[.01 .01],[.01 .01],[.01 .01]);
for i=order
    regionName=ItiaList{i};
    if i>10
        idx_temp=LinReg.(regionName).KmeansIdx_merge;
        GoodBet_temp=LinReg.(regionName).GoodBeta_merge;
    else
        idx_temp=LinReg.(regionName).KmeansIdx_select_AVG;
        GoodBet_temp=LinReg.(regionName).GoodBetas_AVG_final;
    end
    counter=counter2;
    if i==5
        counter=counter+1;
    end
    for j=1:length(GoodBet_temp)
        idx_temp2=find(idx_temp==GoodBet_temp(j));
        %subplot(xplot,yplot,counter+(j-1));
        axes(ha(counter+(j-1)));
        meanToPlot=mean(LinReg.(regionName).ZS_rsq(idx_temp2,:),1);
        std_95=std(LinReg.(regionName).ZS_rsq(idx_temp2,:),1,1);%confidence_intervals(LinReg.(regionName).ZS_rsq(idx_temp2,:),95);
        %std(LinReg.(regionName).ZS_rsq(idx_temp2,:),1,1);
        
        H=shadedErrorBar(x(:),meanToPlot-mean(meanToPlot(1:4)), std_95);axis([0 190 -5 5]);
        H.mainLine.Color=colors{i}(j,:);
        H.patch.FaceColor=colors{i}(j,:)/2;
        H.edge(1).Color=colors{i}(j,:)/2;
        H.edge(2).Color=colors{i}(j,:)/2;
        set(gca,'XTickLabel',[]);set(gca,'YTickLabel',[]);
    end
    counter2=counter2+yplot;
end
clearvars i j k counter counter2 H idx_temp GoodBet_temp std_95 meanToPlot

EyeMvmt_intensity={};
EyeMvmt_timing={};
for idx_fish=1:size(matrixeyes,1)
    if ismember(idx_fish,idx_ELO)
        [EyeMvmt_intensity{idx_fish},EyeMvmt_timing{idx_fish}] = findpeaks(ZS_eyes(idx_fish,:),'MinPeakDistance',20,'MinPeakHeight',1);
    else
        [EyeMvmt_intensity{idx_fish},EyeMvmt_timing{idx_fish}] = findpeaks(1.01*max(ZS_eyes(idx_fish,:))-ZS_eyes(idx_fish,:),'MinPeakDistance',20,'MinPeakHeight',1);
    end    
end

figure;
for idx_fish=1:size(matrixeyes,1)
    plot(ZS_eyes(idx_fish,:));title(num2str(idx_fish));
    pause;
end

i=1;EyeMvmt_timing{i,1}=[43 83 122 162 203 243 281 322 362 401 441 480 521 561 601];
i=2;EyeMvmt_timing{i,1}=[43 82 122 162 201 242 281 321 361 401 441 481 522 568 626];
i=2;EyeMvmt_timing{i,2}=[349 542 597 680];
i=3;EyeMvmt_timing{i,1}=[41 82 121 145 161 201 241 272 281 321 362 401 442 482 521 562 601];
i=4;EyeMvmt_timing{i,1}=[41 83 94 121 161 202 244 282 323 361 401 442 481];
i=5;EyeMvmt_timing{i,1}=[43 81 122 143 161 201 221 241 281 303 321 361 444 598 660];
i=5;EyeMvmt_timing{i,2}=[351];
i=6;EyeMvmt_timing{i,1}=[43 84 122 162 203 242 282 322 362 405 441 482 542];
i=6;EyeMvmt_timing{i,2}=[94 297 348 494];
i=7;EyeMvmt_timing{i,1}=[42 82 122 162 177 205 222 245 284 322 362 407 442 574 617 653 682 723];
i=7;EyeMvmt_timing{i,2}=[58 748 753];
i=8;EyeMvmt_timing{i,1}=[43 83 123 167 202 242 282 322 373 563];
i=9;EyeMvmt_timing{i,2}=[97 378 448 673];
i=9;EyeMvmt_timing{i,1}=[42 82 122 129 202 241 281 329 359 374 389 403 477 575 661 672 702 733];
i=10;EyeMvmt_timing{i,2}=[138 168 233 257 344 391 407 473 488 548];
i=10;EyeMvmt_timing{i,1}=[42 84 122 162 202 241 283 293 321 341 364 442 536 651];
i=11;EyeMvmt_timing{i,1}=[43 83 94 122 161 189 202 245 281 322 361];
i=11;EyeMvmt_timing{i,2}=[142 170 213 234 250 275 306 339 376 487 570 590 644];
i=12;EyeMvmt_timing{i,1}=[43 82 122 165 205 243 283 324 346 363 403 442 465 481 542 662 671];
i=12;EyeMvmt_timing{i,2}=[210 333 627 717];
i=13;EyeMvmt_timing{i,1}=[41 89 124 205 243 324 343 364 402 444 481];
i=13;EyeMvmt_timing{i,2}=[132 177 231 254 294 431 531 570];

Stimuli=zeros(6,size(ZS2,2));
EyeMvmt_stimulus={};
GCaMP6=[0,1.69644104899772,5.13796058542217,8.27886020152244,10.3756715204800,11.8173714529814,12.2425184714093,10.8571417354877,8.80831829681196,6.91339112244670,5.46959264663869,4.30868766622567,3.42533619066766,2.75378443486879,2.18017250852183,1.72816235135824,1.32732537295463,1.00684435500268,0.730210038304555,0.530242444093118,0.362253250339685,0.227668255288566,0.0869242416152502,0.000718266708050853,-0.0828334873368325]';
for idx_fish=1:size(matrixeyes,1)
    Stimuli=zeros(2,size(ZS2,2));
    GCaMP6;
    for idx_mvmt=1:length(EyeMvmt_timing{idx_fish,1});
        Stimuli(1,EyeMvmt_timing{idx_fish,1}(idx_mvmt):EyeMvmt_timing{idx_fish,1}(idx_mvmt)+size(GCaMP6,1)-1)=GCaMP6;
    end
    if ~isempty(EyeMvmt_timing{idx_fish,2})
        for idx_mvmt=1:length(EyeMvmt_timing{idx_fish,2});
            Stimuli(2,EyeMvmt_timing{idx_fish,2}(idx_mvmt):EyeMvmt_timing{idx_fish,2}(idx_mvmt)+size(GCaMP6,1)-1)=GCaMP6;
        end
    else
        Stimuli=Stimuli(1,:);
    end
    EyeMvmt_stimulus{idx_fish}=Stimuli;
end

for idx_fish=1:length(PerBrainRegions)
    Stimuli=EyeMvmt_stimulus{idx_fish};Stimuli=Stimuli(:,1:760);
    for region_nb=1:length(ItiaList)
        region_name=ItiaList{region_nb};
        ZS_temp=PerBrainRegions(idx_fish).(region_name).ZS;
        coefficients=struct();
        rsquared=zeros(1,length(ZS_temp));
        parfor i=1:size(ZS_temp,1)
            mdl=fitlm(Stimuli',ZS_temp(i,:));%,'interactions');
            coefficients(i).coef=mdl.Coefficients;
            rsquared(i)=mdl.Rsquared.Adjusted;
        end
        MovReg(idx_fish).(region_name).coef=coefficients;
        MovReg(idx_fish).(region_name).rsquared=rsquared;
    end
end

ZS_move=struct();
ZS_pool=struct();
for idx_fish=1:length(PerBrainRegions)
    for region_nb=1:length(ItiaList)
        region_name=ItiaList{region_nb};
        idx_rsq=find(MovReg(idx_fish).(region_name).rsquared>0.1);
        ZS_move(idx_fish).(region_name)=PerBrainRegions(idx_fish).(region_name).ZS(idx_rsq,:);
        if idx_fish==1
            ZS_pool.(region_name)=ZS_move(idx_fish).(region_name);
        else
            ZS_pool.(region_name)=[ZS_pool.(region_name) ; ZS_move(idx_fish).(region_name)];        
        end
    end
end

%May be good then do look at eye movement + trap and trap only responses in
%the brain regions. But only if reviewers ask it


for idx_fish=1:length(PerBrainRegions)
    region_name=ItiaList{2};
    ZS_temp=PerBrainRegions(idx_fish).(region_name).ZS;
    coefficients=struct();
    rsquared=zeros(1,length(ZS_temp));
    parfor i=1:size(ZS_temp,1)
        mdl=fitlm(Stimuli',ZS_temp(i,:));%,'interactions');
        coefficients(i).coef=mdl.Coefficients;
        rsquared(i)=mdl.Rsquared.Adjusted;
    end
    RegTrap(idx_fish).(region_name).coef=coefficients;
    RegTrap(idx_fish).(region_name).rsquared=rsquared;
end

ZS_cereb={};
for idx_fish=1:length(PerBrainRegions)
    region_name=ItiaList{2};    
    idx_rsq=find(RegTrap(idx_fish).(region_name).rsquared>0.1);
    ZS_cereb{idx_fish}=PerBrainRegions(idx_fish).(region_name).ZS(idx_rsq,:);
    if idx_fish==1
        ZS_pool_cereb=ZS_cereb{idx_fish};
    else
        ZS_pool_cereb=[ZS_pool_cereb ; ZS_cereb{idx_fish}];
    end    
end

idx_fish=4;region_name=ItiaList{2};  
idx_rsq_mov=find(MovReg(idx_fish).(region_name).rsquared>0.1);
idx_rsq_trap=find(RegTrap(idx_fish).(region_name).rsquared>0.1);
idx_trap_only=ismember(idx_rsq_trap,idx_rsq_mov);
figure;imagesc(PerBrainRegions(idx_fish).(region_name).ZS(idx_rsq_trap(~idx_trap_only),:),[-1 3]);

