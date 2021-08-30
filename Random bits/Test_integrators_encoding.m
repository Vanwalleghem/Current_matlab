int_quant=nan(length(Fish_list),9);
Int_fish_temp=idx_Fish(Test_Int_idx);
good_cluster=2;
delay=10;
Nb_neurons=5;
for fish_nb=1:length(Fish_list)
    idx_fish_temp=find(Int_fish_temp==Fish_list(fish_nb));        
    idx_fish_temp=intersect(find(Int_idxKmeans==good_cluster),idx_fish_temp);
    if length(idx_fish_temp)>Nb_neurons                
        Int_data_temp_mean=mean(Int_data_temp(idx_fish_temp,:),1);        
        int_quant(fish_nb,1)=trapz(Int_data_temp_mean(49-delay:176)-min(Int_data_temp_mean(49:176)));
        int_quant(fish_nb,2)=trapz(Int_data_temp_mean(249-delay:440)-min(Int_data_temp_mean(249:440)));
        int_quant(fish_nb,3)=trapz(Int_data_temp_mean(549-delay:800)-min(Int_data_temp_mean(549:800)));
        int_quant(fish_nb,4)=trapz(Int_data_temp_mean(1000-delay:1120)-min(Int_data_temp_mean(1000:1120)));
        int_quant(fish_nb,9)=trapz(Int_data_temp_mean(1460-delay:1570)-min(Int_data_temp_mean(1460:1570)));
        int_quant(fish_nb,5)=trapz(Int_data_temp_mean(1120-delay:1272)-min(Int_data_temp_mean(1120:1272)));
        int_quant(fish_nb,6)=trapz(Int_data_temp_mean(1825-delay:2005)-min(Int_data_temp_mean(1825:2005)));
        int_quant(fish_nb,7)=trapz(Int_data_temp_mean(2100-delay:2245)-min(Int_data_temp_mean(2100:2245)));
        int_quant(fish_nb,8)=length(idx_fish_temp>Nb_neurons);
    end
end

delay=10;
Area_flow=[];
Area_flow(1)=trapz(Speed_flow(1,49-delay:176));
Area_flow(5)=trapz(Speed_flow(1,1100-delay:1272));
Area_flow(2)=trapz(Speed_flow(1,249-delay:440));
Area_flow(3)=trapz(Speed_flow(1,549-delay:800));
Area_flow(4)=trapz(Speed_flow(1,1008-delay:1100));
Area_flow(6)=trapz(Speed_flow(1,1825-delay:2005));
Area_flow(7)=trapz(Speed_flow(1,2080-delay:2245));

Area_flow=[];
Area_flow(1)=trapz(Speed_flow(1,49-delay:176));
Area_flow(5)=trapz(Speed_flow(1,1100-delay:1272));
Area_flow(2)=trapz(Speed_flow(1,249-delay:440));
Area_flow(3)=trapz(Speed_flow(1,549-delay:800));
Area_flow(4)=trapz(Speed_flow(1,1008-delay:1100));
Area_flow(6)=trapz(Speed_flow(1,1825-delay:2005));
Area_flow(7)=trapz(Speed_flow(1,2080-delay:2245));

Sp_flow=[];
Sp_flow(1)=max(Speed_flow(1,49-delay:176));
Sp_flow(5)=max(Speed_flow(1,1100-delay:1272));
Sp_flow(2)=max(Speed_flow(1,249-delay:440));
Sp_flow(3)=max(Speed_flow(1,549-delay:800));
Sp_flow(4)=max(Speed_flow(1,1008-delay:1100));
Sp_flow(6)=max(Speed_flow(1,1825-delay:2005));
Sp_flow(7)=max(Speed_flow(1,2080-delay:2245));

Time_flow=[];
Time_flow(1)=sum(Speed_flow(1,49-delay:176)>0);
Time_flow(5)=sum(Speed_flow(1,1100-delay:1272)>0);
Time_flow(2)=sum(Speed_flow(1,249-delay:440)>0);
Time_flow(3)=sum(Speed_flow(1,549-delay:800)>0);
Time_flow(4)=sum(Speed_flow(1,1008-delay:1100)>0);
Time_flow(6)=sum(Speed_flow(1,1825-delay:2005)>0);
Time_flow(7)=sum(Speed_flow(1,2080-delay:2245)>0);

Coding_Volume=[];
Coding_Volume=int_quant(:,1:7);
Coding_Volume(19,:)=Area_flow;
[Coding_Volume, idx_sort]=sortrows(Coding_Volume',19);
Coding_Volume=single(round(Coding_Volume,2));

figure;
for i=1:7
    scatter(squeeze(int_quant(:,1)),Area_flow(1));
end

ZS=Int_data_temp(Int_idxKmeans==good_cluster,:);
save('Integrators_ZS.mat','ZS');

Spikes=[];
for i=715:size(ZS,1)
    [~,~,~,~,~,Spikes(i,:)]=constrained_foopsi(ZS(i,:));
end
save('Integrators_Spikes.mat','Spikes');

Features=zeros(4,length(Speed_flow));
Features(1,6:end)=Speed_flow(1,1:length(Speed_flow)-5)==1;
Features(2,6:end)=Speed_flow(1,1:length(Speed_flow)-5)==2;

temp=zeros(1,length(Speed_flow));
back=    [56 256 557 1006 1106 1466 1827 2086]+5; %Infuse
back_off=[106 356 706 1056 1206 1516 1926 2176]+5;
for i=1:length(back)
    length_temp=back_off(i)-back(i);
    temp(1,back(i):back_off(i))=[1:1:length_temp+1];
end

Features(3,:)=temp;

temp=zeros(1,length(Speed_flow));
for i=1:length(back)
    length_temp=back_off(i)-back(i);
    temp(1,back(i):back_off(i))=cumsum(Speed_flow(1,back(i):back_off(i)));
end

Features(4,:)=temp;
save('Integrators_Features.mat','Features');