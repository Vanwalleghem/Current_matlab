figure;
plot(Coeff(1,:),'r');hold on;
plot(Coeff(2,:),'g');
plot(Coeff(3,:),'b');hold off;

figure;
for i=1:size(Scores,1)
   image(squeeze(Scores(i,:,:,:)));
   pause; 
end

for i=1:3    
    outputFileName = strcat('Ca8_loom_bin2_Thunder_',num2str(i),'.tif');
    min_score=min(min(min(Scores(:,:,:,i))));
    max_score=max(max(max(Scores(:,:,:,i))));
    for Plane=1:size(Scores,1)
        image=squeeze(Scores(Plane, :, :,i));
        image=(image-min_score)/(max_score-min_score);image=image*256;image=uint16(image);
        imwrite(image, outputFileName, 'WriteMode', 'append');
    end
end

imData=bigread2('CA8_AllButFlow_F2_range250_step5_exposure20_power40_1.tif');
dims=size(Betas); %%%to make dimensions out of the Beta values
temp=reshape(Betas,dims(1)*dims(2)*dims(3),dims(4)); %%%% to create a 3D variable instead of 4D
rsq_tmp=reshape(rsquared,1,dims(1)*dims(2)*dims(3)); %%%%Dont remeber what this is for...
idx_rsq=find(rsq_tmp>0.1); %%% to identify the rsquare values above 0.1
max_betas=max(temp,[],1); %%%%to lo for the maximum values of the Betas
mean_betas=mean(temp,1);  %%%%to lo for the mean values of the Betas
std_betas=std(temp,1,1);  %%%%to lo for the standar deviation values of the Betas
outputFileName = strcat('Ca8_AllButFlow_Thunder_Fish1.tif');
min_Beta1=min(min(min(Betas(:,:,:,2))));
max_Beta1=max(max(max(Betas(:,:,:,2))));
min_Beta2=min(min(min(Betas(:,:,:,3))));
max_Beta2=max(max(max(Betas(:,:,:,3))));

outputFileName=('Itia_test.tif');
for Plane=1:size(rsquared,1)
    image_avg=double(imData(:,:,Plane));image_avg=image_avg/max(max(image_avg));image_avg=image_avg*128;
    image_avg=repmat(image_avg,1,1,3);image_avg=uint8(image_avg);
    image_1=squeeze(Betas(Plane, :, :,2));
    image_1(image_1<(mean_betas(2)+5*std_betas(2)))=0;
    image_1=image_1/max_Beta1;image_1=image_1*128;image_1=uint8(image_1);
    image_avg(:,:,1)=image_avg(:,:,1)+image_1;
    image_2=squeeze(Betas(Plane, :, :,3));
    image_2(image_2<(mean_betas(3)+5*std_betas(3)))=0;
    image_2=image_2/max_Beta2;image_2=image_2*128;image_2=uint8(image_2);
    image_avg(:,:,2)=image_avg(:,:,2)+image_2;
    image_3=squeeze(rsquared(Plane, :, :));
    image_3(image_3<0.05)=0;
    image_3=image_3/max(max(image_3));image_3=image_3*128;image_3=uint8(image_3);
    image_avg(:,:,3)=image_avg(:,:,3)+image_3;
    imwrite(image_avg, outputFileName, 'WriteMode', 'append');
end

figure;
for Plane=1:size(rsquared,1)
    subplot(1,2,1);imagesc(squeeze(rsquared(Plane,:,:)),[0.05 0.2]);colormap jet;
    subplot(1,2,2);imagesc(squeeze(rsq2(Plane,:,:)),[0.05 0.2]);colormap jet;
    pause;
end

outputFileName=('Itia_test.tif');
for Plane=1:size(rsquared,1)
    image_3=squeeze(rsquared(Plane, :, :));
    image_3(image_3<0.05)=0;
    image_3=image_3*128;image_3=uint8(image_3);    
    imwrite(image_3, outputFileName, 'WriteMode', 'append');
end
