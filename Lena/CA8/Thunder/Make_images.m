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
    for K=1:size(Scores,1)
        image=squeeze(Scores(K, :, :,i));
        image=(image-min_score)/(max_score-min_score);image=image*256;image=uint16(image);
        imwrite(image, outputFileName, 'WriteMode', 'append');
    end
end