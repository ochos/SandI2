%pca

dlist = dir('trainingImages/*.jpg');
allImages = [];

for k = 1:length(dlist)
    img = imread(['trainingImages/' dlist(k).name]);
%     img = rgb2gray(img);
%     allImages = cat(2,allImages,double(img(:)));
    
    % colour  
    allImages = cat(2,allImages,double(img(:)));
    
end

%mu = mean, P = modes of variation, v = most variestion
[mu,P,v] = pca(allImages,0.95);


%gray
% muImage = reshape(mu,size(img));
% imshow(reshape(uint8(mu),size(img));

%colour
imshow(reshape(uint8(mu),size(img)));



mimg = P(:,1);
mimg = mimg - min(mimg);
mimg = 255 * (mimg / max(mimg));
imshow(reshape(uint8(mimg),size(img)));




b = P' * (double(img(:)) - mu);

%gray
% img3 = mu + P * b;
% img3 = uint8(reshape(img3,size(img)));
% imshow(img3);
% imagesc(reshape(P(:,1),size(img)));
% imagesc(reshape(P(:,2),size(img)));
% imagesc(reshape(P(:,3),size(img)));

%colour
imgr = mu + P * b;
imr = reshape(uint8(imgr),size(img));
imshow(imr);


% fileName = sprintf('ASR/MFCCs/Training/DCT/%s.mfc', fileList(p).name(1:end-4));
% writehtk(fileName,coeff_inv,samPeriod,9);
    


%Inspecting a shape model