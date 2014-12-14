function writeMFCCs_PCA

    %write down each speeches images for training and testing
%   fileList = dir('speechMAT/training/*.mat');  
%   for p = 1:length(fileList)  
%     imageCount = 1;
%     load(['speechMAT/training/' fileList(p).name]);
%     numberOfImages = size(imgs,2);
%     idx = randperm(numberOfImages);
%     for k = 1:numberOfImages
%         img = imgs{k};
%         imwrite(img,sprintf(['images/Training/'  fileList(p).name(1:end-4) '/image_%.2d.jpg'],imageCount),'jpg','Quality',100);
%         imageCount = imageCount + 1;
%     end
%     
%   end
%   fileList = dir('speechMAT/testing/*.mat');  
%   for p = 1:length(fileList)  
%     imageCount = 1;
%     load(['speechMAT/testing/' fileList(p).name]);
%     numberOfImages = size(imgs,2);
%     idx = randperm(numberOfImages);
%     for k = 1:numberOfImages
%         img = imgs{k};
%         imwrite(img,sprintf(['images/Testing/'  fileList(p).name(1:end-4) '/image_%.2d.jpg'],imageCount),'jpg','Quality',100);
%         imageCount = imageCount + 1;
%     end
%     
%   end
    
type = 1;
quality = 0.95;
upRate = 1;
adjustment = 0;
componentsNo = 50;

img = imread('trainingImages/image_01.jpg');

fileList = dir('speechMAT/training/*.mat');
audioList = dir('Data/Training Data wav/*.wav');
for p = 1:length(fileList)
    
   
    imageList = dir(['images/Training/' fileList(p).name(1:end-4) '/*.jpg']);
    allImages = [];

    for k = 1:length(imageList);
        img = imread(['images/Training/' fileList(p).name(1:end-4) '/' imageList(k).name]);
    %     img = rgb2gray(img);
    %     allImages = cat(2,allImages,double(img(:)));
    
        % colour 
        allImages = cat(2,allImages,double(img(:)));
    end

    %mu = mean, P = modes of variation, v = most variestion
    [mu,P,v] = pca(allImages,quality);

    %gray
%     muImage = reshape(mu,size(img));
%     imshow(reshape(uint8(mu),size(img));

    %colour
%     imshow(reshape(uint8(mu),size(img)));

%     mimg = P(:,1);
%     mimg = mimg - min(mimg);
%     mimg = 255 * (mimg / max(mimg));
%     imshow(reshape(uint8(mimg),size(img)));


    coeff = [];
    for k = 1:length(imageList)
        img = imread(['images/Training/' fileList(p).name(1:end-4) '/' imageList(k).name]);
        b = P' * (double(img(:)) - mu);
        coeff = cat(2,coeff,double(b(:)));
    end

   

%gray
% img3 = mu + P * b;
% img3 = uint8(reshape(img3,size(img)));
% imshow(img3);
% imagesc(reshape(P(:,1),size(img)));
% imagesc(reshape(P(:,2),size(img)));
% imagesc(reshape(P(:,3),size(img)));

%colour
% imgr = mu + P * b;
% imr = reshape(uint8(imgr),size(img));
% imshow(imr);
    

    vidreader = VideoReader(['Data/Training Data/' fileList(p).name(1:end-4) '.mov']);
    samPeriod = 1/(upRate * get(vidreader, 'FrameRate')); 
  
    
    if(type == 1)
        
         %adjustment = componentsNo / size(coeff,2);
        % coefff = spline(1:size(coeff,2), coeff, linspace(1,size(coeff,1),adjustment * size(coeff,1)));
         
         coeff = spline(1:size(coeff,2), coeff, linspace(1,size(coeff,2),upRate * size(coeff,2)));
         coeff_inv = coeff';
    
    elseif(type == 2)
        % adding audio featuers - name of wav file, length of shortframe, 4 =
        % MFCC_E_D_A?14 + 1 + 15 + 15 = 45 components
        audioFeature = featureExtraction(['Data/Training Data wav/' audioList(p).name],20,4);
        % re-cal upRate value for fitting number of rows
        upRate = size(audioFeature,1) / size(coeff,2);
        % up sample
        coeff = spline(1:size(coeff,2), coeff, linspace(1,size(coeff,2),upRate * size(coeff,2)));
        % add visual and audio featuers 
        coeff_inv = [coeff' audioFeature];
    end
    
    
    fileName = sprintf('ASR/MFCCs/Training/PCA/%s.mfc', fileList(p).name(1:end-4));
    writehtk(fileName,coeff_inv,samPeriod,9);
end

fileList = dir('speechMAT/testing/*.mat');
audioList = dir('Data/Testing Data wav/*.wav');
for p = 1:length(fileList)
    
   
    imageList = dir(['images/Testing/' fileList(p).name(1:end-4) '/*.jpg']);
    allImages = [];

    for k = 1:length(imageList);
        img = imread(['images/Testing/' fileList(p).name(1:end-4) '/' imageList(k).name]);
    %     img = rgb2gray(img);
    %     allImages = cat(2,allImages,double(img(:)));
    
        % colour 
        allImages = cat(2,allImages,double(img(:)));
    end

    %mu = mean, P = modes of variation, v = most variestion
    [mu,P,v] = pca(allImages,quality);

    %gray
%     muImage = reshape(mu,size(img));
%     imshow(reshape(uint8(mu),size(img));

    %colour
%     imshow(reshape(uint8(mu),size(img)));

%     mimg = P(:,1);
%     mimg = mimg - min(mimg);
%     mimg = 255 * (mimg / max(mimg));
%     imshow(reshape(uint8(mimg),size(img)));


    coeff = [];
    for k = 1:length(imageList)
        img = imread(['images/Testing/' fileList(p).name(1:end-4) '/' imageList(k).name]);
        b = P' * (double(img(:)) - mu);
        coeff = cat(2,coeff,double(b(:)));
    end

   

%gray
% img3 = mu + P * b;
% img3 = uint8(reshape(img3,size(img)));
% imshow(img3);
% imagesc(reshape(P(:,1),size(img)));
% imagesc(reshape(P(:,2),size(img)));
% imagesc(reshape(P(:,3),size(img)));

%colour
% imgr = mu + P * b;
% imr = reshape(uint8(imgr),size(img));
% imshow(imr);
    

    vidreader = VideoReader(['Data/Test Data/' fileList(p).name(1:end-4) '.mov']);
    samPeriod = 1/(upRate * get(vidreader, 'FrameRate')); 
  
    
    if(type == 1)
         adjustment = componentsNo / size(coeff,2);
         coeff = spline(1:size(coeff,2), coeff, linspace(1,adjustment * size(coeff,2),upRate * size(coeff,2)));
         coeff_inv = coeff';
    
    elseif(type == 2)
        % adding audio featuers - name of wav file, length of shortframe, 4 =
        % MFCC_E_D_A?14 + 1 + 15 + 15 = 45 components
        audioFeature = featureExtraction(['Data/Testing Data wav/' audioList(p).name],20,4);
        % re-cal upRate value for fitting number of rows
        upRate = size(audioFeature,1) / size(coeff,2);
        % up sample
        coeff = spline(1:size(coeff,2), coeff, linspace(1,size(coeff,2),upRate * size(coeff,2)));
        % add visual and audio featuers 
        coeff_inv = [coeff' audioFeature];
    end
    
    
    fileName = sprintf('ASR/MFCCs/Testing/PCA/%s.mfc', fileList(p).name(1:end-4));
    writehtk(fileName,coeff_inv,samPeriod,9);
end


end