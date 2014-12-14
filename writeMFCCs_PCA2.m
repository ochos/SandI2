function writeMFCCs_PCA2

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
    
type = 2;
quality = 0.95;
upRate = 1;
%adjustment = 0;
%componentsNo = 30;

img = imread('trainingImages/image_0001.png');

load trainingImages/allImages.mat;
[mu,P,v] = pca(allImages,quality);

fileList = dir('speechMAT/training/*.mat');
audioList = dir('Data/Training Data wav/*.wav');
for n = 1:length(fileList)
    
   
    load (['speechMAT/training/' fileList(n).name]);

    coeff = [];
    for k = 1:size(imgs,2)
        img = imgs{k};
        b = P' * (double(img(:)) - mu);
        coeff = cat(2,coeff,double(b(:)));
    end

    vidreader = VideoReader(['Data/Training Data/' fileList(n).name(1:end-4) '.mov']);
    samPeriod = 1/(upRate * get(vidreader, 'FrameRate')); 
  
    
    if(type == 1)
         
         coeff = spline(1:size(coeff,2), coeff, linspace(1,size(coeff,2),upRate * size(coeff,2)));
         coeff_inv = coeff';
    
    elseif(type == 2)
        
       
        % adding audio featuers - name of wav file, length of shortframe, 4 =
        % MFCC_E_D_A?14 + 1 + 15 + 15 = 45 components
        audioFeature = featureExtraction(['Data/Training Data wav/' audioList(n).name],20,4);
        % re-cal upRate value for fitting number of rows
        upRate = size(audioFeature,1) / size(coeff,2);
        % up sample
        coeff = spline(1:size(coeff,2), coeff, linspace(1,size(coeff,2),upRate * size(coeff,2)));
        % add visual and audio featuers 
        coeff_inv = [coeff' audioFeature];
    end
    
    fileName = sprintf('ASR/MFCCs/Training/PCA/%s.mfc', fileList(n).name(1:end-4));
    writehtk(fileName,coeff_inv,samPeriod,9);
 
end

fileList = dir('speechMAT/testing/*.mat');
audioList = dir('Data/Testing Data wav/*.wav');
for n = 1:length(fileList)
    
   
    load (['speechMAT/testing/' fileList(n).name]);

    coeff = [];
    for k = 1:size(imgs,2)
        img = imgs{k};
        b = P' * (double(img(:)) - mu);
        coeff = cat(2,coeff,double(b(:)));
    end

    vidreader = VideoReader(['Data/Test Data/' fileList(n).name(1:end-4) '.mov']);
    samPeriod = 1/(upRate * get(vidreader, 'FrameRate')); 
  
    
    if(type == 1)
         
         coeff = spline(1:size(coeff,2), coeff, linspace(1,size(coeff,2),upRate * size(coeff,2)));
         coeff_inv = coeff';
    
    elseif(type == 2)
        
       
        % adding audio featuers - name of wav file, length of shortframe, 4 =
        % MFCC_E_D_A?14 + 1 + 15 + 15 = 45 components
        audioFeature = featureExtraction(['Data/Testing Data wav/' audioList(n).name],20,4);
        % re-cal upRate value for fitting number of rows
        upRate = size(audioFeature,1) / size(coeff,2);
        % up sample
        coeff = spline(1:size(coeff,2), coeff, linspace(1,size(coeff,2),upRate * size(coeff,2)));
        % add visual and audio featuers 
        coeff_inv = [coeff' audioFeature];
    end
    
    fileName = sprintf('ASR/MFCCs/Testing/PCA/%s.mfc', fileList(n).name(1:end-4));
    writehtk(fileName,coeff_inv,samPeriod,9);
    
    count = n
end


end