function write_NoisyPCA

    
type = 2;
quality = 0.95;
upRate = 1;
%adjustment = 0;
%componentsNo = 30;

img = imread('trainingImages/image_0001.png');


load trainingImages/allImages.mat;
[mu,P,v] = pca(allImages,quality);


path = 'Data/Babble Speech/10dB/';

fileList = dir('speechMAT/testing_Noisy/*.mat');
audioList = dir([path '*.wav']);
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
        audioFeature = featureExtractionNoisy([path audioList(n).name],20,4);
        % re-cal upRate value for fitting number of rows
        upRate = size(audioFeature,1) / size(coeff,2);
        % up sample
        coeff = spline(1:size(coeff,2), coeff, linspace(1,size(coeff,2),upRate * size(coeff,2)));
        % add visual and audio featuers 
        coeff_inv = [coeff' audioFeature];
    end
    
    fileName = sprintf('ASR/MFCCs/Testing/Noisy/%s.mfc', fileList(n).name(1:end-4));
    writehtk(fileName,coeff_inv,samPeriod,9);
    
    count = n
end
