function write_testingMFCCs_NoisyDCT

img = imread('trainingImages/image_0001.png');

type = 2;
n = 10;
upRate = 1;

path = 'Data/Babble Speech/20dB/';

fileList = dir('speechMAT/testing_Noisy/*.mat');
audioList = dir([path '*.wav']);
for p = 1:length(fileList)
    
    %load speech sequence
    load(['speechMAT/testing/' fileList(p).name]);

    coeff_k = [];
    coeff = [];
    for k = 1:length(imgs)
        img = rgb2gray(imgs{k});
        coeff_k = dct2(img);
        %coeff_k = reshape(coeff_k(1:n,1:n),[],1); %100x1 double 
        coeff_k = coeff_k(1:n,1:n);
        coeff(:,k) = coeff_k(:);
    
    end

  
    %coeff = coeff(1:10,1:10);
    vidreader = VideoReader(['Data/Test Data/' fileList(p).name(1:end-4) '.mov']);
    samPeriod = 1/get(vidreader, 'FrameRate'); 
    
    if(type == 1)
         coeff = spline(1:size(coeff,2), coeff, linspace(1,size(coeff,2),upRate * size(coeff,2)));
         coeff_inv = coeff';
    
    elseif(type == 2)
        % adding audio featuers - name of wav file, length of shortframe, 4 =
        % MFCC_E_D_A
        audioFeature = featureExtractionNoisy([path audioList(p).name],20,4);
        % re-cal upRate value for fitting number of rows
        upRate = size(audioFeature,1) / size(coeff,2);
        % up sample
        coeff = spline(1:size(coeff,2), coeff, linspace(1,size(coeff,2),upRate * size(coeff,2)));
        % add visual and audio featuers 
        coeff_inv = [coeff' audioFeature];
    end
    
    
    fileName = sprintf('ASR/MFCCs/Testing/Noisy/%s.mfc', fileList(p).name(1:end-4));
    writehtk(fileName,coeff_inv,samPeriod,9);
end