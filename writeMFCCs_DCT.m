function writeMFCCs_DCT(type)
%   
%   type 1 = extract just visual featuers
%   type 2 = extract visual and audio featuers
%   default value is 1.
%

    
img = imread('trainingImages/image_0001.png');

n = 10;
upRate = 1;

fileList = dir('speechMAT/training/*.mat');
audioList = dir('Data/Training Data wav/*.wav');
for p = 1:length(fileList)
    
    %load speech sequence
    load(['speechMAT/training/' fileList(p).name]);

    coeff_k = [];
    coeff = [];
    for k = 1:length(imgs)
        img = rgb2gray(imgs{k});
        coeff_k = dct2(img);
        %coeff_k = reshape(coeff_k(1:n,1:n),[],1); %100x1 double 
        coeff_k = coeff_k(1:n,1:n);
        coeff(:,k) = coeff_k(:);
    
    end

     
   

    vidreader = VideoReader(['Data/Training Data/' fileList(p).name(1:end-4) '.mov']);
    samPeriod = 1/(upRate * get(vidreader, 'FrameRate')); 
  
    
    if(type == 1)
         coeff = spline(1:size(coeff,2), coeff, linspace(1,size(coeff,2),upRate * size(coeff,2)));
         coeff_inv = coeff';
    
    elseif(type == 2)
        
       % coeff = spline(1:size(coeff,2), coeff, linspace(1,size(coeff,2),upRate * size(coeff,2)));
        
        % adding audio featuers - name of wav file, length of shortframe, 4 =
        % MFCC_E_D_A?14 + 1 + 15 + 15 = 45 components
        audioFeature = featureExtraction(['Data/Training Data wav/' audioList(p).name],20,4);
        % re-cal upRate value for fitting number of rows
        upRate = size(audioFeature,1) / size(coeff,2);
        % up sample
        coeff = spline(1:size(coeff,2), coeff, linspace(1,size(coeff,2),upRate * size(coeff,2)));
        % add visual and audio featuers 
        coeff_inv = [coeff' audioFeature];

%        % re-cal upRate value for fitting number of rows
%         upRate = size(coeff,2) / size(audioFeature,1);
%         % up sample
%         audioFeature_inv = audioFeature';
%         audioFeature_inv = spline(1:size(audioFeature_inv,2), audioFeature_inv, linspace(1,size(audioFeature_inv,2),upRate * size(audioFeature_inv,2)));
%         % add visual and audio featuers 
%         coeff_inv = [coeff'  audioFeature_inv'];
%        
    end
    
    
    fileName = sprintf('ASR/MFCCs/Training/DCT/%s.mfc', fileList(p).name(1:end-4));
    writehtk(fileName,coeff_inv,samPeriod,9);
end


 fileList = dir('speechMAT/testing/*.mat');
 audioList = dir('Data/Testing Data wav/*.wav');
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
        % MFCC_E_D_A?14 + 1 + 15 + 15 = 45 components
        audioFeature = featureExtraction(['Data/Training Data wav/' audioList(p).name],20,4);
        % re-cal upRate value for fitting number of rows
        upRate = size(audioFeature,1) / size(coeff,2);
        % up sample
        coeff = spline(1:size(coeff,2), coeff, linspace(1,size(coeff,2),upRate * size(coeff,2)));
        % add visual and audio featuers 
        coeff_inv = [coeff' audioFeature];

%        % re-cal upRate value for fitting number of rows
%         upRate = size(coeff,2) / size(audioFeature,1);
%         % up sample
%         audioFeature_inv = audioFeature';
%         audioFeature_inv = spline(1:size(audioFeature_inv,2), audioFeature_inv, linspace(1,size(audioFeature_inv,2),upRate * size(audioFeature_inv,2)));
%         % add visual and audio featuers 
%         coeff = [coeff ; audioFeature_inv'];
%         coeff_inv = coeff';
    end
    
    
    fileName = sprintf('ASR/MFCCs/Testing/DCT/%s.mfc', fileList(p).name(1:end-4));
    writehtk(fileName,coeff_inv,samPeriod,9);
end

end