function pre_process()

fileList = dir('Data/Training Data/*.mov');
startW = 'speechMAT/training/speech';
%startM = 'ASR/MFCCs/Training/speech';
number1to9 = '00';
number10to20 = '0';
extensionMovie = '.mov';
extensionMat = '.mat';


for f = 1:length(fileList)
    
    if n >= 1 && n <=9
        matname = strcat(startW, number1to9, int2str(n), extensionMat);
       % mfcName = strcat(startM, number1to9, int2str(n), extensionM);
    else
        matname = strcat(startW, number10to20, int2str(n), extensionMat);
       % mfcName = strcat(startM, number10to20, int2str(n), extensionM);
    end
    
    % Labelling the training set of images

    imageCount = 1;
    load(matname);
    numberOfImages = size(imgs,2);
    idx = randperm(numberOfImages);
    
    for k = 1:numberOfImages
        img = imgs{idx(k)};
        imwrite(img,sprintf('trainingImages/image_%.2d.jpg',imageCount),'jpg','Quality',100);
        imageCount = imageCount + 1;
    end

    % Building the AAM from the labelled data

    load('speech001_aam.mat');
    aamobj = build(aamobj,'trainingImages/','jpg',0.95,0.8);
    save('speech001_aam1.mat','aamobj');


    % Fitting the model to image sequence

    load('speech001_aam1.mat');
    load('trainingImages/landmarks/image_01_lm.mat');
    bs = get_s(aamobj,pts);
    save('initial_guess.mat','bs');

    load('speech001_aam1.mat');
    load('speechMAT/speech001.mat');  %'imageSequence.mat' 
    load('initial_guess.mat');
    fit = fitp(aamobj,imgs,bs,1,1);
    save('imageSequence_fit.mat','fit');
end

end