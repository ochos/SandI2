%-------------- visual lab1 -------------------%
% load('rect.mat')
% vidreader = VideoReader('Data/Training Data/speech001.mov');
% 
% % vidFrames = read(vidreader);
% % too large, not recommended 
% 
% nframes = get(vidreader, 'NumberOfFrames');
% % for frame = 1:nframes
% %     img = read(vidreader,frame);
% % end
% 
% img = read(vidreader,10);
% %imshow(img);
% 
% %roirect = getrect; % specify by dragging
% cimg = imcrop(img,r);
% imshow(cimg);

load('Data mat/Training Data/videoList.mat');
fileList = dir('Data/Training Data/*.mov');
load('Data mat/Training Data/videoList.mat');
startW = 'speechMAT/training/speech';
startM = 'ASR/MFCCs/Training movie/speech';
startMov = 'Data mat/Training Data/v';
number1to9 = '00';
number10to20 = '0';
extensionMfcc = '.mfc';
extensionMat = '.mat';


for f = 1:length(fileList)
    
    if f >= 1 && f <=9
        matname = strcat(startW, number1to9, int2str(f), extensionMat);
        mfcName = strcat(startM, number1to9, int2str(f));
        vidreader = strcat(startMov, number1to9, int2str(f), extensionMat);
    else
        matname = strcat(startW, number10to20, int2str(f), extensionMat);
        mfcName = strcat(startM, number10to20, int2str(f));
        vidreader = strcat(startMov, number10to20, int2str(f), extensionMat);
    end
    
    % Labelling the training set of images
    imageCount = 1;
    load(matname);
    numberOfImages = size(imgs,2);
    cimg = imgs{10};
    coeff = zeros(size(cimg,1)*size(cimg,2),numberOfImages);
    
    for frames = 1:numberOfImages
        coeff_frame = dct2(rgb2gray(imgs{frames}));
        imshow(rgb2gray(imgs{frames}));
        %imshow(log(abs(J)),[]), colormap(jet(64)), colorbar
        coeff(:,frames) = coeff_frame(:);
    end


    samPeriod = 1/get(vidreader, 'FrameRate');  % FrameRate is not constant by QuickTime
    %htk_write_mfc(mfcName,size(coeff,2),samPeriod,4*size(coeff,1),9,coeff);
    writehtk(mfcName,coeff,samPeriod,9)
end


coeff = zeros(size(img,1)*

% coeff = zeros(size(cimg,1)*size(cimg,2),nframes);
% for frames = 1:nframes
%     img = read(vidreader,frames);
%     cimg = rgb2gray(imcrop(img,r));
%     coeff_frame = dct2(cimg);
%     %imshow(log(abs(J)),[]), colormap(jet(64)), colorbar
%     coeff(:,frames) = coeff_frame(:);
% end

% coeff(abs(coeff) < 10) = 0;
% re_coeff = idct2(coeff);
% % imshow(I)
% figure,imshow((re_coeff,10),[0 255]);

% samPeriod = 1/get(vidreader, 'FrameRate');  % FrameRate is not constant by QuickTime
% htk_write_mfc('speech001.mfc',size(coeff,2),samPeriod,4*size(coeff,1),9,coeff);