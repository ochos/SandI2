% dct%% 
img = imread('trainingImages/image_203.jpg');
img = rgb2gray(img);
figure(1)
imshow(img);

coeff = dct2(img);

img2 = uint8(idct2(coeff));
err = imsubtract(img,img2);
figure(1);
imshow(img2);
figure(2);
imagesc(coeff);

waitforbuttonpress;

coeff(1,1)= 0;
figure(2);
imagesc(coeff);


coeff = dct2(img);
coeff(:,16:end) = 0; coeff(16:end,:) = 0;
img2 = uint8(idct2(coeff));
figure(1);
imshow(img2);
figure(2);
imagesc(coeff);

%load speech sequence
load('speechMAT/training/speech001.mat');

n = 10;
mov = [];
for k = 1:length(imgs)
    img = rgb2gray(imgs{k});
    coeff = dct2(img);
    coeff(:,n+1:end) = 0; coeff(n+1:end,:) = 0;
    img2 = uint8(idct2(coeff));
    mov(k).cdata = img2;
    mov(k).colormap = gray(256);
end

fps = 25;
movie(mov,1,fps);


% truncate for HTK file. 10 by 10 numbers
%coeff = reshape(coeff(1:n,1:n),[],1); %100x1 double 
coeff = coeff(1:10,1:10);
vidreader = VideoReader('Data/Training Data/speech001.mov');
samPeriod = 1/get(vidreader, 'FrameRate'); 
writehtk('ASR/MFCCs/Training/speech001_DCT.mfc',coeff,samPeriod,9);