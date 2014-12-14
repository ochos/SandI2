fileList = dir('Data/Training Data/*.mov');
%vidobj = VideoReader(strcat('Data/Training Data/',fileList(1).name));
%mov = read(vidobj);
%imshow(mov(:,:,:,300));
%r = getrect;
% load('rect.mat');
% for f = 1:length(fileList)
%    vidobj = VideoReader(strcat('Data/Training Data/',fileList(f).name));
%    mov = read(vidobj);
%    imgs = {};
%    for k = 1:size(mov,4);
%       imgs{k} = imcrop(mov(:,:,:,k), r);
%    end
%    save(sprintf('speechMAT/training/%s.mat', fileList(f).name(1:end-4)), 'imgs');
% end

% % for test data
% fileList = dir('Data/test Data/*.mov');
% load('rect.mat');
% for f = 1:length(fileList)
%    vidobj = VideoReader(strcat('Data/test Data/',fileList(f).name));
%    mov = read(vidobj);
%    imgs = {};
%    for k = 1:size(mov,4);
%       imgs{k} = imcrop(mov(:,:,:,k), r);
%    end
%    save(sprintf('speechMAT/test/%s.mat', fileList(f).name(1:end-4)), 'imgs');
% end
%------------------------------- DONE -------------------------------------%

% Defining the shape of the AAM
aamdefine

% Labelling the training set of images

imageCount = 1;
load('speechMAT/training/speech001.mat');
numberOfImages = size(imgs,2);
idx = randperm(numberOfImages);
for k = 1:numberOfImages
    img = imgs{idx(k)};
    imwrite(img,sprintf('trainingImages/image_%.2d.jpg',imageCount),'jpg','Quality',100);
    imageCount = imageCount + 1;
end



aamplace


% Building the AAM from the labelled data

load('defined_aam.mat');
aamobj = build(aamobj,'trainingImages/','jpg',0.95,0.8);
save('defined_aam1.mat','aamobj');


% Fitting the model to image sequence

load('defined_aam1.mat');
load('trainingImages/landmarks/image_01_lm.mat');
bs = get_s(aamobj,pts);
save('initial_guess.mat','bs');



load('defined_aam1.mat');
load('speechMAT/training/speech002.mat');  %'imageSequence.mat' 
load('initial_guess.mat');
fit = fitp(aamobj,imgs,bs,1,1);
save('imageSequence_fit.mat','fit');