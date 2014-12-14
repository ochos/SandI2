fileList = dir('Data/Test Data/*.mov');
%vidobj = VideoReader(strcat('Data/Training Data/',fileList(1).name));
%mov = read(vidobj);
%imshow(mov(:,:,:,300));
%r = getrect;
load('rect.mat');
for f = 1:length(fileList)
   vidobj = VideoReader(strcat('Data/Test Data/',fileList(1).name));
   mov = read(vidobj);
   imgs = {};
   for k = 1:size(mov,4);
      imgs{k} = imcrop(mov(:,:,:,k), r);
   end
   save(sprintf('speechMAT/testing%s.mat', fileList(f).name(1:end-4)), 'imgs');
end
