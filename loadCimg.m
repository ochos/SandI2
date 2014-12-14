function [cimg] = loadCimg(filename,rect)

    vidreader = VideoReader(filename);
    nframes = get(vidreader, 'NumberOfFrames');

    img = read(vidreader,10);
    cimg = imcrop(img,roirect);
    
    coeff = zeros(size(cimg,1)*size(cimg,2),nframes);
    for frames = 1:nframes
        img = read(vidreader,frame);
        cimg = rgb2gray(imcrop(img,rect));
        coeff_frame = dct2(cimg);
        coeff(:,frame) = coeff_frame(:);
    end

end