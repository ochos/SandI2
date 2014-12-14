% load 'defined_aam.mat';
% aamobj = build(aamobj, 'trainingImagesJPG', 'png', 0.95,0.8);
% save('definedBuilt_amm.mat','aamobj');
% 
% load('definedBuilt_amm.mat');
% load('trainingImagesJPG/landmarks/image_0001_lm.mat');
% bs = get_s(aamobj, pts);
% save('initial_guess.mat','bs');
% 
% load('definedBuilt_amm.mat');
% load('speechMat/training/speech001.mat');
% load('initial_guess.mat');
% fit = fitp(aamobj,imgs,bs,1,1);
% save('speechMat/trainingFit/speech001_fit.mat','fit');
% 
% load('CORNERSwith2Loops_aam.mat');
% load('trainingImagesJPG/landmarks/image_0001_lm.mat');
% bs = get_s(aamobj, pts);
% save('initial_guess.mat','bs');
% 
% load('CORNERSwith2Loops_aam.mat');
% load('speechMat/training/speech002.mat');
% load('initial_guess.mat');
% fit = fitp(aamobj,imgs,bs,1,1);


for k = 1:20
    load('CORNERSwith2Loops_aam.mat');
    load(sprintf('speechMat/training/speech%.3d.mat',k));
    load('initial_guess.mat');
    fit = fitp(aamobj,imgs,bs,0,1);
    save(sprintf('speechMat/trainingFitCORNERSwith2Loops/speech%.3d_fit.mat',k),'fit');
end