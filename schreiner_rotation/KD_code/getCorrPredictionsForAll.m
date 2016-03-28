matlabmail('kderosier6@gmail.com','starting rn4','getCorrPredictionsForAll');

% load rn4
load('141215_234435-site6-800um-20db-rn4-fs20000-A-spk-strf-raster');

% run predict corrs for rn4
%newPred=realPrediction(rasterSort);

% save rn4
%save('KD_RasterData/site6/rn4/new_predictions4.mat');
clear all;

% email me to say done with rn4 and moving on to rn8
matlabmail('kderosier6@gmail.com','done with rn4, starting rn8','getCorrPredictionsForAll');

% load rn 8
load('141215_234435-site6-800um-20db-rn8-fs20000-A-spk-strf-raster');

% run predictPairCorrs for rn8
newPred=realPrediction(rasterSort);

% save rn8
save('KD_RasterData/site6/rn8/new_predictions8.mat');
clear all;

% email me to say done with rn8 and moving on to rn16
matlabmail('kderosier6@gmail.com','done with rn8, starting rn16','getCorrPredictionsForAll');

% load rn16
load('141215_234435-site6-800um-20db-rn16-fs20000-A-spk-strf-raster');

% run predictPairCorrs for rn16
newPred=realPrediction(rasterSort);

% save rn16
save('KD_RasterData/site6/rn16/new_predictions16.mat');
clear all;

% email me to say done with everything
matlabmail('kderosier6@gmail.com','done with rn16, come back!','getCorrPredictionsForAll');