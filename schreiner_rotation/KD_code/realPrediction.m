function corrs = realPrediction(raster)
% function to predict pairwise correlations using Okun cRMM.m code
%
% (KD 2015/05/27)

n=size(raster,1); % real size of population

% use Okun code from README.m to get inputs for cRMM.m
s = sum(raster,2)';
prd = histc(sum(raster), 0:n);
c = sum(raster) * raster';
CRMM = cRMM(s, prd, c);

% compute pairwise correlations of CRMM
corrs=getPairCorrs(CRMM);




