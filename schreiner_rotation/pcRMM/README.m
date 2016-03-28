% A script that demonstrates how the provided functions can be used.

% Copyright (c) 2015 Michael Okun, michael.okun@mail.huji.ac.il
% License: GPLv2

% Load a raster representing ~5 min. of activity 
% (using 2 ms bins) on 10 electrodes
%load rn8_raster_small_logical.mat;
raster=raster_smaller;
N = size(raster,1);

% The word distribution in the actual physiological data:
dist = raster2dist(raster);

% Compute raster and distribution under the assumption of spike trains being 
% independent. The mean firing rate (MFR) of each spike train is preserved.
rasterMFR = false(size(raster));
s = sum(raster,2)';
for i = 1:N
  p = randperm(size(raster,2));
  rasterMFR(i,p(1:s(i))) = true;
end;
distMFR = raster2dist(rasterMFR);

% Compute the distribution according to the raster marginals model:
prd = histc(sum(raster), 0:N);
rasterRMM = RMM(s, prd);
distRMM = raster2dist(rasterRMM);

% Compute the distribution according to the raster marginals model with coupling:
c = sum(raster) * raster'; % population coupling of each unit: its inner product with population rate
rasterCRMM = cRMM(s, prd, c);
distCRMM = raster2dist(rasterCRMM);

FontSize = 12;
% Plot the actual word distribution vs. prediction by MFRs only:
WorddistSimilarityPlot(dist, distMFR);
title('MFR only', 'FontWeight', 'bold', 'FontSize', FontSize);
xlabel('observed word probability', 'FontWeight', 'bold', 'FontSize', FontSize);
ylabel('predicted word probability', 'FontWeight', 'bold', 'FontSize', FontSize);

% Plot the actual word distribution vs. prediction by raster marginals model:
WorddistSimilarityPlot(dist, distRMM);
title('raster marginals', 'FontWeight', 'bold', 'FontSize', FontSize);
xlabel('observed word probability', 'FontWeight', 'bold', 'FontSize', FontSize);
ylabel('predicted word probability', 'FontWeight', 'bold', 'FontSize', FontSize);

% Plot the actual word distribution vs. prediction by raster marginals with coupling model:
WorddistSimilarityPlot(dist, distCRMM);
title('raster marginals+coupling terms', 'FontWeight', 'bold', 'FontSize', FontSize);
xlabel('observed word probability', 'FontWeight', 'bold', 'FontSize', FontSize);
ylabel('predicted word probability', 'FontWeight', 'bold', 'FontSize', FontSize);
% Since the data in raster comes from multi-unit, rather than single-unit
% recording, for this particular example we do not expect a major difference between 
% raster marginals and raster marginals with population coupling terms.
