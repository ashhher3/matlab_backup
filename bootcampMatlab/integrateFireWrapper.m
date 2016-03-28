function integrateFireWrapper(data,thresholds)
% Wrapper function for integrate and fire demonstration; takes a matrix of 
% data and an array of thresholds for firing decision, calls 
% integrateAndFire for each threshold, and plots rasters of each spike train.
%
% 8/3/2015 KD & NH

% Initialize space for spike train results
spTrains = cell(1,length(thresholds));

% Loop over threshold array and generate the spike train for each
% entry
for i=1:length(thresholds)
    thresh = thresholds(i); % define single threshold for this iteration
    spTrains{i} = integrateAndFireSolution(data, thresh); % THIS IS THE FUNCTION YOU HAVE TO WRITE
end

% Generate plots of the spike trains
raster(spTrains,thresholds);

end

