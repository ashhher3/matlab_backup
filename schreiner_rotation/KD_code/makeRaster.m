function raster= makeRaster(spk)
% takes a spk structure from a spk-strf.mat file and turns it into a 
% raster matrix for use with Okun code
%
% (KD 2015/05/05)

n = size(spk,2); % number of units

% find time of first and last spikes
first=NaN;
last = -1.0;
for i=1:n
   first=min(first,spk(i).spiketimes(1));
   last=max(last, spk(i).spiketimes(end));
end

% initialize raster using 1ms bins
m=fix((last-first))+1;
raster=false(n,m);
for i=1:n % loop over all units (rows)
    for k=1:size(spk(i).spiketimes,2) % loop over all spiketimes
        timebin = fix((spk(i).spiketimes(k)-first))+1;
        % raster(i,timebin) = raster(i,timebin) + 1;
        raster(i,timebin)= true;
    end
end


end