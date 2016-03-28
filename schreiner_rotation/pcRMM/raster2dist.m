function dist = raster2dist(raster)
% dist = raster2dist(raster)
%
% Input: raster is a 0/1 matrix of dimensions 
% N by (no. of time bins), where N is the no. of spike trains.
%
% Output: dist is a vector of length 2^N,
% where the value in position k indicates the number of times 
% the word dec2bin(k-1, N) appears in the raster

% Copyright (c) 2012 Michael Okun, michael.okun@mail.huji.ac.il
% License: GPLv3

raster = logical(raster);
N = size(raster,1); 

dist = histc( (2.^(0:N-1))*raster, 0:2^N-1 );