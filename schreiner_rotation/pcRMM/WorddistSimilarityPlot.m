function WorddistSimilarityPlot(dist1, dist2)
% WorddistSimilarityPlot(dist1, dist2) 
%
% Plots the probability of each pattern (word) in the first distribution vs. 
% its probability in the other distribution.
%
% Input: dist1, dist2 - the two distributions on binary words of length N,
% i.e., dist1 and dist2 are vectors of length 2^N,  where the value 
% in position k (1<=k<=2^N) indicates the number of times  
% the word dec2bin(k-1, N) was observed. 
%
% Output: a new figure.

% Copyright (c) 2012 Michael Okun, michael.okun@mail.huji.ac.il
% License: GPLv3


MarkerSize = 5;
figure; hold on;
N = log2(length(dist1));

if length(dist1) ~= length(dist2) || N ~= floor(N)
  error('WorddistSimilarityPlot: bad input arguments');
end;

S = sum(dec2bin(0:length(dist1)-1)=='1', 2);
% S counts the numbers of 1s in each binary word of length N

dist1(dist1 == 0) = 1;
dist2(dist2 == 0) = 1; % 0s cannot be seen on log-scale

plot(dist1(1)/sum(dist1), dist2(1)/sum(dist2), 'o', ...
  'MarkerFaceColor', [0 0 0], ...
  'MarkerEdgeColor', 'none', ...
  'MarkerSize', MarkerSize); % 0 word is in black
for spk = 1:N
  plot(dist1(S==spk)/sum(dist1), dist2(S==spk)/sum(dist2), 'o' , ...
    'MarkerFaceColor', [spk/N,  1-abs(spk-N/2)/(N/2), 1-spk/N], ...
    'MarkerEdgeColor', 'none', ...
    'MarkerSize', MarkerSize);  
end;

xlim([0.5/max(sum(dist1), sum(dist2)) 1]);
ylim([0.5/max(sum(dist1), sum(dist2)) 1]);
plot(xlim, ylim, '--', 'Color', 0.5*[1 1 1]);
set(gca, 'XScale', 'log');
set(gca, 'YScale', 'log');

