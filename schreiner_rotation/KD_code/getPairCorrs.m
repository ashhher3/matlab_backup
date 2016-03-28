function corrs = getPairCorrs(raster)
% function to compute pairwise spike correlations using Pearson correlation
% as given by formula in Cutts & Eglen 2014 J Neurosci
%
% rows in raster should be units, columns should be timebins
%
% (KD 2015/05/27)

n = size(raster,1);

corrs=zeros(n); %initialize empty correlation matrix

for i=1:n % loop over all neurons
    corrs(i,i)=1.0; % we know diagonal is 1
    iAvg=mean(raster(i,:)); % avg FR of neuron i
    for j=i+1:n % loop over all new pairs
        jAvg=mean(raster(j,:)); % avg FR of neuron j
        corrs(i,j)=(sum((raster(i,:)-iAvg).*(raster(j,:)-jAvg)))/(sqrt(sum((raster(i,:)-iAvg).*(raster(i,:)-iAvg)))*sqrt(sum((raster(j,:)-jAvg).*(raster(j,:)-jAvg))));
        corrs(j,i)=corrs(i,j); % pearson corr is symmetric
    end
end

end

