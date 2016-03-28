function corrs = predictPairCorrs(raster)
% function to predict pairwise correlations using Okun cRMM.m code
%
% can't run cRMM on a full population, so this is based on choosing random
% subpopulations and averaging.
%
% (KD 2015/05/27)

n=size(raster,1); % real size of population
k=25; % size of subpopulations, KD laptop can do 20, banshee can do 25
m = 10*n; % how many times to choose a subpop

corrs=zeros(n); % corr values will go here
denoms=zeros(n); %counts how many times each neuron pair appears in subpops

for h=1:m % loop to choose random subpopulation
    perm=randsample(n,k); %indices of the subpop
    copy=raster(perm,:); %copy(i,:)=raster(perm(i),:);
    
    % use Okun code from README.m to get inputs for cRMM.m
    s = sum(copy,2)';
    prd = histc(sum(copy), 0:k);
    c = sum(copy) * copy';
    CRMM = cRMM(s, prd, c);
    
    % compute pairwise correlations and add to corrs and denoms
    for i=1:k % loop over all neurons
    corrs(perm(i),perm(i))=1.0; % we know diagonal is 1
    iAvg=mean(CRMM(i,:)); % avg FR of neuron i
        for j=i+1:k % loop over all new pairs
            jAvg=mean(CRMM(j,:)); % avg FR of neuron j
            Pcorr=(sum((CRMM(i,:)-iAvg).*(CRMM(j,:)-jAvg)))/(sqrt(sum((CRMM(i,:)-iAvg).*(CRMM(i,:)-iAvg)))*sqrt(sum((CRMM(j,:)-jAvg).*(CRMM(j,:)-jAvg))));
            corrs(perm(i),perm(j))=corrs(perm(i),perm(j))+Pcorr; % add to corrs
            denoms(perm(i),perm(j))=denoms(perm(i),perm(j))+1; % increment denoms
        end
    end
end

% more looping to get final correlation values
for i=1:n
    for j=i+1:n
        corrs(i,j)=corrs(i,j)/denoms(i,j);
        corrs(j,i)=corrs(i,j); % Pearson corr is symmetric
    end
end

end

