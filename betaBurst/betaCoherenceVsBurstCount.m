function [bCoh,bCount] = betaCoherenceVsBurstCount(area1, area2)
%% betaCoherenceVsBurstCount.m
% 
% Compares cross-stucture coherence and number of beta bursts for two brain
% areas. Function assumes data is in Starr lab style data structures.
%
%%

% make sure just pass one pair of channels
% if (length(area1.contact)>1) || (length(area2.contact)>1)
%     error('you must pick a single channel for each region')
% end

% get fs
if (area1.Fs(1)==area2.Fs(1)) && (area1.FsB(1)==area2.FsB(1))
    fs = area1.Fs(1);
    fsB = area1.FsB(1);
else
    error('all channels must have same sampling frequency')
end

% get end time
if ((length(area1.contact(1).signal)/fs) == (length(area1.contact(1).beta)/fsB))
    timeEnd = length(area1.contact(1).signal)/fs;
else
    error('recordings end at different times')
end

% get number of channels
numA1 = size(area1.contact,2);
numA2 = size(area2.contact,2);


% okay now actually step through the recording and calculate things
wSig=10*fs; % 10sec window for complete signal
wBeta=10*fsB; % 10sec window for beta bursts
bCoh = NaN(numA1,numA2,floor(timeEnd/10));
bCount = NaN(numA1,numA2,floor(timeEnd/10));
for sec = 0:10:timeEnd-10
    % indicies of this window for complete signal
    sSig=floor(sec*fs) +1;
    eSig=sSig + wSig - 1;
    
    %indices of this window for beta bursts
    sBeta=floor(sec*fsB) +1;
    eBeta=sBeta + wBeta - 1;
    
    % make data the right shape for cohmatrixc.m
    data = NaN(wSig, numA1+numA2);
    for i=1:numA1
        data(:,i) = [area1.contact(i).signal(sSig:eSig)];
    end
    for i=1:numA2
        data(:,numA1+i) = [area2.contact(i).signal(sSig:eSig)];
    end
    
    % call cohmatrixc.m
    params.Fs = fs;
    params.fpass = [13 30];
    C = cohmatrixc(data,params);
    Cavg = mean(C); % should average right dimension but stay 3D
    
    % record average beta coherence & num beta bursts for each channel pair
    for i=1:numA1
        for j=1:numA2
            bCoh(i,j, sec/10+1) = Cavg(1, i,j+numA1);
 
            bCount(i,j,sec/10+1) = sum(area1.contact(i).burst(sBeta:eBeta)) ...
                + sum(area2.contact(j).burst(sBeta:eBeta));
        end
    end
    
end


end