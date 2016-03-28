%% betaburstcoherence.m
% - model LFP as gaussian noise plus trains of beta bursts
% - there's 2 areas, area1 & area2, and 2 patterns of beta bursts, A & B
% - pick the number of bursts in pattern A & pattern B, and the duration of
%       the burst trains. bursts are evenly spaced within a train, and 
%       area2 has a phase shift from area1. it would be better to place the
%       bursts randomly.
% - the fake data is divided into 10 second blocks. some fraction of the
%       blocks are randomly chosen to be burst blocks. burst trains are
%       randomly placed within the burst blocks. random jitter is chosen 
%       and data is rotated to shift the blocks a little.
% - coherence is calculated for 10 second non-overlapping blocks.
% - make plots: fake data, coherence heat map, beta coherence histogram
%
% Kate Derosier 12/17/2015
%
%% sampling rate and time points
fs=500;
m=50; % how many minutes of fake data do we want
timeSteps=0:(1/fs):m*60; % time steps in seconds

clear m

%% make trains of beta bursts to add to fake data

% make 2-cycle beta burst for each area
betaF = 20;
pshift = pi/12;
amp = 1;

burst1 = amp*sin(2*pi*betaF*timeSteps(1:200));
burst2 = amp*sin(2*pi*betaF*timeSteps(1:200) + pshift);

% make A and B beta burst trains for each area
trainDurationA = 2.75*fs; % length of burst train in samples
trainDurationB = 2.75*fs;
burstLength = length(burst1); % length of burst in samples
numA = 16; % number of bursts in a train
numB = 7;

% make A trains
beta1A = zeros(1,trainDurationA);
beta2A = zeros(1,trainDurationA);
intervalA = floor((trainDurationA - numA*burstLength)/(numA-1));
for i=1:burstLength+intervalA:trainDurationA-burstLength+1
    beta1A(i:i+burstLength-1)=burst1;
    beta2A(i:i+burstLength-1)=burst2;
end

% make B trains
beta1B = zeros(1,trainDurationB);
beta2B = zeros(1,trainDurationB);
intervalB = floor((trainDurationB - numB*burstLength)/(numB-1));
for i=1:burstLength+intervalB:trainDurationB-burstLength+1
    beta1B(i:i+burstLength-1)=burst1;
    beta2B(i:i+burstLength-1)=burst2;
end

clear i pshift amp burstLength numA numB intervalA intervalB burst1 burst2

%% make fake data that is uncorrelated except for these bursts

rng default

% make some gaussian noise of the right length
area1A = randn(size(timeSteps)); 
area1B = randn(size(timeSteps)); 
area2A = randn(size(timeSteps));
area2B = randn(size(timeSteps));

% pick out blocks to add burst trains to
burstBlocksA = [];
burstBlocksB = [];
blockLength = 10*fs;
pr=5/8;
for b=1:(size(timeSteps,2)-1)/blockLength
    if rand()>pr
        burstBlocksA = [burstBlocksA, b];
    end
    if rand()>pr
        burstBlocksB = [burstBlocksB, b];
    end
end

% place a burst train within each burst block
for b=1:length(burstBlocksA)
    s = (burstBlocksA(b)-1)*blockLength+1;
    e = burstBlocksA(b)*blockLength;
    
    i = randi([s,e-(trainDurationA-1)],1); % random burst start
    area1A(i:i+trainDurationA-1)=area1A(i:i+trainDurationA-1)+beta1A;
    area2A(i:i+trainDurationA-1)=area2A(i:i+trainDurationA-1)+beta2A;
end
 for b=1:length(burstBlocksB)
    s = (burstBlocksB(b)-1)*blockLength+1;
    e = burstBlocksB(b)*blockLength;
    
    i = randi([s,e-(trainDurationB-1)],1);
    area1B(i:i+trainDurationB-1)=area1B(i:i+trainDurationB-1)+beta1B;
    area2B(i:i+trainDurationB-1)=area2B(i:i+trainDurationB-1)+beta2B;
end

% rotate the data a little so that blocks aren't aligned
i = randi([1,blockLength/2],1);
area1A=circshift(area1A,i,2);
area2A=circshift(area2A,i,2);
area1B=circshift(area1B,i,2);
area2B=circshift(area2B,i,2);

clear burstBlocksA burstBlocksB pr b s e i beta1A beta1B beta2A beta2B betaF trainDurationA trainDurationB w

%% compute coherence for 10sec windows
windowLength=10*fs; % 10sec window
blockStep=5*fs; % slide the window by 5 sec each time
CA=[];
FA=[];
CB=[];
FB=[];
timeBlocks=[];
CA_beta = [];
CB_beta = [];
for i = 1:blockStep:numel(timeSteps)-windowLength+1
    s=i;
    e=i+windowLength-1;
    
    [CAnew,FAnew]=mscohere(area1A(s:e),area2A(s:e),[],[],[],fs);
    [CBnew,FBnew]=mscohere(area1B(s:e),area2B(s:e),[],[],[],fs);
    
    CA = [CA; CAnew'];
    FA = [FA; FAnew'];
    CB = [CB; CBnew'];
    FB = [FB; FBnew'];
    
    timeBlocks = [timeBlocks, timeSteps(s)];
    
    % compute avg over 13-30Hz, indices were figured out the stupid way
    CA_beta = [CA_beta, mean(CAnew(55:127))];
    CB_beta = [CB_beta, mean(CBnew(55:127))];
end

% a=1;
% b=[0.2 0.2 0.2 0.2 0.2];
% CA_beta = filter(b,a,CA_beta);
% CB_beta = filter(b,a,CB_beta);

clear i s e blockLength

%% plot the fake data
figure(10)
clf
set(gcf,'Color',[1 1 1])

subplot(2,1,1)
hold on
plot(timeSteps,area1A,'b')
plot(timeSteps,area2A,'c')
xlim([0 timeSteps(15*60*fs)])
title('Fake LFP data with burst pattern A')
ylim([-5 5])
hold off

subplot(2,1,2)
hold on
plot(timeSteps,area1B,'r')
plot(timeSteps,area2B,'m')
xlim([0 timeSteps(15*60*fs)])
title('Fake LFP data with burst pattern B')
ylim([-5 5])
hold off

%% make coherence plots
figure(11)
clf
set(gcf,'Color',[1 1 1])

subplot(1,2,1)
pcolor(timeBlocks,FA',CA')
shading flat
colorbar
ylim([0 250])
title('Coherence between area1 and area2 for beta burst pattern A')
hold off

subplot(1,2,2)
pcolor(timeBlocks,FB',CB')
shading flat
colorbar
ylim([0 250])
title('Coherence between area1 and area2 for beta burst pattern B')
hold off

%% make a histogram of beta coherence values
figure(12)
clf
set(gcf,'Color',[1 1 1])

hold on
histogram(CA_beta,fix(numel(timeBlocks)/20))
histogram(CB_beta,fix(numel(timeBlocks)/20))

xlim([0 0.5])
title('Histogram of beta band coherence values')
hold off


