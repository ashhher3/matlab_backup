%% 0 - make filters for later
Fs = lfp.Fs(1); % ASSUMING ALL CHANNELS HAVE SAME SAMPLING FREQ
butter100low = designfilt('lowpassiir','FilterOrder',8,'HalfPowerFrequency',100,'SampleRate',round(Fs),'DesignMethod','butter');
butterBeta = designfilt('bandpassiir','FilterOrder',4, ...
    'HalfPowerFrequency1',13,'HalfPowerFrequency2',30, ...
    'SampleRate',220);

clear Fs
%% 1 - lowpass filter @ 100 Hz

for i=1:length(ecog.contact)
    ecog.contact(i).beta = filtfilt(butter100low,ecog.contact(i).signal);
end

for i=1:length(lfp.contact)
    lfp.contact(i).beta = filtfilt(butter100low,lfp.contact(i).signal);
end

clear i

%% 2 - downsample to 220Hz

newF = 220; 

for i=1:length(ecog.contact)
    ecog.contact(i).beta = resample(ecog.contact(i).beta,newF,round(ecog.Fs(i)));
    ecog.FsB(i) = newF;
end

for i=1:length(lfp.contact)
    lfp.contact(i).beta = resample(lfp.contact(i).beta,newF,round(lfp.Fs(i)));
    lfp.FsB(i) = newF;
end

clear i newF 

%% 3 - bandpass filter @ 13-30 Hz

for i=1:length(ecog.contact)
    ecog.contact(i).beta = filtfilt(butterBeta,ecog.contact(i).beta);
end

for i=1:length(lfp.contact)
    lfp.contact(i).beta = filtfilt(butterBeta,lfp.contact(i).beta);
end

clear i

%% CONVERT TO VOLTS

scale = 1e-6; %1e-2;

for i=1:length(ecog.contact)
    for j=1:length(ecog.contact(i).beta)
        ecog.contact(i).beta(j) = scale*ecog.contact(i).beta(j);
    end
end

for i=1:length(lfp.contact)
    for j=1:length(lfp.contact(i).beta)
        lfp.contact(i).beta(j) = scale*lfp.contact(i).beta(j);
    end
end

clear i j scale

%% 4 - square the value of each sample
% 
% 
% for i=1:length(ecog.contact)
%     for j=1:length(ecog.contact(i).beta)
%         ecog.contact(i).betaPwr(j) = ecog.contact(i).beta(j)^2;
%     end
% end
% 
% for i=1:length(lfp.contact)
%     for j=1:length(lfp.contact(i).beta)
%         lfp.contact(i).betaPwr(j) = lfp.contact(i).beta(j)^2;
%     end
% end
% 
% clear i j

%% 4.5 - OR, USE HILBERT POWER

for i=1:length(ecog.contact)
    ecog.contact(i).betaPwr = abs(hilbert(ecog.contact(i).beta));
end

for i=1:length(lfp.contact)
    lfp.contact(i).betaPwr = abs(hilbert(lfp.contact(i).beta));
end

clear i j

%% 5 - smooth with Hanning window 35 samples wide

% w = hann(35, 'symmetric');
% 
% for i=1:length(ecog.contact)
%     ecog.contact(i).betaPwr = conv(ecog.contact(i).betaPwr,w,'same');
% end
% 
% for i=1:length(lfp.contact)
%     lfp.contact(i).betaPwr = conv(lfp.contact(i).betaPwr,w,'same');
% end
% 
% clear w i

%% 6 - calculate median for whole trace
for i=1:length(ecog.contact)
    ecog.contact(i).medB = median(ecog.contact(i).betaPwr);
    ecog.contact(i).thresh = 3*ecog.contact(i).medB;
end

for i=1:length(lfp.contact)
    lfp.contact(i).medB = median(lfp.contact(i).betaPwr);
    lfp.contact(i).thresh = 3*lfp.contact(i).medB;
end

clear i

%% 6 B - OR, IF STIM FILE, only calculate median for pre-stim samples

% stimStart = 4140; %13750; %14520; % estimate sample number in betaPwr file
% 
% for i=1:length(ecog.contact)
%     ecog.contact(i).medB = median(ecog.contact(i).betaPwr(1:stimStart));
%     ecog.contact(i).thresh = 3*ecog.contact(i).medB;
% end
% 
% for i=1:length(lfp.contact)
%     lfp.contact(i).medB = median(lfp.contact(i).betaPwr(1:stimStart));
%     lfp.contact(i).thresh = 3*lfp.contact(i).medB;
% end
% 
% clear i stimStart

%% 6 C - ALTERNATE THRESHOLD, 3x std dev
% for i=1:length(ecog.contact)
%     ecog.contact(i).stdDevB = std(ecog.contact(i).betaPwr);
%     ecog.contact(i).thresh = ecog.contact(i).stdDevB;
% end
% 
% for i=1:length(lfp.contact)
%     lfp.contact(i).stdDevB = std(lfp.contact(i).betaPwr);
%     lfp.contact(i).thresh = 3*lfp.contact(i).stdDevB;
% end
% 
% clear i
%% 7 - mark a burst every time go above threshold

for i=1:length(ecog.contact)
    for j=1:length(ecog.contact(i).betaPwr)
        if ecog.contact(i).betaPwr(j)>ecog.contact(i).thresh
            ecog.contact(i).burst(j)=1;
        else
            ecog.contact(i).burst(j)=0;
        end
    end
end

for i=1:length(lfp.contact)
    for j=1:length(lfp.contact(i).betaPwr)
        if lfp.contact(i).betaPwr(j)>lfp.contact(i).thresh
            lfp.contact(i).burst(j)=1;
        else
            lfp.contact(i).burst(j)=0;
        end
    end
end

clear i j

%% 8 - length of each burst
Fs=ecog.FsB(2);
betaTime = 0:1/Fs:length(ecog.contact(1).betaPwr)/Fs-1/Fs;

for i=1:length(ecog.contact)
samp=1;
burst=false;
ecog.contact(i).startTimes=[];
ecog.contact(i).endTimes=[];
while samp<length(ecog.contact(i).betaPwr)
    while burst==false && (samp<length(ecog.contact(i).betaPwr))
        if ecog.contact(i).burst(samp)==1
            ecog.contact(i).startTimes=[ecog.contact(i).startTimes, betaTime(samp)];
            burst=true;
        end
        samp=1+samp;
    end
    while burst==true (samp<length(ecog.contact(i).betaPwr))
        if ecog.contact(i).burst(samp)~=1
            ecog.contact(i).endTimes=[ecog.contact(i).endTimes, betaTime(samp-1)];
            burst=false;
        end
        samp=1+samp;
    end
end
end

for i=1:length(lfp.contact)
samp=1;
burst=false;
lfp.contact(i).startTimes=[];
lfp.contact(i).endTimes=[];
while samp<length(lfp.contact(i).betaPwr)
    while burst==false && (samp<length(lfp.contact(i).betaPwr))
        if lfp.contact(i).burst(samp)==1
            lfp.contact(i).startTimes=[lfp.contact(i).startTimes, betaTime(samp)];
            burst=true;
        end
        samp=1+samp;
    end
    while burst==true && (samp<length(lfp.contact(i).betaPwr))
        if lfp.contact(i).burst(samp)~=1
            lfp.contact(i).endTimes=[lfp.contact(i).endTimes, betaTime(samp-1)];
            burst=false;
        end
        samp=1+samp;
    end
end
end


clear i j samp burst