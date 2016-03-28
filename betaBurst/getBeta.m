function data = getBeta( data )
% Takes Starr lab style data structure and computes beta trace for each
% channel:
%       1. convert from uV to V
%       2. lowpass filter @ 100 Hz
%       3. downsample to 220Hz
%       4. bandpass filter @ 13-30 Hz

%% make filters for later
Fs = data.Fs(1); % ASSUMING ALL CHANNELS HAVE SAME SAMPLING FREQ
butter100low = designfilt('lowpassiir','FilterOrder',8,...
    'HalfPowerFrequency',100,'SampleRate',round(Fs),'DesignMethod','butter');
butterBeta = designfilt('bandpassiir','FilterOrder',4, ...
    'HalfPowerFrequency1',13,'HalfPowerFrequency2',30, ...
    'SampleRate',220);

%% CONVERT TO VOLTS
scale = 1e-6; %1e-2;

for i=1:length(data.contact)
    for j=1:length(data.contact(i).signal)
        data.contact(i).beta(j) = scale*data.contact(i).signal(j);
    end
end

%% lowpass filter @ 100 Hz

for i=1:length(data.contact)
    data.contact(i).beta = filtfilt(butter100low,data.contact(i).beta);
end

%% downsample to 220Hz

newF = 220; 

for i=1:length(data.contact)
    data.contact(i).beta = resample(data.contact(i).beta,newF,round(data.Fs(i)));
    data.FsB(i) = newF;
end


%% bandpass filter @ 13-30 Hz

for i=1:length(data.contact)
    data.contact(i).beta = filtfilt(butterBeta,data.contact(i).beta);
end

end

