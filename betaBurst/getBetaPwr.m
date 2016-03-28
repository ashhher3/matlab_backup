function [ data ] = getBetaPwr(data)

n=size(data.signal,1);

butterBeta = designfilt('bandpassiir','FilterOrder',4, ...
    'HalfPowerFrequency1',13,'HalfPowerFrequency2',30, ...
    'SampleRate',data.fs);

data.beta = [];
data.betaPwr = [];

for j=1:n
    data.beta(j,:) = filtfilt(butterBeta, data.signal(j,:));
    data.betaPwr(j,:) = abs(hilbert(data.beta(j,:)));
end

end