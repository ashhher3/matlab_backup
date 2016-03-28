function [ ecog ] = combineFiles( foldername, idx )

h5files=dir(sprintf('%s/*.h5', foldername));

ecog=[];
r=0;

for f=1:length(h5files)
    filename=sprintf('%s/%s',foldername, char(h5files(f).name));
    
    new = h5read(filename, '/ECoG Array');
    
    ecog = [ecog, new(idx(1):idx(2),mod(1+r,4):4:end)];
    
    l=size(new,2);
    clear new
    
    r=4-mod(l-r,4);
    
    fprintf('done with %s\n', filename);
    
end

end

%%


ofc3_filt=NaN(size(ofc3));
fprintf('running eegfilt...\n')
for i=1:10
    ofc3_filt(i,:)=eegfilt(ofc3(i,:),256,1,[]);
end
fprintf('done with eegfilt, subtracting mean...\n')
for i=1:10
    offset=mean(ofc3_filt(i,:));
    ofc3_filt(i,:)=ofc3_filt(i,:)-offset;
end
fprintf('saving ofc3_filt...\n')
save('ofc3_filt','ofc3_filt');
clear ofc3 offset

fprintf('done with ofc3_filt, filtering for beta...\n')
ofc3_beta=NaN(size(ofc3_filt));
butterBeta = designfilt('bandpassiir','FilterOrder',4, ...
    'HalfPowerFrequency1',13,'HalfPowerFrequency2',30, ...
    'SampleRate',256);
for i=1:10
    ofc3_beta(i,:)=filtfilt(butterBeta,ofc3_filt(i,:));
end
fprintf('saving ofc3_beta...\n')
save('ofc3_beta','ofc3_beta');
clear ofc3_filt butterBeta

fprintf('done with ofc3_beta, calculating power...\n')
ofc3_betaPwr=NaN(size(ofc3_beta));
for i=1:10
    ofc3_betaPwr(i,:)=abs(hilbert(ofc3_beta(i,:)));
end
fprintf('saving ofc3_betaPwr...\n')
save('ofc3_betaPwr','ofc3_betaPwr')
fprintf('done!\n')
