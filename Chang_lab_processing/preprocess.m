function preprocess(keep)

fprintf('Initializing...\n')

% SET RESAMPLED FS HERE
fs = 512;

% get the prefix for this set of files
[~,prefix,~] = fileparts(pwd);

% make dirs to store .mat files
mkdir('mat_allChannels')

% get list of .h5 files
h5files=dir('*.h5');

% get start and end times
times=h5read(char(h5files(1).name), '/timestamp vector');
globalStart=times(1)+60600; %%%% CHANGED THIS! PUT IT BACK AFTERWARDS!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
times=h5read(char(h5files(end).name), '/timestamp vector');
globalEnd=globalStart + 40*60; % times(end);  %%%% CHANGED THIS! PUT IT BACK AFTERWARDS!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
clear times

fprintf('Checking the sampling rate...\n')

% loop over h5 files
fs_orig = h5readatt(char(h5files(1).name), '/ECoG Array', 'Sampling Rate');
for j=2:length(h5files)
    % check that sampling rates match and set Fs
    x=h5readatt(char(h5files(j).name), '/ECoG Array', 'Sampling Rate');
    if x~= fs_orig
        error('preprocess.m: files have different sampling rates\n')
    end
end
clear x

% check number of channels
info=h5info(char(h5files(1).name),'/ECoG Array');
n=info.Dataspace.Size(1);
clear info

% loop over time
idx_save = 101; % savefile index %%% CHANGED THIS!!! PUT IT BACK AFTERWARDS!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
idx_load = 515-314; % 1; % load file index  %%% CHANGED THIS!!! PUT IT BACK AFTERWARDS!!!!!!!!!!!!!!!!!!!!!!!!
startTime = globalStart;
endTime = startTime + 12*60;


while startTime < globalEnd
    tic
    
    if idx_load > length(h5files)
        break
    end
    
    data=NaN(n, fs_orig*60*15); % make an array bigger than we need it to be
    idx_nan = 1; % index of first column to rewrite
    fprintf('Loading new data for file %03d...\n',idx_save)
    
    % load in data
    t = startTime;
    while t <= endTime && idx_load <= length(h5files)
        new=h5read(char(h5files(idx_load).name), '/ECoG Array');
        
        % zero mean each file
        new = bsxfun(@minus, new, mean(new,2));
        
        % trim beginning of if it was saved in the previous file
        times=h5read(char(h5files(idx_load).name), '/timestamp vector');
        s=find(times>=startTime,1,'first');
        new = new(:,s:end);
%         if idx_nan ==1
%             t_start = times(s+fs_orig*60); % fix this!!!
%         end
        
        % add new to data
        data(:, idx_nan:(idx_nan+size(new,2)-1) ) = new;
        idx_nan = idx_nan+size(new,2);
        
        % move bookmarks
        idx_load = 1+idx_load;
        if idx_load > length(h5files)
            break
        end
        times=h5read(char(h5files(idx_load).name), '/timestamp vector');
        t=times(1);
    end
    clear times new s
    
    
    % find new start time and end time
    startTime = startTime + 10*60;
    endTime = startTime + 12*60;
    while t>startTime % rewind idx_load to new startTime
        idx_load = idx_load-1;
        times=h5read(char(h5files(idx_load).name), '/timestamp vector');
        t=times(1);
    end
    clear t new 
    
    % drop grid channels and empty channels, make it 12 min long
    last = min((12*60*fs_orig), size(data,2));
    data = data(keep(1):keep(2), 1:last);
    
    
    % loop over channels
    fprintf('Filtering and downsampling each channel...\n')
    smallData = NaN(size(data,1), fs*60*12); % make matrix to hold 12 min of data
    for k=1:size(data,1)
        
        % bandpass filter 1-256 Hz
        data(k,:) = eegfilt(data(k,:), fs_orig, 1, fs/2);
        
        % downsample to 512 Hz
        smallData(k,:) = resample(data(k,:), fs, fs_orig);
        
        % notch filter 60, 120, 180, 240
        [b, a]=butter(3,2*[59 61]/fs,'stop'); %60hz
        smallData(k,:)=filtfilt(b, a, smallData(k,:)); %notch out at 60
        
        [b, a]=butter(3,2*[119 121]/fs,'stop'); %120hz
        smallData(k,:)=filtfilt(b, a, smallData(k,:)); %notch out at 120
        
        [b, a]=butter(3,2*[179 181]/fs,'stop'); %180hz
        smallData(k,:)=filtfilt(b, a, smallData(k,:)); %notch out at 180
        
        [b, a]=butter(3,2*[239 241]/fs,'stop'); %240hz
        smallData(k,:)=filtfilt(b, a, smallData(k,:)); %notch out at 240
          
    end
    smallData(:,isnan(smallData(1,:)))=[]; % delete any extra columns
    data = smallData;
    clear b a smallData
    
    % drop first and last minute (to reduce edge effects)
    last = min((11*60*fs), size(data,2));
    data=data(:, (1*60*fs+1):last);
    
    % common mean reference
    data=bsxfun(@minus,data,mean(data));
    
    % save .mat file (now 10 min long and preprocessed)
    fprintf('Saving file...\n')
    save(sprintf('mat_allChannels/%s_all_%03d.mat',prefix,idx_save),'data', 'fs'); %, 't_start');
    clear t_start
    idx_save=1+idx_save;
    
    
    toc
    
    
    
end
clear data idx k h5files keep n

fprintf('All done!\n')

end