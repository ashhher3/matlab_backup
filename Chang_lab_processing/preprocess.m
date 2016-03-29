function preprocess(channelInds,filedir)

%-------------------------------------------------------------------------%
% Revised: 03/28/16 (JGM)
%   -functionized
% Revised: 03/25/16 (JGM)
%   -only load the channels that are requested with the keep vector, rather
%   than loading all and then deleting some
%   -added argument filedir so that this file can be run from other
%   directories
%   -cleaned up, rearranged
% Created: 02/??/16
%   by KD
%-------------------------------------------------------------------------%

fprintf('Initializing...\n')

% PARAMETERS TO SET (might want to make these arguments to this fxn)
fs = 512;
Nminperfile = 12;
Nmintrim = 1;
Flow = 1;
Fhigh = fs/2;
foldername = 'mat_jgmChannels';


% get the prefix for this set of files
if isempty(filedir)
    [~,prefix,~] = fileparts(pwd);
else
    [~,prefix,~] = fileparts(filedir(1:end-1));
end

% make dirs to store .mat files
mkdir([filedir,foldername])

% get list of .h5 files
h5files = dir([filedir,'*.h5']);

% get start and end times
globalStart = h5read([filedir,h5files(1).name], '/timestamp vector',1,1);
info = h5info([filedir,h5files(end).name],'/timestamp vector');
globalEnd = h5read([filedir,h5files(end).name], '/timestamp vector',info.Dataspace.Size,1);


% loop over h5 files
fprintf('Checking the sampling rate...\n')
fs_orig = h5readatt([filedir,h5files(1).name], '/ECoG Array', 'Sampling Rate');
for j = 2:length(h5files)
    % check that sampling rates match and set Fs
    x = h5readatt([filedir,h5files(j).name], '/ECoG Array', 'Sampling Rate');
    if x ~= fs_orig
        error('preprocess.m: files have different sampling rates\n')
    end
end
clear x

% check number of channels
%%%info = h5info([filedir,h5files(1).name],'/ECoG Array');
%%%Nchannels = info.Dataspace.Size(1);
%%%Nsamples = info.Dataspace.Size(2);
%%%clear info


% init indices
idx_save = 1; % savefile index
idx_load = 1; % load file index
startTime = globalStart;
endTime = startTime + Nminperfile*60;

% malloc
data = NaN(diff(channelInds)+1, fs_orig*60*ceil(Nminperfile*1.25));
smallData = NaN(diff(channelInds)+1, fs*60*Nminperfile);

% loop over time
while (startTime < globalEnd) && (idx_load <= length(h5files))
    tic
    
    % init indices
    fprintf('Loading new data for file %03d...\n',idx_save)
    idx_nan = 1; % index of first column to rewrite
    t0 = startTime;
    
    % loop over h5 files until Nminperfile minutes of data or no more files
    while t0 <= endTime && idx_load <= length(h5files)
        
        % add new to data
        sampleInds = getNewDataInds([filedir,h5files(idx_load).name],startTime);
        idx_nanF = idx_nan + diff(sampleInds);
        data(:, idx_nan:idx_nanF) = getNextH5array(...
            channelInds,sampleInds,[filedir,h5files(idx_load).name]);
        idx_nan = idx_nanF + 1;
        
        % move bookmarks
        idx_load = idx_load + 1;
        if idx_load <= length(h5files)
            t0 = h5read([filedir,h5files(idx_load).name], '/timestamp vector',1,1);
        end
    end
    
    % set final index at Nminperfile mark
    last = min((Nminperfile*60*fs_orig), size(data,2));
    
    % filter each channel consecutively
    fprintf('Filtering and downsampling each channel...\n')
    for k = 1:size(data,1)
        smallData(k,:) = standardFilters(data(k,1:last),...
            fs_orig, fs, Flow, Fhigh);
    end
    
    % trim, save, increment the file counter
    trimAndSave(smallData,startTime,Nminperfile,Nmintrim,fs,...
        filedir,foldername,prefix,idx_save);
    idx_save = 1 + idx_save;
    
    
    % find *next* start and end times
    startTime = startTime + (Nminperfile - 2*Nmintrim)*60;
    endTime = startTime + Nminperfile*60;
    while t0 > startTime % rewind idx_load to new startTime
        fprintf('rewinding...\n');
        idx_load = idx_load - 1;
        t0 = h5read([filedir,h5files(idx_load).name], '/timestamp vector',1,1);
    end
    toc

end
%clear data idx k h5files keep n
fprintf('All done!\n')

end
%-------------------------------------------------------------------------%
%-------------------------------------------------------------------------%



%-------------------------------------------------------------------------%
function sampleInds = getNewDataInds(filename,startTime)

times = h5read(filename, '/timestamp vector');
i0 = find(times>=startTime,1,'first');
if isempty(i0)
    fprintf('\nwarning!!! time skips -- jgm\n\n');
    sampleInds(1) = 1;
else
    sampleInds(1) = i0;
end
info = h5info(filename,'/ECoG Array');
sampleInds(2) = info.Dataspace.Size(2);

end
%-------------------------------------------------------------------------%

%-------------------------------------------------------------------------%
function new = getNextH5array(channelInds,sampleInds,filename)

% for this file, load in all samples of just the KEEP channels
new = h5read(filename, '/ECoG Array',...
    [channelInds(1),sampleInds(1)],[diff(channelInds)+1,diff(sampleInds)+1]);

% zero mean each file
new = bsxfun(@minus, new, mean(new,2));


end
%-------------------------------------------------------------------------%

%-------------------------------------------------------------------------%
function channelData = standardFilters(channelData,fs_orig,fs,Flow,Fhigh)

% bandpass filter 1-256 Hz
channelData = eegfilt(channelData, fs_orig, Flow, Fhigh,...
    0, 3*fix(fs_orig/Flow),0,'fir1');

% downsample to 512 Hz
channelData = resample(channelData, fs, fs_orig);

% notch filter 60, 120, 180, 240
[b, a]=butter(3,2*[59 61]/fs,'stop');           % 60hz
channelData = filtfilt(b, a, channelData);      % notch out at 60

[b, a]=butter(3,2*[119 121]/fs,'stop');         % 120hz
channelData = filtfilt(b, a, channelData);      % notch out at 120

[b, a]=butter(3,2*[179 181]/fs,'stop');         % 180hz
channelData = filtfilt(b, a, channelData);      % notch out at 180

[b, a]=butter(3,2*[239 241]/fs,'stop');         % 240hz
channelData = filtfilt(b, a, channelData);      % notch out at 240


end
%-------------------------------------------------------------------------%


%-------------------------------------------------------------------------%
function trimAndSave(smallData,startTime,Nminperfile,Nmintrim,fs,...
    filedir,foldername,prefix,idx_save)

%%%% presumably this is for screwed up channelData.....
% delete any extra columns
smallData(:,isnan(smallData(1,:)))=[];
%%%%

% drop first and last minute (to reduce edge effects)
last = min((Nminperfile-Nmintrim)*60*fs, size(smallData,2));
smallData = smallData(:, (Nmintrim*60*fs+1):last);

% common mean reference
smallData = bsxfun(@minus,smallData,mean(smallData));

% save .mat file (now Nminperfile long and preprocessed)
fprintf('Saving file...\n')
filename = sprintf('%s%s/%s_all_%03d.mat',filedir,foldername,prefix,idx_save);
save(filename,'smallData', 'fs','startTime');

end
%-------------------------------------------------------------------------%











