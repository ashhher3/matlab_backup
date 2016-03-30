function preprocess(channelInds,fs,filedir,savefolder,Nminperfile,...
    Nmintrim,Flow,Fhigh,medianFlag)
%% preprocess.m
% Takes HDF5 files with neural data stored in '/ECoG Array' and timestamps 
% stored in '/timestamp vector', preprocesses, and stores as .mat files.
%
% INPUTS:
%   channelInds - 2 vector containing first and last channel to keep
%   fs - new sampling rate to downsample to
%   filedir - OPTIONAL, name of the directory containing HDF5 files to
%           preprocess. Defaults to current directory.
%   savefolder - OPTIONAL, name of the subdirectory where processed data
%           will be saved. Defaults to 'mat_allChannels'
%   Nminperfile - OPTIONAL, number of minutes of data to preprocess at
%           once. Defaults to 12.
%   Nmintrim - OPTIONAL, number of minutes to trim off at the beginning and
%           ending after preprocessing, to avoid edge effects. Saved .mat 
%           files will be Nminperfile - 2*Nmintrim long. Trimmed sections
%           overlap, so that final files can be concatenated with no loss.
%   Flow - OPTIONAL, low frequency for the bandpass filter
%   Fhigh - OPTIONAL, high frequency for the bandpass filter
%   medianFlag - OPTIONAL, if this is set to 1, performs common median 
%           reference instead of common mean reference
%
% OUTPUTS: Saves a sequence of .mat files. Filenames are of the form 
%   'savefolder/filedir_all_%03d.mat'. Each file contains 3 variables:
%   'data', 'fs', and 'startTime'.
%
% PREPROCESSING STEPS:
%   1 - Gather Nminperfile minutes of data by sequentially loading HDF5
%       files (zero-meaning each file as we load it in)
%   2 - Bandpass filter each channel, at Flow-Fhigh Hz
%   3 - Downsample each channel to fs Hz
%   4 - Notch filter at 60Hz and harmonics
%   5 - Trim the first Nmintrim minutes and last Nmintrim minutes
%   6 - Perform common mean reference at each time point
%

%% CHANGE LOG
%-------------------------------------------------------------------------%
% Revised: 03/29/16 (KHPD)
%   -tweaked parameters and added some input checking
%   -changed some of the print statements
%   -removed the tic toc
%   -reverted the name of the saved data to 'data' from 'smallData'
% Revised: 03/28/16 (JGM)
%   -functionized
% Revised: 03/25/16 (JGM)
%   -only load the channels that are requested with the keep vector, rather
%   than loading all and then deleting some
%   -added argument filedir so that this file can be run from other
%   directories
%   -cleaned up, rearranged
% Created: 03/07/16
%   -first version of preprocessing script
%   -turn h5 files of unknown (high) sampling rate and unknown length into 
%    mat files sampled at 512 Hz and 10 min long
%   -take overlapping windows and trim edges to avoid edge effects
%   by KHPD
%-------------------------------------------------------------------------%
%%

fprintf('preprocess.m: Initializing...\n')

% fill in missing parameters
n=nargin;
if n<2
    error('channelInds and fs are required inputs\n')
end
if n<9 || isempty(medianFlag)
    medianFlag = 0;
end
if n<8 || isempty(Fhigh)
    Fhigh = fs/2;
end
if n<7 || isempty(Flow)
    Flow = 1;
end
if n<6 || isempty(Nmintrim)
    Nmintrim = 1;
end
if n<5 || isempty(Nminperfile)
    Nminperfile = 12;
end
if n<4 || isempty(savefolder)
    savefolder = 'mat_allChannels';
end
if n<3
    filedir = [];
end
clear n

% check that parameters make sense
if 2*Fhigh>fs
    error('Fhigh must be <= fs/2\n');
end
if ~(Flow<Fhigh)
    error('Flow must be < Fhigh\n');
end
if ~(Nminperfile-2*Nmintrim>0)
    error(['Cannot trim 2x%d min off of a %d min array,'...
        ' try smaller Nmintrim or larger Nminperfile\n'],Nmintrim,Nminperfile)
end


% get the prefix for this set of files
if isempty(filedir)
    [~,prefix,~] = fileparts(pwd);
else
    [~,prefix,~] = fileparts(filedir(1:end-1));
end

% make dirs to store .mat files
mkdir([filedir,savefolder])

% get list of .h5 files
h5files = dir([filedir,'*.h5']);

% get start and end times
globalStart = h5read([filedir,h5files(1).name], '/timestamp vector', 1,1);
info = h5info([filedir,h5files(end).name],'/timestamp vector');
globalEnd = h5read([filedir,h5files(end).name], '/timestamp vector',info.Dataspace.Size,1,1);


% loop over h5 files
fprintf('\tChecking the sampling rate...\n')
fs_orig = h5readatt([filedir,h5files(1).name], '/ECoG Array', 'Sampling Rate');
for j = 2:length(h5files)
    % check that sampling rates match and set Fs
    x = h5readatt([filedir,h5files(j).name], '/ECoG Array', 'Sampling Rate');
    if x ~= fs_orig
        error('preprocess.m: h5 files have different sampling rates\n')
    end
end
clear x


% init indices
idx_save = 1; % savefile index
idx_load = 1; % load file index
startTime = globalStart;
endTime = startTime + Nminperfile*60;

% malloc
try
    data = NaN(diff(channelInds)+1, fs_orig*60*ceil(Nminperfile*1.25));
    smallData = NaN(diff(channelInds)+1, fs*60*Nminperfile);
catch
    error('Requested matrices are too large, try smaller fs or smaller Nminperfile\n')
end


% loop over time
while (startTime < globalEnd) && (idx_load <= length(h5files))
    
    
    % init indices
    fprintf('\tLoading new data for file %03d...\n',idx_save)
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
    fprintf('\tFiltering and downsampling each channel...\n')
    for k = 1:size(data,1)
        smallData(k,:) = standardFilters(data(k,1:last),...
            fs_orig, fs, Flow, Fhigh);
    end
    
    % trim, save, increment the file counter
    trimAndSave(smallData,startTime,Nminperfile,Nmintrim,fs,...
        filedir,savefolder,prefix,idx_save, medianFlag);
    idx_save = 1 + idx_save;
    
    
    % find *next* start and end times
    startTime = startTime + (Nminperfile - 2*Nmintrim)*60;
    endTime = startTime + Nminperfile*60;
    while t0 > startTime % rewind idx_load to new startTime
        fprintf('\trewinding...\n');
        idx_load = idx_load - 1;
        t0 = h5read([filedir,h5files(idx_load).name], '/timestamp vector',1,1);
    end
    

end
%clear data idx k h5files keep n
fprintf('\tAll done!\n')

end

%-------------------------------------------------------------------------%
%-------------------------------------------------------------------------%


%% HELPER FUNCTIONS
%-------------------------------------------------------------------------%
function sampleInds = getNewDataInds(filename,startTime)

times = h5read(filename, '/timestamp vector');
i0 = find(times>=startTime,1,'first');
if isempty(i0)
    fprintf('\nwarning!!! time skips in file %s\n\n',filename);  %hey joe i removed your '--jgm' because i wanted to add the filename and i figured it'd be long
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
function trimAndSave(data,startTime,Nminperfile,Nmintrim,fs,...
    filedir,savefolder,prefix,idx_save,medianFlag)

%%%% presumably this is for screwed up channelData.....
% delete any extra columns
data(:,isnan(data(1,:)))=[];
%%%%

% drop first and last minute (to reduce edge effects)
last = min((Nminperfile-Nmintrim)*60*fs, size(data,2));
data = data(:, (Nmintrim*60*fs+1):last);

% common mean reference
if medianFlag==1
    data = bsxfun(@minus,data,median(data));
else
    data = bsxfun(@minus,data,mean(data));
end

% save .mat file (now Nminperfile long and preprocessed)
fprintf('\tSaving file %03d...\n',idx_save)
filename = sprintf('%s%s/%s_all_%03d.mat',filedir,savefolder,prefix,idx_save);
save(filename,'data', 'fs','startTime');

end
%-------------------------------------------------------------------------%











