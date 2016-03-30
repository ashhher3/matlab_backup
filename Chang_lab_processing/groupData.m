function groupData(filedir,preprocessfolder,eNames,eNums,channelInds)

fprintf('groupData.m: Initializing...\n')

% get list of the files made by preprocess.m
shortMatFiles=dir(sprintf('%s/%s/*.mat',filedir,preprocessfolder));


% make dirs to store regrouped files
mkdir(sprintf('%s/mat_oneStruct',filedir));

% loop over all channel .mat files
for j=1:length(shortMatFiles)
    
    fprintf('\tLoading %s/%s/%s_all_%03d...\n',filedir,preprocessfolder, filedir, j)
    
    % load in data
    d1=importdata(sprintf('%s/%s/%s_all_%03d.mat',filedir,preprocessfolder,filedir,j));
 
    % reorganize data
    data.signal = d1.data;
    data.fs = d1.fs;
    data.startTime = d1.startTime;
    clear d1
    
    data.channelNames = eNames(channelInds(1):channelInds(2));
    data.channelNums = eNums(channelInds(1):channelInds(2));
    
    % save .mat file
    fprintf('\tSaving mat_oneStruct/%s_%03d.mat\n',filedir, j)
    
    save(sprintf('%s/mat_oneStruct/%s_%03d.mat',filedir,filedir,j), 'data');
    
end

fprintf('\tAll done!\n')

end