function preprocessWrapper(montageFile,folderList,...
    newFs, langGridFlag, ekgFlag, medianFlag, preprocessfolder,Nminperfile,...
    Nmintrim,Flow,Fhigh)
%% preprocessWrapper.m

%% CHANGE LOG

% Created 03/29/2016
%   - run a whole patient's worth of files through preprocess.m, then 
%     reorganize using groupData.m
%   - pulls out channel names and keep indices from given montage
%     spreadsheet
%   by KHPD
%
% 20160512 KHPD - edited to handle electrode montages stored as .mat files
%

%%

% fill in missing parameters
n=nargin;
if n<1 || isempty(montageFile)
    error('montageFile is a required parameter')
end
if n<11
    Fhigh = []; % let preprocess.m set default
end
if n<10
    Flow = []; % let preprocess.m set default
end
if n<9
    Nmintrim = []; % let preprocess.m set default
end
if n<8
    Nminperfile = []; % let preprocess.m set default
end
if n<7 || isempty(preprocessfolder)
    preprocessfolder = 'mat_allChannels'; % groupData needs to know this
end
if n<6 
    medianFlag=[]; % let preprocess.m set default
end
if n<5 || isempty(ekgFlag) % these ones need defaults to use here
    ekgFlag=0;
end
if n<4 || isempty(langGridFlag)
    langGridFlag = 0;
end
if n<3 || isempty(newFs)
    newFs = 512;
end
if n<2 || isempty(folderList)
    folderList = dir('EC*');
    folderList = folderList([folderList(:).isdir]);
end
clear n

% read in montageSheet and get variable names
[eNames, eNums] = getElectrodeNames(montageFile);

% figure out 'keep' indices for preprocess.m
channelInds = getKeepRange(eNames, langGridFlag, ekgFlag);

% loop over folders
for j=1:length(folderList)
    fprintf('processing folder %s\n',char(folderList(j).name))
    
    % run preprocess.m
    tic
    preprocess(channelInds,newFs,[char(folderList(j).name),'/'],preprocessfolder,...
        Nminperfile,Nmintrim,Flow,Fhigh,medianFlag)
    toc
    
    % run groupData.m
    tic
    groupData(char(folderList(j).name),preprocessfolder,eNames,eNums,channelInds)
    toc
end

end

%% HELPER FUNCTIONS

function channelInds = getKeepRange(eNames, langGridFlag, ekgFlag)
% Assumes that channel organization is
%       Language Grid
%       Other neural channels
%       EKG channels
%       Empty channels

if langGridFlag == 1
    first = 1;
else
    first = max(find(~strcmp(eNames,'L64Contact'),1,'first'),find(~strcmp(eNames,'64'),1,'first'));
end

if ekgFlag
    last = size(eNames,1);
else
    last = find(~strcmp(eNames,'EKG'),1,'last');
end

channelInds = [first last];

end

function [eNames, eNums] = getElectrodeNames(montageFile)
% figure out if montage saved as spreadsheet or .mat
[~, ~, ext]=fileparts(montageFile);

if strcmp(ext, '.xlsx')
    % call xlsread to get the spreadsheet as text
    [~, txt, ~] = xlsread(montageFile);
elseif strcmp(ext, '.mat')
    % then it lives in the 'anatomy' variable
    load(montageFile,'anatomy');
    txt=anatomy;
    clear anatomy
else
    error('montageFile must be .xlsx or .mat')
end

% initialize arrays to hold electrode info
n=size(txt,1);
eNums=zeros(n,1);
eNames=cell(n,1);

% loop over electrodes
for j=1:n
    % remove extra words and add spaces
    new = regexprep(char(txt(j,2)), 'Electrode', ' ');
    new = regexprep(new, 'Strip', ' ');
    new = regexprep(new, 'Depth', ' ');
    new = regexprep(new, 'Grid', ' ');
    new = regexprep(new, 'EKG', 'EKG ');
    
    % break out the string and number portions
    A = textscan(new, '%s %d', 'MultipleDelimsAsOne', 1);
    
    % save them to the lists
    eNames(j) = A{1};
    if ~isempty(A{2}) % if this channel doesn't have a # it will stay 0
        eNums(j) = double(A{2});
    end
end

end