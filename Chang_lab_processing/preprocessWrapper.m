function [eNames, eNums]= preprocessWrapper(montageSheet, folderList, newFs, ...
    langGridFlag, ekgFlag, medianFlag)

%% preprocessWrapper.m

%% CHANGE LOG

% Created 03/29/2016
%   - run a whole patient's worth of files through preprocess.m, then 
%     reorganize using groupData. Automate as much as possible.
%   by KHPD
%%

% fill in missing parameters
n=nargin;
if n<1 || isempty(montageSheet)
    error('montageSheet is a required parameter')
end
if n<6 || isempty(medianFlag)
    medianFlag=0;
end
if n<5 || isempty(ekgFlag)
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
[eNames, eNums] = getElectrodeNames(montageSheet);

% figure out 'keep' indices for preprocess.m

% loop over folders
for j=1:length(folderList)
    % run preprocess.m

    % run groupData.m
end

end

%% HELPER FUNCTIONS

function [eNames, eNums] = getElectrodeNames(montageSheet)
% call xlsread to get the spreadsheet as text
[~, txt, ~] = xlsread(montageSheet);

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