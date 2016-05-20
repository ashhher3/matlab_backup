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