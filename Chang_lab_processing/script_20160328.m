%% trim the extra minute off of mat files for 9ae...

% get list of mat files
matfiles=dir('*.mat');

% loop over mat files
for j=1:length(matfiles)
    
    if j~=102 && j~=103 && j~=114
        % load in next file
        load(char(matfiles(j).name))
        
        % trim off last minute
        data = data(:, 1:10*60*fs);
        
        % save
        save(char(matfiles(j).name), 'data', 'fs');
        
    end
    
    % print
    fprintf('done with %s\n',char(matfiles(j).name))
end