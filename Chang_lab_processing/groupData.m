function groupData(ofc_idx, am_idx, hp_idx, acc_idx)
% make dirs to store .mat files
mkdir('mat_AmHpOFC')

% get the prefix for this set of files
[~,prefix,~] = fileparts(pwd);


% get list of the files we just made
shortMatFiles=dir('mat_allChannels/*.mat');

% loop over all channel .mat files
for j=1:3:length(shortMatFiles)
    
    fprintf('Loading mat_allChannels/%s_all_%03d, %03d, %03d...\n',prefix, j, j+1, j+2)
    
    % load in three files
    d2=[]; d3=[];   
    d1=importdata(sprintf('mat_allChannels/%s_all_%03d.mat',prefix,j));
    try
        d2=importdata(sprintf('mat_allChannels/%s_all_%03d.mat',prefix,j+1));
        try
            d3=importdata(sprintf('mat_allChannels/%s_all_%03d.mat',prefix,j+2));
        catch
            fprintf('mat_AmHpOFC/%s_%03d.mat will only be 20 min long.\n',prefix,ceil(j/3))
        end
    catch
        fprintf('mat_AmHpOFC/%s_%03d.mat will only be 10 min long.\n', prefix, ceil(j/3))
    end
    
    
    % concatenate, keeping only Am, Hp, OFC
    data = [d1.data, d2.data, d3.data];
    fs = d1.fs;
    clear d1 d2 d3
    
    am=data(am_idx(1):am_idx(2),:);
    hp=data(hp_idx(1):hp_idx(2),:);
    ofc=data(ofc_idx(1):ofc_idx(2),:);
    acc=data(acc_idx(1):acc_idx(2),:);
    clear data
    
    % save .mat file
    fprintf('Saving mat_AmHpOFCACC/%s_%03d.mat\n',prefix, ceil(j/3))
    
    save(sprintf('mat_AmHpOFCACC/%s_%03d.mat',prefix,ceil(j/3)), 'am', 'hp', 'ofc', 'acc','fs');
    clear am hp ofc
end

fprintf('All done!\n')

end