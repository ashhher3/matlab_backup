%% get data for 24 hour period

% ASSUMING WE START IN DATA FOLDER

% make dir to save to
mkdir('processed_fewChannelsAllTime');

% set channels we're interested in and number of files to look at
acc_chan = 4;
am_chan = 10;
hp_chan = 10;
ofc_chan = 19;
lastfile = (22*60)/30; % language mapping happens during last 2 hours

% set up mapping from quantities of interest to matrix
vars = {'ACC beta power', 'Am beta power', 'Hp beta power', ...
    'OFC beta power', 'Acc-Am coherence', 'Acc-Hp coherence', ...
    'Acc-OFC coherence', 'Am-Hp coherence', 'Am-OFC coherence', ...
    'Hp-OFC coherence'};

x = NaN(size(vars,2), lastfile*30*60/10);
nan_idx = 1;

% loop over files
for file = 1:lastfile
    % load next file
    load(sprintf('mat_AmHpOFCACC/EC108_2fd7f1ac_%03d.mat',file))
    
    % restructure data
    acc_old = acc;
    acc = struct();
    acc.signal = acc_old(acc_chan,:);
    acc.fs = fs;
    clear acc_old
    
    am_old = am;
    am = struct();
    am.signal = am_old(am_chan,:);
    am.fs = fs;
    clear am_old
    
    hp_old = hp;
    hp = struct();
    hp.signal = hp_old(hp_chan,:);
    hp.fs = fs;
    clear hp_old
    
    ofc_old = ofc;
    ofc = struct();
    ofc.signal = ofc_old(ofc_chan,:);
    ofc.fs = fs;
    clear ofc_old
    
    % compute beta for each channel of interest
    acc=getBetaPwr(acc);
    am=getBetaPwr(am);
    hp=getBetaPwr(hp);
    ofc=getBetaPwr(ofc);
    
    % compute 10s beta power for each channel of interest
    accB = getPwr10Sec(acc.betaPwr,fs);
    n = length(accB);
    x(1,nan_idx:nan_idx + n -1) = accB;
    x(2,nan_idx:nan_idx + n -1) = getPwr10Sec(am.betaPwr,fs);
    x(3,nan_idx:nan_idx + n -1) = getPwr10Sec(hp.betaPwr,fs);
    x(4,nan_idx:nan_idx + n -1) = getPwr10Sec(ofc.betaPwr,fs);
    
    % compute 10s coherences between channels of interest
    x(5,nan_idx:nan_idx + n -1) = getCoh10Sec(acc.beta,am.beta,fs);
    x(6,nan_idx:nan_idx + n -1) = getCoh10Sec(acc.beta,hp.beta,fs);
    x(7,nan_idx:nan_idx + n -1) = getCoh10Sec(acc.beta,ofc.beta,fs);
    x(8,nan_idx:nan_idx + n -1) = getCoh10Sec(am.beta,hp.beta,fs);
    x(9,nan_idx:nan_idx + n -1) = getCoh10Sec(am.beta,ofc.beta,fs);
    x(10,nan_idx:nan_idx + n -1) = getCoh10Sec(hp.beta,ofc.beta,fs);
    
    nan_idx = nan_idx + n - 1;
    
    fprintf('Done with file %03d\n',file)
    
end

x = x(:, 1:end-1);

save('processed_fewChannelsAllTime/beta.mat');

%% clean workspace

clear accB lastfile n acc am hp ofc

%% ACC-Am coh vs ACC-Hp coh
a = 1;
b = 3;

for h=0.5:0.5:22
    start = 1; %(h-0.5)*60*60/10;
    last = h*60*60/10;
    
    last = min(last, size(x,2));
    
    indep = x(a,start:last)';
    dep = x(b,start:last)';
    
    fig=figure(1);
    set(fig,'Color',[1 1 1])
    clf
    hold on
    c=linspace(0,10*length(indep),length(indep));
    scatter(indep,dep,[],c)
    colormap('jet')
    cb = colorbar;
    cb.Label.String = 'Time since start of recording (s)';
% 
%     xlim([2 22])
%     ylim([10 60])
    
    xlabel(vars(a))
    ylabel(vars(b))
    title(sprintf('%s vs %s over %0.1f hours', char(vars(a)), char(vars(b)), h))
    
    hold off
    pause(2)
end

clear a b dep indep l xplot yplot