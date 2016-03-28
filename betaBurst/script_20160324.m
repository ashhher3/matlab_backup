%% get data for 24 hour period

% ASSUMING WE START IN DATA FOLDER

% ASSUMING WE ALREADY DID script_20160322 to get acc, am, hp, ofc

rand.idx = [116, 131, 137, 141, 145, 173];
rand.idx = rand.idx-64;
rand.vars = {'Lateral frontal beta power', 'Temporal pole beta power', ...
    'Middle Subtemporal beta power', 'Posterior Subtemporal beta power',...
    'Inferior Temporal beta power', 'Insula beta power'};

lastfile = (22*60)/10+2; % language mapping happens during last 2 hours

y = NaN(size(rand.vars,2), lastfile*10*60/10);
nan_idx = 1;

% loop over files
for file = 1:lastfile
    % load next file
    load(sprintf('mat_allChannels/EC108_2fd7f1ac_all_%03d.mat',file))
    
    % restructure data
    rand.signal = data(rand.idx,:);
    rand.fs = 512;
    fs = 512;
    
    % compute beta for each channel of interest
    rand = getBetaPwr(rand);

    % compute 10s beta power for each channel of interest
    first = getPwr10Sec(rand.betaPwr(1,:),fs);
    n = length(first);
    y(1,nan_idx:nan_idx + n -1) = first;
    for j=2:size(rand.signal,1)
        y(j,nan_idx:nan_idx + n -1) = getPwr10Sec(rand.betaPwr(j,:),fs);
    end
    
    
    nan_idx = nan_idx + n - 1;
    
    fprintf('Done with file %03d\n',file)
    
end

y = y(:, 1:end-1);

clear lastfile nan_idx file fs first n j data

save('processed_fewChannelsAllTime/beta_rand.mat');

%% clean workspace

clear accB lastfile n acc am hp ofc

%% ACC-Am coh vs ACC-Hp coh

fig=figure(1);

for a = 1:4
    for b=1:6
        set(fig, 'Color', [0 0 0])
        clf
        text(0.5, 0.5,sprintf('%s \n vs. \n %s', char(vars(a)), ...
            char(rand.vars(b))), 'Color', 'white', 'FontSize', 18, ...
            'Units', 'normalized', 'HorizontalAlignment', 'center')
        set(gca, 'Color', [0 0 0], 'XTick', [], 'YTick', [], ...
            'XColor', 'none', 'YColor', 'none')
        pause(3)
        for h=0.5:0.2:22
            start = max(1, fix((h-5)*60*60/10));
            last = fix(h*60*60/10);
            
            last = min(last, size(y,2));
            
            indep = x(a,start:last)';
            dep = y(b,start:last)';
            
            fig=figure(1);
            set(fig,'Color',[1 1 1])
            clf
            hold on
            c=linspace(0,10*length(indep),length(indep));
            scatter(indep,dep,[],c)
            colormap('cool')
            cb = colorbar;
            cb.Label.String = 'Time since start of recording (s)';%
            xlim([2 22])
            ylim([0 60])
            
            xlabel(vars(a))
            ylabel(rand.vars(b))
            title(sprintf('%s vs %s over %0.1f hours', char(vars(a)), char(rand.vars(b)), h))
            
            hold off
            pause(0.00001)
        end
    end
end

clear a b dep indep l xplot yplot

%% ACC-Am coh vs ACC-Hp coh

fig=figure(1);

for a = 1:4
    for b=a+1:4
        set(fig, 'Color', [0 0 0])
        clf
        text(0.5, 0.5,sprintf('%s \n vs. \n %s', char(vars(a)), ...
            char(vars(b))), 'Color', 'white', 'FontSize', 18, ...
            'Units', 'normalized', 'HorizontalAlignment', 'center')
        set(gca, 'Color', [0 0 0], 'XTick', [], 'YTick', [], ...
            'XColor', 'none', 'YColor', 'none')
        pause(3)
        for h=0.5:0.2:22
            start = max(1, fix((h-5)*60*60/10));
            last = fix(h*60*60/10);
            
            last = min(last, size(y,2));
            
            indep = x(a,start:last)';
            dep = x(b,start:last)';
            
            fig=figure(1);
            set(fig,'Color',[1 1 1])
            clf
            hold on
            c=linspace(0,10*length(indep),length(indep));
            scatter(indep,dep,[],c)
            colormap('cool')
            cb = colorbar;
            cb.Label.String = 'seconds';%
            xlim([2 22])
            ylim([0 60])
            
            xlabel(vars(a))
            ylabel(vars(b))
            title(sprintf('%s vs %s from %0.1f hr to %0.1f hr', char(vars(a)), char(vars(b)), max(0,h-5), h))
            
            hold off
            pause(0.0000001)
        end
    end
end

clear a b dep indep l xplot yplot

%%

for i=1:190
new = regexprep(char(txt(i,2)), 'Electrode', ' ');
new = regexprep(new, 'Strip', ' ');
new = regexprep(new, 'Depth', ' ');
new = regexprep(new, 'Grid', ' ');
A = textscan(new, '%s %d', 'MultipleDelimsAsOne', 1);
if ~isempty(A{2})
names(i) = A{1};
numbers(i) = double(A{2});
end
end