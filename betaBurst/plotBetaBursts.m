function [fig_Bpwr, fig_Rast] = plotBetaBursts(area1, bad1, area2, bad2)

% set some variables about drawing things nicely
set(0,'units','pixels');
pixSS = get(0,'screensize');
ctxC = [0.94 0.66 0.02];
bgC = [0.29 0.62 0.96];

% get number of all good channels
a1Num=[];
for i=1:length(area1.contact)
    if ~any(bad1==i)
        a1Num = [a1Num, i];
    end
end
a2Num=[];
for i=1:length(area2.contact)
    if ~any(bad2==i)
        a2Num = [a2Num, i];
    end
end

%% PLOT BETA RASTER OF BETA BURSTS
figHt = 25*(length(a1Num)+length(a2Num))+100; % Set figure height based on number of thresh
figPos = pixSS(4)-(0.2*pixSS(4)+figHt); % Set position based on figure size
fig_Rast= figure('Position',[600,figPos,500,figHt],'Color',[1 1 1]); % Make new figure

hold on; % Plot all data on same axes 

channelLabels = cell(0);
for i = 1:length(a1Num) % plot the ecog channels
    j=a1Num(i);
    step=1/area1.FsB(j);
	train = area1.contact(j).burst; 
	allTimes = [0:step:length(train)*step-step]; % make timestamps
	bTimes = allTimes(logical(train)); % Grab the timestamps when a burst happened
	times = ones(1,numel(bTimes)); % Generate dummy values so plot will look like a raster
	times = times*i; % Transform dummy values so each threshold gets its own line
	plot(bTimes,times,'bs','markersize',1.5,'markerfacecolor',ctxC,'markeredgecolor',ctxC)
    channelLabels = [channelLabels, sprintf('ECoG %d',j)];
end

for i = 1:length(a2Num) % plot the lfp channels
    j=a2Num(i);
    step=1/area2.FsB(j);
	train = area2.contact(j).burst; 
	allTimes = [0:step:length(train)*step-step]; % make timestamps
	bTimes = allTimes(logical(train)); % Grab the timestamps when a spike happened
	times = ones(1,numel(bTimes)); % Generate dummy values so plot will look like a raster
	times = times*i+length(a1Num); % Transform dummy values so each threshold gets its own line
	plot(bTimes,times,'bs','markersize',1.5,'markerfacecolor',bgC,'markeredgecolor',bgC)
    channelLabels = [channelLabels, sprintf('LFP %d',j)]; 
end
hold off



% Set labels and titles for graph
set(gca,'ylim',[0.5 length(a1Num)+length(a2Num)+0.5]); % Scale the y axis
set(gca,'ytick',[1:length(a1Num)+length(a2Num)]); % Set the y axis tick marks
set(gca,'YTickLabel',channelLabels); % Set the y axis tick labels
ylabel('Channel'); % Set the y axis label

set(gca,'xlim',[allTimes(1)-1 allTimes(numel(allTimes))+1]); % Scale the x axis
xlabel('Time (s)'); % Set the x axis label

title('Beta Burst Raster','FontWeight','bold','FontSize',12); % Set the plot title
%set(gcf,'name','Beta Burst Raster','numbertitle','off')  % Set the title of the figure window

%% PLOT BETA POWER OVER THE WHOLE RECORDING
fig_Bpwr = figure('Position',[600,figPos,500,figHt],'Color',[1 1 1]); % Make new figure
hold on; % Plot all data on same axes 

ampScale=10e6; % set this to 1/ roughly the amplitude beta power on ecog channels
for i = 1:length(a1Num) % plot the ecog channels
    j=a1Num(i);
    step=1/area1.FsB(j);
	trace = ampScale*area1.contact(j).betaPwr; 
	time = [0:step:length(trace)*step-step]; % make timestamps
    % want to put the traces so their median power values are evenly spaced & labels are by medians
    shift = i - ampScale*area1.contact(j).medB;
	p = plot(time,trace+shift);
    plot(time, ones(size(time))*ampScale*area1.contact(j).medB+shift,'Color',get(p,'Color'),'LineStyle','--');
end

ampScale=10e6; % set this to 1/ roughly the amplitude beta power on lfp channels
for i = 1:length(a2Num) % plot the lfp channels
    j=a2Num(i);
    step=1/area2.FsB(j);
    trace = ampScale*area2.contact(j).betaPwr;
    time = [0:step:length(trace)*step-step]; % make timestamps
    % want to put the traces so their median power values are evenly spaced & labels are by medians
    shift = (i+length(a1Num)) - ampScale*area2.contact(j).medB;
	p = plot(time,trace+shift);
    plot(time, ones(size(time))*ampScale*area2.contact(j).medB+shift,'Color',get(p,'Color'),'LineStyle','--');
end
hold off



% Set labels and titles for graph
set(gca,'ylim',[0 length(a1Num)+length(a2Num)+1]); % Scale the y axis
set(gca,'ytick',[1:length(a1Num)+length(a2Num)]); % Set the y axis tick marks
set(gca,'YTickLabel',channelLabels); % Set the y axis tick labels
ylabel('Channel'); % Set the y axis label

set(gca,'xlim',[allTimes(1) allTimes(end)]); % Scale the x axis
xlabel('Time (s)'); % Set the x axis label

title('Beta Power','FontWeight','bold','FontSize',12); % Set the plot title
%set(gcf,'name','Beta Burst Raster','numbertitle','off')  % Set the title of the figure window

end