function raster(spTrains,thresholds)
% Takes a cell array of spike trains and a vector of threshold values and
% plots a raster of spikes for each threshold value.  There are no outputs.
%
% 8/3/2015 KD & NH

% Get the size of the screen so the plot will show up in a reasonable place
set(0,'units','pixels');
pixSS = get(0,'screensize');

figHt = 25*numel(thresholds)+100; % Set figure height based on number of thresh
figPos = pixSS(4)-(0.2*pixSS(4)+figHt); % Set position based on figure size
fig = figure('Position',[600,figPos,500,figHt]); % Make new figure

hold on; % Plot all data on same axes 

% Make sure thresholds are in ascending order, if not, flip
if thresholds(1) > thresholds(numel(thresholds))
	thresholds = fliplr(thresholds);
	spTrains = fliplr(spTrains);
end

for i = 1:size(spTrains,2) % Make plot for each threshold

	spTrain = spTrains{i}; % Use correct spike train
	allTimes = [1:numel(spTrain)]; % Generate array of timestamps for spike train
	spTimes = allTimes(logical(spTrain)); % Grab the timestamps when a spike happened
	times = ones(1,numel(spTimes)); % Generate dummy values so plot will look like a raster
	times = times*i; % Transform dummy values so each threshold gets its own line
	plot(spTimes,times,'bs','markersize',1.5,'markerfacecolor','b','markeredgecolor',[0 0 .7])
	
end

hold off

% Set labels and titles for graph
set(gca,'ylim',[0.5 numel(thresholds)+0.5]); % Scale the y axis
set(gca,'ytick',[1:numel(thresholds)]); % Set the y axis tick marks
set(gca,'YTickLabel',strread(num2str(thresholds),'%s')); % Set the y axis tick labels
ylabel('Threshold'); % Set the y axis label

set(gca,'xlim',[allTimes(1)-1 allTimes(numel(allTimes))+1]); % Scale the x axis
xlabel('Time'); % Set the x axis label

title('Raster for Each Threshold','FontWeight','bold','FontSize',12); % Set the plot title
set(gcf,'name','Your Integrate and Fire Model','numbertitle','off')  % Set the title of the figure window


end