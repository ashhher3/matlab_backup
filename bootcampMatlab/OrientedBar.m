function OrientedBar()
% Models a population of neurons that respond to oriented bars. 
% Simulates data, calculates the response of neurons with different
% preferred orientations, makes a raster plot. No inputs or outputs. Uses
% raster.m from NS Bootcamp files and the helper functions in this file.
%
% 9/16/2015 KD

data=makeData(); % generate input data
prefs=0:pi/4:3*pi/4; % generate a set of preferred orientations
m=100; % max firing rate
s=0.3; % std. dev. for the gaussian tuning curves

pop=zeros(length(prefs),length(data)); % initialize matrix to hold responses
% compute the firing rates for each of the preferred orientations
for i=1:length(prefs)
    pop(i,:)=orientedNeuron(prefs(i),data,m,s);   
end

if length(prefs)==4 % if there's 4 neurons in the pop, use basic plot below
    plot4Neur(data,pop)
end

end


function firingRates=orientedNeuron(pref,stim,m,s)
% Models a single neuron with Gaussian tuning curve
% INPUTS: pref - the preferred orientation of the neuron,
%                given as a number in [0, pi]
%         stim - the stimulus, ie a moving oriented bar,
%         m,s  - mean and std dev to use in Gaussian
%                given as vector of orientations over time (in [0,pi])
% RETURNS: firingRates - the response of a neuron, as a vector of firing
%                        rates over time

n=size(stim,2); % how many columns in stim? these are our time bins
firingRates=zeros(1,n); % initialize FR

% calculate alternate pref for wrapping around [0 pi]
if pref < (pi/2)
    altpref=pref+pi;
else
    altpref=pref-pi;
end

for i=1:n % for each time point
    % compute difference between stimulus and preferred orientation,
    % accounting for the fact that 0 and pi are the same orientation
    diff=min(abs(pref-stim(i)),abs(altpref-stim(i)));
    firingRates(i)=gaussTuning(diff,m,s); % compute FR
end

end

function y=gaussTuning(x,m,s)
y=m*exp(-1*(x^2)/(2*s^2)); % from the wikipedia page for gaussian function
end

function data=makeData()
steps=-1*pi/2:pi/200:pi/2; % make a bunch of tiny steps over one cycle of sine
data=(pi/2)*sin(steps)+(pi/2); % scale so it fits in [0 pi]

% note there's no reason why your input should be sine shaped and not
% linear, I just thought that might be more interesting. I think it works
% out to be a bar that is rotating at a non-constant speed this way.
end

function plot4Neur(data,pop)
% plot the population response
figure(10) % maybe you're doing other stuff in figures 1-9 so leave them alone
clf % clear anything that's already there
hold on % draw everything in this figure

x=1:length(data); % this is just index numbers for our data

subplot(2,1,1) % plot the stimulus up top
plot(x,data,'k') % plot stimulus in black
ylim([0 pi])
xlim([0 201])
ylabel('stimulus orientation (rad)')

subplot(2,1,2) % plot the neurons on the bottom in different colors
plot(x,pop(1,:),'r',x,pop(2,:),'g',x,pop(3,:),'c',x,pop(4,:),'b')
xlim([0 201])
ylabel('percent of maximum firing rate')
xlabel('time bin #')

end
