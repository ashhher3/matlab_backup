%% load data
clear all
plotDate = '20151021';
load(['EffortVSReward_MP_data_' plotDate '_1.mat'])

%% ...load individual trials
clear all
plotDate = '20160105';
top=cd('trials_1');
files = dir;
cd(top);
for ii=1:(numel(files)-2)
    tmp = load(sprintf('trials_1/DATA_%03.0f.mat',ii));
    DATA(ii) = tmp.(sprintf('DATA_%03.0f',ii));
end
%save(['EffortVSReward_MP_data_' plotDate '_1.mat'], 'DATA');
clear ii tmp files top

%% Make 'TrialChoiceID' = -1 for No choice
for ii = 1:numel(DATA)
    if isempty(DATA(ii).TrialChoiceID)
        DATA(ii).TrialChoiceID = -1;
    end
end
clear ii

%% get just the sucesses
success=DATA([DATA.OutcomeID]==0);
%success=DATA([DATA.OutcomeID]==0 | [DATA.OutcomeID]==6);

%% time axis
timeAxis = [DATA.TimeStart];
timeAxis = timeAxis - timeAxis(1);

%% reward tracking summary figure
% x=1:numel(anyChoice);
% 
% 
% figure(9)
% clf
% plot(x,[anyChoice.RecentAvgChoice],'b',x,[anyChoice.UpReward]/100,'c');
% title(['Reward tracking for ' plotDate])
% legend({'Chose top target','Top reward multiplier'})
% ylabel('Probability of choosing target || Reward multiplier value')
% xlabel('Successful trial #')
% ylim([0 1])


%% effort tracking summary figure
n=numel(success);
x=1:n;
Up=zeros(size(x));
Down=zeros(size(x));
choseTop15=zeros(size(x));
Ratio=zeros(size(x));
stdR=success(1).Params.StdReward;

for i=1:n
    Up(i) = (1.0/success(i).UpEffort);
    Down(i) = (1.0/success(i).DownEffort);
    Ratio(i)=log((Down(i)*success(i).UpReward)/(Up(i)*success(i).DownReward));
    
    if i>7 && i<(n-6)
        choseTop15(i) = sum([success(i-7:i+7).TrialChoiceID]==1)/15;
    else
        choseTop15(i)=NaN;
    end
end
clear i n

figure(10)
clf

subplot(4,1,1)
plot(x,[success.UpReward]/stdR,'g',x,[success.DownReward]/stdR,'r');
%legend({'Top reward multiplier', 'Bottom reward multiplier'})
ylabel('% of standard reward')


subplot(4,1,2)
plot(x,Up,'g',x,Down,'r');
%legend({'Top effort multiplier', 'Bottom effort multiplier'})
ylabel('% of standard effort')

subplot(4,1,3)
plot(x, Ratio, 'c');
hline(0,'k--');
m=max(abs(Ratio));
ylim([-1*m m])
clear m
ylabel('log((Ed*Ru)/(Eu*Rd))')

subplot(4,1,4)
plot(x,choseTop15,'b');
hline(0.5,'k--')
ylabel('% chose top target')
xlabel('Successful trial #')

clear Up Down Ratio choseTop15 x

%% plot trial times
% 
% topIdx=[success.TrialChoiceID]==1;
% bottomIdx=[success.TrialChoiceID]==0;
% 
% figure(11)
% clf
% plot(timeAxis_top,[[success(topIdx).TimeEnd]-[success(topIdx).TimeStart]],'bo',timeAxis_bottom,[[success(bottomIdx).TimeEnd]-[success(bottomIdx).TimeStart]],'ro')
% title(['Trial times for ', plotDate ])
% legend({'Chose top target', 'Chose bottom target'})
% ylabel('Trial duration')
% xlabel('Time since experiment start')
% 
% clear topIdx bottomIdx

%% plot %top vs top reward multiplier
% figure(12)
% clf
% plot(topReward,choseTop,'.b')
% title(['Top reward multiplier vs. top choice %, ' plotDate])
% xlabel('Top reward multiplier')
% ylabel('Top choice % (averaged over 10 trials)')

%% Plot rewards
figure(13)
clf
hold on
tmpIdx = [DATA.OutcomeID]==0 & [DATA.TrialChoiceID]==0;
plot(timeAxis(tmpIdx)/3600,[DATA(tmpIdx).ActualReward],'ro')

tmpIdx = [DATA.OutcomeID]==0 & [DATA.TrialChoiceID]==1;
plot(timeAxis(tmpIdx)/3600,[DATA(tmpIdx).ActualReward],'bo')

ylabel('Reward amount received (msec)')
xlabel('Time since experiment start (hours)')

clear tmpIdx

%% Plot outcomes


UpTarget=[DATA.TrialChoiceID]==1;
DownTarget=[DATA.TrialChoiceID]==0;
NoTarget=[DATA.TrialChoiceID]==-1;

figure(14)
clf
hold on
plot(timeAxis(UpTarget)/3600, [DATA(UpTarget).OutcomeID],'bo')
plot(timeAxis(DownTarget)/3600,[DATA(DownTarget).OutcomeID],'ro')
plot(timeAxis(NoTarget)/3600,[DATA(NoTarget).OutcomeID],'ko')
set(gca, 'ytick', [0, 1, 2, 3, 4, 5, 6])
set(gca, 'ytickLabel', {'Success', 'Failed to reach start', 'Failed to hold start', 'Failed to begin reach','Failed to reach target','Wrong way', 'Failed to hold target'})
title(['Trial outcomes for ' plotDate])
xlabel('Time since experiment start (hours)')

clear topIdx bottomIdx UpTarget DownTarget NoTarget

%% Plot total force trace
figure(15)
clf
hold on
subplot(2,1,1)
hold on
for ii = 1:numel(DATA)
    tmpPosY = DATA(ii).ForceTrace(:,5);
    newPosY = tmpPosY*(50);
    tmpTime = DATA(ii).ForceTrace(:,1);
    if DATA(ii).TrialChoiceID==1 && DATA(ii).OutcomeID==0 
        plot(tmpTime/3600-timeAxis(1), newPosY, '.b')
    elseif DATA(ii).TrialChoiceID==0 && DATA(ii).OutcomeID==0
        plot(tmpTime/3600-timeAxis(1), newPosY, '.r')
    elseif DATA(ii).TrialChoiceID==1 && DATA(ii).OutcomeID==6 
        plot(tmpTime/3600-timeAxis(1), newPosY, '.c')
    elseif DATA(ii).TrialChoiceID==0 && DATA(ii).OutcomeID==6
        plot(tmpTime/3600-timeAxis(1), newPosY, '.m')
    elseif DATA(ii).OutcomeID==2 %~(DATA(ii).TrialType==3 || DATA(ii).TrialType==6)
        plot(tmpTime/3600-timeAxis(1), newPosY, '.', 'Color',[0.5 0.5 0.5])
    elseif (DATA(ii).TrialType==6) && (DATA(ii).OutcomeID==5) && (DATA(ii).TrialChoiceID==1)
        plot(tmpTime/3600-timeAxis(1), newPosY, '.', 'Color',[0.75 0.75 1])
    elseif (DATA(ii).TrialType==6) && (DATA(ii).OutcomeID==5) && (DATA(ii).TrialChoiceID==0)
        plot(tmpTime/3600-timeAxis(1), newPosY, '.', 'Color',[1 0.75 0.75])
    else
        plot(tmpTime/3600-timeAxis(1), newPosY, '.k')
    end
end
ylabel('vertical force')
hline(-1)
hline(1,'b:')
hline(3,'b:')
hline(-3)
hline(0,'g--')
ylim([-150,150]/30)



subplot(2,1,2)
hold on
for ii = 1:numel(DATA)
    tmpPosX = DATA(ii).ForceTrace(:,4);
    newPosX = tmpPosX*(50);
    tmpTime = DATA(ii).ForceTrace(:,1);
    if DATA(ii).TrialChoiceID==1 && DATA(ii).OutcomeID==0 
        plot(tmpTime/3600-timeAxis(1), newPosX, '.b')
    elseif DATA(ii).TrialChoiceID==0 && DATA(ii).OutcomeID==0
        plot(tmpTime/3600-timeAxis(1), newPosX, '.r')
    elseif DATA(ii).TrialChoiceID==1 && DATA(ii).OutcomeID==6 
        plot(tmpTime/3600-timeAxis(1), newPosX, '.c')
    elseif DATA(ii).TrialChoiceID==0 && DATA(ii).OutcomeID==6
        plot(tmpTime/3600-timeAxis(1), newPosX, '.m')
    elseif DATA(ii).OutcomeID==2 %~(DATA(ii).TrialType==3 || DATA(ii).TrialType==6)
        plot(tmpTime/3600-timeAxis(1), newPosX, '.', 'Color',[0.5 0.5 0.5])
    elseif (DATA(ii).TrialType==6) && (DATA(ii).OutcomeID==5) && (DATA(ii).TrialChoiceID==1)
        plot(tmpTime/3600-timeAxis(1), newPosX, '.', 'Color',[0.75 0.75 1])
    elseif (DATA(ii).TrialType==6) && (DATA(ii).OutcomeID==5) && (DATA(ii).TrialChoiceID==0)
        plot(tmpTime/3600-timeAxis(1), newPosX, '.', 'Color',[1 0.75 0.75])
    else
        plot(tmpTime/3600-timeAxis(1), newPosX, '.k')
    end
end
ylabel('horizontal force')
ylim([-150,150]/30)
hline(0,'g--')
hline(-70/30)
hline(70/30)
hline(-140/30,'b:')
hline(140/30,'b:')


clear tmpTime newPosX newPosY tmpPosX tmpPosY ii
