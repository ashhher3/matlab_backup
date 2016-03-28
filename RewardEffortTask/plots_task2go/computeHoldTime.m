%% get just the sucesses
success=DATA([DATA.OutcomeID]==0);
attempt=DATA([DATA.OutcomeID]==0 | [DATA.OutcomeID]==6);

%% figure out how long each trial was in box

% old = cd('/home/motorlab/SUBNETS/Task codes/EffortVSReward');

upPos = [0,60];
upScale = [120,60]/2;
downPos = [0,-60];
downScale = upScale;
cScale = [35 17.5];
sPos = [0, 0];
sScale = [75,40]/2;

for i=1:numel(success) % for each trial
    timeOut=success(i).TimeEnd; %initialize to end of trial if don't leave
    
    %first have to get start target
    for s=1:size(success(i).ForceTrace,1) % for each sample in force trace
        x = 300*success(i).ForceTrace(s,4)*(50/10);
        y = 300*success(i).ForceTrace(s,5)*(50/10);
        
        if y>0
            y = y*success(i).ActualEffort(1)*success(i).Params.UpEScale;
            x = x*0.4;
        elseif y<0
            y = y*success(i).ActualEffort(2);
        end
        
        if TrialInBox([x;y],cScale,sPos,sScale);
            idxS = s ;
            break
        end
    end
    
    %find time when we first enter box
    for j=idxS+1:size(success(i).ForceTrace,1) % for each sample in force trace
        
        x = 300*success(i).ForceTrace(j,4)*(50/10);
        y = 300*success(i).ForceTrace(j,5)*(50/10);
        
        if y>0
            y = y*success(i).ActualEffort(1)*success(i).Params.UpEScale;
            x = x*0.4;
        elseif y<0
            y = y*success(i).ActualEffort(2);
        end
        
        if success(i).TrialChoiceID==1
            inBox = TrialInBox([x;y],cScale,upPos,upScale);
        else
            inBox = TrialInBox([x;y],cScale,downPos,downScale);
        end
        
        if inBox
            idxIn = j ;
            timeIn = success(i).ForceTrace(j,1);
            break
        end
    end
    
    %find time when we first leave box
    for k=idxIn+1:size(success(i).ForceTrace,1) % for each sample in force trace
        
        x = 300*success(i).ForceTrace(k,4)*(50/10);
        y = 300*success(i).ForceTrace(k,5)*(50/10);
        
        if y>0
            y = y*success(i).ActualEffort(1)*success(i).Params.UpEScale;
            x = x*0.4;
        elseif y<0
            y = y*success(i).ActualEffort(2);
        end
        
        if success(i).TrialChoiceID==1
            inBox = TrialInBox([x;y],cScale,upPos,upScale);
        else
            inBox = TrialInBox([x;y],cScale,downPos,downScale);
        end
        
        if ~(inBox)
            idxOut = k ;
            timeOut = success(i).ForceTrace(k,1);
            break
        end
    end
    
    %record duration in box
    success(i).TimeInBox = timeOut - timeIn;
    
end

% do same thing for all attempts
for i=1:numel(attempt) % for each trial
    idxIn=size(attempt(i).ForceTrace,1);
    timeIn=attempt(i).TimeEnd;
    timeOut=attempt(i).TimeEnd; %initialize to end of trial if don't leave
    
    %first have to get start target
    for s=1:size(attempt(i).ForceTrace,1) % for each sample in force trace
        x = 300*attempt(i).ForceTrace(s,4)*(50/10);
        y = 300*attempt(i).ForceTrace(s,5)*(50/10);
        
        if y>0
            y = y*attempt(i).ActualEffort(1)*attempt(i).Params.UpEScale;
            x = x*0.4;
        elseif y<0
            y = y*attempt(i).ActualEffort(2);
        end
        
        if TrialInBox([x;y],cScale,sPos,sScale);
            idxS = s;
            break
        end
    end
    
    %find time when we first enter box
    for j=idxS+1:size(attempt(i).ForceTrace,1) % for each sample in force trace
        
        x = 300*attempt(i).ForceTrace(j,4)*(50/10);
        y = 300*attempt(i).ForceTrace(j,5)*(50/10);
        
        if y>0
            y = y*attempt(i).ActualEffort(1)*attempt(i).Params.UpEScale;
            x = x*0.4;
        elseif y<0
            y = y*attempt(i).ActualEffort(2);
        end
        
        if attempt(i).TrialChoiceID==1
            inBox = TrialInBox([x;y],cScale,upPos,upScale);
        else
            inBox = TrialInBox([x;y],cScale,downPos,downScale);
        end
        
        if inBox
            idxIn = j;
            timeIn = attempt(i).ForceTrace(j,1);
            break
        end
    end
    
    %find time when we first leave box
    for k=idxIn+1:size(attempt(i).ForceTrace,1) % for each sample in force trace
        
        x = 300*attempt(i).ForceTrace(k,4)*(50/10);
        y = 300*attempt(i).ForceTrace(k,5)*(50/10);
        
        if y>0
            y = y*attempt(i).ActualEffort(1)*attempt(i).Params.UpEScale;
            x = x*0.4;
        elseif y<0
            y = y*attempt(i).ActualEffort(2);
        end
        
        if attempt(i).TrialChoiceID==1
            inBox = TrialInBox([x;y],cScale,upPos,upScale);
        else
            inBox = TrialInBox([x;y],cScale,downPos,downScale);
        end
        
        if ~(inBox)
            idxOut = k;
            timeOut = attempt(i).ForceTrace(k,1);
            break
        end
    end
    
    %record duration in box
    attempt(i).TimeInBox = timeOut - timeIn;

    
end

% cd(old);
clear old cScale downPos downScale i idxIn idxOut idxS inDown inUp j k s timeIn timeOut upPos upScale x y yMult xMult sPos sScale inBox
%% save info about dist
holdTimes(1).MinS   = min([success([success.TrialChoiceID]==0).TimeInBox]);
holdTimes(2).MinS   = min([success([success.TrialChoiceID]==1).TimeInBox]);
holdTimes(1).MeanS  = mean([success([success.TrialChoiceID]==0).TimeInBox]);
holdTimes(2).MeanS  = mean([success([success.TrialChoiceID]==1).TimeInBox]);
holdTimes(1).MinA   = min([attempt([attempt.TrialChoiceID]==0).TimeInBox]);
holdTimes(2).MinA   = min([attempt([attempt.TrialChoiceID]==1).TimeInBox]);
holdTimes(1).MeanA  = mean([attempt([attempt.TrialChoiceID]==0).TimeInBox]);
holdTimes(2).MeanA  = mean([attempt([attempt.TrialChoiceID]==1).TimeInBox]);
holdTimes(1).SNum   = length(success([success.TrialChoiceID]==0));
holdTimes(2).SNum   = length(success([success.TrialChoiceID]==1));
holdTimes(1).SRate  = length(success([success.TrialChoiceID]==0))/length(attempt([attempt.TrialChoiceID]==0));
holdTimes(2).SRate  = length(success([success.TrialChoiceID]==1))/length(attempt([attempt.TrialChoiceID]==1));
holdTimes(1).SRAvg  = length(success)/length(attempt);

%%
figure(70)
subplot(1,2,1)
boxplot([success(:).TimeInBox],[success.TrialChoiceID])
xlabel('successful trials')
ylim([0 1.4])
hline(success(1).Params.HoldUp,'c--')
hline(success(1).Params.HoldDown,'m--')


subplot(1,2,2)
boxplot([attempt(:).TimeInBox],[attempt.TrialChoiceID])
xlabel('successful trials + target timeout trials')
ylim([0 1.4])
hline(success(1).Params.HoldUp,'c--')
hline(success(1).Params.HoldDown,'m--')

%%
figure(71)
hold on
plot([attempt([attempt.TrialChoiceID]==0).TimeStart],[attempt([attempt.TrialChoiceID]==0).TimeInBox],'mo')
plot([attempt([attempt.TrialChoiceID]==1).TimeStart],[attempt([attempt.TrialChoiceID]==1).TimeInBox],'co')
plot([success([success.TrialChoiceID]==0).TimeStart],[success([success.TrialChoiceID]==0).TimeInBox],'ro')
plot([success([success.TrialChoiceID]==1).TimeStart],[success([success.TrialChoiceID]==1).TimeInBox],'bo')
hline(success(1).Params.HoldUp,'c--')
hline(success(1).Params.HoldDown,'m--')

