function [states, actions] = DataToStateSpace(data)

success = data([data.OutcomeID]==0);

n = length(success);

states=NaN(2, n+1);
actions=NaN(1, n+1);

states(:,1)=[0;0]; % initial state is no reward, no effort

stdR=success(1).Params.StdReward;

for i=1:n % convert each trial into action i and resulting state i+1
    
    if i==1
        actions(i)=1;
    else
        actions(i)=abs(success(i).TrialChoiceID - success(i-1).TrialChoiceID);
    end
    % now 0 means same target, 1 means switched targets
    
    states(1,i+1)=success(i).ActualReward/stdR;
    if success(i).TrialChoiceID==0
        idx=2; %down effort is at pos 2
    else
        idx=1; %  up effort is at pos 1
    end
    states(2,i+1)=1/success(i).ActualEffort(idx);
end


end