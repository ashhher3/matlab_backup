function [C,R,E,i] = processData(DATA)
success=DATA([DATA.OutcomeID]==0);
n=numel(success);
i=1:n;

stdR=success(1).Params.StdReward;

% true choice, reward, effort (rows are trials)
C=[success.TrialChoiceID]'; % 1 for U, 0 for D
R=[success.UpReward; success.DownReward]'; % clm 1 is U, clm 2 is D
R=R/stdR;
E=1./[success.UpEffort; success.DownEffort]';

end