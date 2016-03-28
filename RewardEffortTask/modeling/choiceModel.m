function [params, pred, fig] = choiceModel(train, test)
% load data

% process data and save as global variables
global Ctrue
global Rtrue
global Etrue
global num
global idx

[Ctrue,Rtrue,Etrue,num,idx] = processData(train);
[Ctest,Rtest,Etest,n,i] = processData(test);

% train params
params = optimizeParams();
paramsInit = [100, 1, 1, 0.5, 1, 4, 0, 1, 1.003];
err1=predictionError(params)
err2=predictionError(paramsInit)

% test params and plot
pred = predictChoice(Rtest,Etest,n,params);
fig = RETplot(Ctest,Rtest,Etest,n,i,pred);

end

function [C,R,E,n,i,success] = processData(DATA)
success=DATA([DATA.OutcomeID]==0);
n=numel(success);
i=1:n;

% true choice, reward, effort (rows are trials)
C=[success.TrialChoiceID]'; % 1 for U, 0 for D
R=[success.UpReward; success.DownReward]'; % clm 1 is U, clm 2 is D
E=[success.UpEffort; success.DownEffort]'; % clm 1 is U, clm 2 is D

end

function params = optimizeParams()
% params:
%  1 - Rod           initial est of D reward
%  2 - Eou           initial est of U effort
%  3 - Eod           initial est of D effort
%  4 - alphaR        learning rate for reward
%  5 - alphaE        learning rate for effort
%  6 - beta          slope for sigmoid
%  7 - gamma         bias for sigmoid
%  8 - delta         effort scale
%  9 - epsilR        exploration rate for reward

% paramsInit = [100, 1, 1, 0.5, 1, 4, 0, 1, 1.003]; % order is listed order above
lb = [0,     0,   0, 0, 0,      0, -10,  0, 0];
ub = [200,   2,   2, 1, 1, 100000,  10, 2, 2];
% a=rand();
% paramsInit=a*lb + (1-a)*ub
paramsInit = [100, 1, 1, 0.5, 1, 4, 0, 1, 1.003];

params = fmincon(@predictionError, paramsInit,[],[],[],[],lb,ub)

end

function err = predictionError(params)
% params:
%  1 - Rod           initial est of D reward
%  2 - Eou           initial est of U effort
%  3 - Eod           initial est of D effort
%  4 - alphaR        learning rate for reward
%  5 - alphaE        learning rate for effort
%  6 - beta          slope for sigmoid
%  7 - gamma         bias for sigmoid
%  8 - delta         effort scale
%  9 - epsilR        exploration rate for reward

global num
global Ctrue
global Rtrue
global Etrue

Rest        = NaN(num,2); % clm 1 is U, clm 2 is D
Eest        = NaN(num,2); % clm 1 is U, clm 2 is D
Rest(1,1)   = sum(Rtrue(1,:))/2; % fixed param
Rest(1,2)   = params(1);
Eest(1,1)   = params(2);
Eest(1,2)   = params(3);
alphaR      = params(4);
alphaE      = params(5);
beta        = params(6);
gamma       = params(7);
delta       = params(8);
epsilR      = params(9);

err = 0; %initialize error to 0

for i=1:num
    % predict the prob of choosing U on i-th choice
    val=log(delta*(Rest(i,1)*Eest(i,2))/(Rest(i,2)*Eest(i,1)));
    predUpProb=sigmoid(val,beta,gamma);
    
    err = err + abs(Ctrue(i) - predUpProb); % calculate error
    
    % update reward estimate - only get info about reward you choose
    uu=(1-alphaR)*Rest(i,1)+alphaR*Rtrue(i,1);
    ud=epsilR*Rest(i,2);
    du=epsilR*Rest(i,1);
    dd=(1-alphaR)*Rest(i,2)+alphaR*Rtrue(i,2);
    
    Rest(i+1,1)=predUpProb*uu + (1-predUpProb)*du;
    Rest(i+1,2)=predUpProb*ud + (1-predUpProb)*dd;
    
    % update effort estimate - can get info about pushing both directions
    Eest(i+1,1)=(1-alphaE)*Eest(i,1)+alphaE*Etrue(i,1);
    Eest(i+1,2)=(1-alphaE)*Eest(i,2)+alphaE*Etrue(i,2);
    
end

end

function Cpred = predictChoice(Rtest,Etest,num,params)
% params:
%  1 - Rod           initial est of D reward
%  2 - Eou           initial est of U effort
%  3 - Eod           initial est of D effort
%  4 - alphaR        learning rate for reward
%  5 - alphaE        learning rate for effort
%  6 - beta          slope for sigmoid
%  7 - gamma         bias for sigmoid
%  8 - delta         effort scale
%  9 - epsilR        exploration rate for reward

Cpred       = NaN(size(Etest)); % clm 1 is prob, clm 2 is choice
Rest        = NaN(num,2); % clm 1 is U, clm 2 is D
Eest        = NaN(num,2); % clm 1 is U, clm 2 is D
Rest(1,1)   = sum(Rtest(1,:))/2; % fixed param
Rest(1,2)   = params(1);
Eest(1,1)   = params(2);
Eest(1,2)   = params(3);
alphaR      = params(4);
alphaE      = params(5);
beta        = params(6);
gamma       = params(7);
delta       = params(8);
epsilR      = params(9);

for i=1:num
    % predict the prob of choosing U on i-th choice
    val=log(delta*(Rest(i,1)*Eest(i,2))/(Rest(i,2)*Eest(i,1)));
    Cpred(i,1)=sigmoid(val,beta,gamma);
    
    % simulate choice & update reward estimates
    if rand()<=Cpred(i,1)
        Cpred(i,2)=1;
        Rest(i+1,1)=(1-alphaR)*Rest(i,1)+alphaR*Rtest(i,1); 
        Rest(i+1,2)=epsilR*Rest(i,2);
    else
        Cpred(i,2)=0;
        Rest(i+1,1)=epsilR*Rest(i,1);
        Rest(i+1,2)=(1-alphaR)*Rest(i,2)+alphaR*Rtest(i,2);
    end
    
    % update effort estimates
    Eest(i+1,1)=(1-alphaE)*Eest(i,1)+alphaE*Etest(i,1); % update Eu
    Eest(i+1,2)=(1-alphaE)*Eest(i,2)+alphaE*Etest(i,2); % update Ed
end


end

function y = sigmoid(x,b,c)

y = 1/(1+exp(-1*b*(x-c)));

end

function fig = RETplot(Ctest,Rtest,Etest,n,idx,Cpred)
testAvg=nan(size(idx));
predAvg=nan(size(idx));
Rtot=Rtest(1,1)+Rtest(1,2);

for i=1:n
    % invert effort values to get percentage
    Etest(i,1) = (1.0/Etest(i,1));
    Etest(i,2) = (1.0/Etest(i,2));
    
    if i>7 && i<(n-6) % compute average choices
        testAvg(i) = sum(Ctest(i-7:i+7))/15;
        predAvg(i) = sum(Cpred(i-7:i+7,2))/15;
    end
end

fig=figure(20);
clf
hold on
subplot(3,1,1)
plot(idx,Rtest(:,1)/Rtot,'g',idx,Rtest(:,2)/Rtot,'r');
legend({'Up', 'Down'})
ylabel('% of standard reward')
xlabel('Successful trial #')

subplot(3,1,2)
plot(idx,Etest(:,1),'g',idx,Etest(:,2),'r');
legend({'Up', 'Down'})
ylabel('% of standard effort')

subplot(3,1,3)
hold on
plot(idx,testAvg,'b');
plot(idx, predAvg, 'c');
legend({'True','Predicted'})
ylabel('% chose top target')

end



