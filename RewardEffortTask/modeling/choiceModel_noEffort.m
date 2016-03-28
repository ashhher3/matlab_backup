function [params, pred, fig] = choiceModel_noEffort(train, test)

% process data and save as global variables
global Ctrue
global Rtrue
global num
global idx
[Ctrue,Rtrue,num,idx] = processData(train);
[Ctest,Rtest,n,i] = processData(test);

% train params
params = optimizeParams();
paramsInit = [125,      0.5,         4,      0,      1.003];
err1=predictionError(params)
err2=predictionError(paramsInit)

% test params and plot
pred = predictChoice(Rtest,n,params);
fig = RETplot(Ctest,Rtest,n,i,pred);

end

function [C,R,n,i] = processData(DATA)
success=DATA([DATA.OutcomeID]==0);
n=numel(success);
i=1:n;

% true choice, reward, effort (rows are trials)
C=[success.TrialChoiceID]'; % 1 for U, 0 for D
R=[success.UpReward; success.DownReward]'; % clm 1 is U, clm 2 is D

end

function params = optimizeParams()
% params:
%  1 - Rod           initial est of D reward
%  2 - alphaR        learning rate for reward
%  3 - beta          slope for sigmoid
%  4 - gamma         bias for sigmoid
%  5 - epsilR        exploration rate for reward

% paramsInit = [125,      0.5,         4,      0,      1.003]; % order is listed order above
lb         = [0,          0,         0, -10,         0];
ub         = [200,        1,    100000,  10,    2];
a=rand();
paramsInit=a*lb + (1-a)*ub
% paramsInit = [125,      0.5,         4,      0,      1.003];

params = fmincon(@predictionError, paramsInit,[],[],[],[],lb,ub)

end

function err = predictionError(params)
% params:
%  1 - Rod           initial est of D reward
%  2 - alphaR        learning rate for reward
%  3 - beta          slope for sigmoid
%  4 - gamma         bias for sigmoid
%  5 - epsilR        exploration rate for reward

global num
global Ctrue
global Rtrue

Rest        = NaN(num,2); % clm 1 is U, clm 2 is D
Rest(1,1)   = sum(Rtrue(1,:))/2; % fixed param
Rest(1,2)   = params(1);
alphaR      = params(2);
beta        = params(3);
gamma       = params(4);
epsilR      = params(5);

err = 0; %initialize error to 0

for i=1:num
    % predict the prob of choosing U on i-th choice
    val=log((Rest(i,1))/(Rest(i,2)));
    predUpProb=sigmoid(val,beta,gamma);
    
    err = err + abs(Ctrue(i) - predUpProb); % calculate error
    
    % update reward
    uu=(1-alphaR)*Rest(i,1)+alphaR*Rtrue(i,1);
    ud=epsilR*Rest(i,2);
    du=epsilR*Rest(i,1);
    dd=(1-alphaR)*Rest(i,2)+alphaR*Rtrue(i,2);
    
    Rest(i+1,1)=predUpProb*uu + (1-predUpProb)*du;
    Rest(i+1,2)=predUpProb*ud + (1-predUpProb)*dd;
    
end
end

function Cpred = predictChoice(Rtest,num,params)
% params:
%  1 - Rod           initial est of D reward
%  2 - alphaR        learning rate for reward
%  3 - beta          slope for sigmoid
%  4 - gamma         bias for sigmoid
%  5 - epsilR        exploration rate for reward

Cpred       = NaN(size(Rtest)); % clm 1 is prob, clm 2 is choice
Rest        = NaN(num,2); % clm 1 is U, clm 2 is D
Rest(1,1)   = sum(Rtest(1,:))/2; % fixed param
Rest(1,2)   = params(1);
alphaR      = params(2);
beta        = params(3);
gamma       = params(4);
epsilR      = params(5);

for i=1:num
    % predict the prob of choosing U on i-th choice
    val=log((Rest(i,1))/(Rest(i,2)));
    Cpred(i,1)=sigmoid(val,beta,gamma);
    
    if rand()<=Cpred(i,1) % if we choose up
        Cpred(i,2)=1; % record choice
        Rest(i+1,1)=(1-alphaR)*Rest(i,1)+alphaR*Rtest(i,1); % update Ru
        Rest(i+1,2)=epsilR*Rest(i,2); % update Rd       
    else % if we chose down
        Cpred(i,2)=0; % record choice
        Rest(i+1,1)=epsilR*Rest(i,1); % update Ru
        Rest(i+1,2)=(1-alphaR)*Rest(i,2)+alphaR*Rtest(i,2); % update Rd     
    end
end
end

function y = sigmoid(x,b,c)

y = 1/(1+exp(-1*b*(x-c)));

end

function fig = RETplot(Ctest,Rtest,n,idx,Cpred)
testAvg=nan(size(idx));
predAvg=nan(size(idx));
Rtot=Rtest(1,1)+Rtest(1,2);

for i=1:n    
    if i>7 && i<(n-6)
        testAvg(i) = sum(Ctest(i-7:i+7))/15;
        predAvg(i) = sum(Cpred(i-7:i+7,2))/15;
    end
end

fig=figure();
clf
hold on
subplot(2,1,1)
plot(idx,Rtest(:,1)/Rtot,'g',idx,Rtest(:,2)/Rtot,'r');
legend({'Up', 'Down'})
ylabel('% of standard reward')
xlabel('Successful trial #')

subplot(2,1,2)
hold on
plot(idx,testAvg,'b');
plot(idx, predAvg, 'c');
legend({'True','Predicted'})
ylabel('% chose top target')

end



