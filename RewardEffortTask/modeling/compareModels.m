%% - name test and training the right thing
train = DATA_27;
test = DATA_26;

%% process data using processData.m
[C,R,E,~]=processData(train);
V=[R E ones(length(R),1)];
clear R E

[Ctest,Rtest,Etest,~]=processData(test);
Vtest=[Rtest Etest ones(length(Rtest),1)];
clear Rtest Etest

%% run logistic regression
X=logistic(V,C);

%% run learning model
[alpha, beta, Xhat, iteration]=learningModel(C,V,X);

%% or, saved coefficients
X = [2.23513121355013;
     -0.440466670990652;
     -4.85532121273346;
     2.44747478811242;
     0.897332925983431];

Xhat = [2.54958818229972;
        -0.387122390851063;
        -8.05300269354427;
        5.13995673147509;
        1.08123295548491];

alpha = 0.518923846567515;

beta = 0.138852010774656;

%% make predictions and stuff
Vhat = estimateValue(alpha,beta,V);
Vhat_test = estimateValue(alpha,beta,Vtest);

lng_pred = 1./(1+exp(-Vhat*Xhat));
lng_pred_test = 1./(1+exp(-Vhat_test*Xhat));

log_pred = 1./(1+exp(-V*X));
log_pred_test = 1./(1+exp(-Vtest*X));

%% logistic model errors
err(1,1)=evaluatePrediction(log_pred,C);
err(2,1)=evaluatePrediction(log_pred_test,Ctest);

%% learning model errors
err(1,2)=evaluatePrediction(lng_pred,C);
err(2,2)=evaluatePrediction(lng_pred_test,Ctest);

%% average errors
err(1,3) = 0.193577720943450; % 10/26
err(2,3) = 0.262381852551984; % 10/27

%% make plots for training data
success=train([train.OutcomeID]==0);
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
    Ratio(i)=log(success(i).UpReward/success(i).DownReward)-log(Up(i)/Down(i));
    
    if i>7 && i<(n-6)
        choseTop15(i) = sum([success(i-7:i+7).TrialChoiceID]==1)/15;
    else
        choseTop15(i)=NaN;
    end
end
clear i n

figure(10)
clf

subplot(3,1,1)
plot(x,[success.UpReward]/stdR,'g',x,[success.DownReward]/stdR,'r');
%legend({'Top reward multiplier', 'Bottom reward multiplier'})
ylabel('% of standard reward')


subplot(3,1,2)
hold on
plot(x,Up,'g',x,Down,'r');
hline(1,'k')
%legend({'Top effort multiplier', 'Bottom effort multiplier'})
ylabel('% of standard effort')

subplot(3,1,3)
hold on
plot(x,choseTop15,'k');
plot(x,log_pred,'c');
plot(x,lng_pred,'m');
ylabel('% chose top target')
xlabel('Successful trial #')

clear Up Down Ratio choseTop15 x success stdR

%% make plots for test data
success=test([test.OutcomeID]==0);
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
    Ratio(i)=log(success(i).UpReward/success(i).DownReward)-log(Up(i)/Down(i));
    
    if i>7 && i<(n-6)
        choseTop15(i) = sum([success(i-7:i+7).TrialChoiceID]==1)/15;
    else
        choseTop15(i)=NaN;
    end
end
clear i n

figure(11)
clf

subplot(3,1,1)
plot(x,[success.UpReward]/stdR,'g',x,[success.DownReward]/stdR,'r');
%legend({'Top reward multiplier', 'Bottom reward multiplier'})
ylabel('% of standard reward')


subplot(3,1,2)
hold on
plot(x,Up,'g',x,Down,'r');
hline(1,'k')
%legend({'Top effort multiplier', 'Bottom effort multiplier'})
ylabel('% of standard effort')

subplot(3,1,3)
hold on
plot(x,choseTop15,'k');
plot(x,log_pred_test,'c');
plot(x,lng_pred_test,'m');
ylabel('% chose top target')
xlabel('Successful trial #')

clear Up Down Ratio choseTop15 x Rtest Etest success stdR

%% clear extra variables leaving just data
clear C Ctest lng_pred lng_pred_test log_pred log_pred_test test train V Vhat Vhat_test Vtest err alpha beta X Xhat X_2 ans