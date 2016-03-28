function [alpha, beta]=fitLearningRates(C,V,X)
% C, R, E must come from processData.m
% X is logistic regression params from logistic_pred.m

global c
global v
global x

c=C;
v=V;
x=X;

params = fmincon(@predictionError, [0.5, 0.5],[],[],[],[],[0,0],[1,1]);

alpha = params(1);
beta=params(2);
end

function err = predictionError(params)
global c
global v
global x

alpha = params(1);
beta = params(2);

Vhat=estimateValue(alpha,beta,v);

pred = 1./(1+exp(-Vhat*x));

err = evaluatePrediction(pred,c);
end