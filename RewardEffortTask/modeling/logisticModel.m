%% first process data using processData.m

d = DATA_26;% SPECIFIC DATA HERE
[C,R,E,~]=processData(d);
V=[R E ones(1187,1)];
clear d R E

%% run logistic regression

X=logistic(V,C);



