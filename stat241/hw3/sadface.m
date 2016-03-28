%% import data and store it in variables of the right shape
clear all

train=importdata('hw3-2-train.data');
test=importdata('hw3-2-test.data');

% from question: x_i = (1, v_i1, v_i2)^T
x_train=nan(3,length(train));
x_test=x_train; % testing & training data are same size
for i=1:length(train)
    x_train(1,i)=train(i,3);
    x_train(2,i)=train(i,1);
    x_train(3,i)=train(i,2);
    x_test(1,i)=test(i,3);
    x_test(2,i)=test(i,1);
    x_test(3,i)=test(i,3);
end

% from question: y~_i = e_yi
y_train=train(:,4);
y_test=test(:,4);

%%
beta = logistic(x_train', y_train);

%% classify training data
p_train=nan(4,length(train));
labels_train=nan(length(train),1);

% E(Y) = 1 ./ (1+exp(-A*X))
for i=1:length(train)
   labels(i,:)=(1 ./ (1+ exp(-1*x_train(:,i)'*beta)));
end

correct_train = [labels_train == train(:,4)];
missclass_train = 1 - sum(correct_train)/length(correct_train);

clear i