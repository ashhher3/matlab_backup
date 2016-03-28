%% import data and store it in variables of the right shape
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
y_train=nan(4,length(train));
y_test=y_train;
for i=1:length(train)
    switch train(i,4)
        case 1
            y_train(:,i) = [1;0;0;0];
        case 2
            y_train(:,i) = [0;1;0;0];
        case 3
            y_train(:,i) = [0;0;1;0];
        case 4
            y_train(:,i) = [0;0;0;1];
    end
    
    switch test(i,4)
        case 1
            y_test(:,i) = [1;0;0;0];
        case 2
            y_test(:,i) = [0;1;0;0];
        case 3
            y_test(:,i) = [0;0;1;0];
        case 4
            y_test(:,i) = [0;0;0;1];
    end
end

clear i

%% solve MLE for training data
beta = (x_train*x_train')\(x_train*y_train');

%% classify training data
p_train=nan(4,length(train));
labels_train=nan(length(train),1);
for i=1:length(train)
    p_train(1,i) = (2*pi)^(-2)*exp( -0.5 * ([1;0;0;0]- beta'*x_train(:,i))'*([1;0;0;0]- beta'*x_train(:,i)) );
    p_train(2,i) = (2*pi)^(-2)*exp( -0.5 * ([0;1;0;0]- beta'*x_train(:,i))'*([0;1;0;0]- beta'*x_train(:,i)) );
    p_train(3,i) = (2*pi)^(-2)*exp( -0.5 * ([0;0;1;0]- beta'*x_train(:,i))'*([0;0;1;0]- beta'*x_train(:,i)) );
    p_train(4,i) = (2*pi)^(-2)*exp( -0.5 * ([0;0;0;1]- beta'*x_train(:,i))'*([0;0;0;1]- beta'*x_train(:,i)) );
    [~,labels_train(i)] = max(p_train(:,i));
end

correct_train = [labels_train == train(:,4)];
missclass_train = 1 - sum(correct_train)/length(correct_train);

clear i

%% classify testing data
p_test=nan(4,length(test));
labels_test=nan(length(test),1);
for i=1:length(test)
    p_test(1,i) = (2*pi)^(-2)*exp( -0.5 * ([1;0;0;0]- beta'*x_test(:,i))'*([1;0;0;0]- beta'*x_test(:,i)) );
    p_test(2,i) = (2*pi)^(-2)*exp( -0.5 * ([0;1;0;0]- beta'*x_test(:,i))'*([0;1;0;0]- beta'*x_test(:,i)) );
    p_test(3,i) = (2*pi)^(-2)*exp( -0.5 * ([0;0;1;0]- beta'*x_test(:,i))'*([0;0;1;0]- beta'*x_test(:,i)) );
    p_test(4,i) = (2*pi)^(-2)*exp( -0.5 * ([0;0;0;1]- beta'*x_test(:,i))'*([0;0;0;1]- beta'*x_test(:,i)) );
    [~,labels_test(i)] = max(p_test(:,i));
end

correct_test = [labels_test == test(:,4)];
missclass_test = 1 - sum(correct_test)/length(correct_test);

clear i

%% make plot

% separate data by true label
idx1=find(train(:,4)==1);
idx2=find(train(:,4)==2);
idx3=find(train(:,4)==3);
idx4=find(train(:,4)==4);

% make more predictions
x1=-7:0.1:2;
x2=-7:0.1:2;
for i=1:length(x1)
    for j=1:length(x2)
        prob(1) = (2*pi)^(-2)*exp( -0.5 * ([1;0;0;0]- beta'*[1;x1(i);x2(j)])'*([1;0;0;0]- beta'*[1;x1(i);x2(j)]) );
        prob(2) = (2*pi)^(-2)*exp( -0.5 * ([0;1;0;0]- beta'*[1;x1(i);x2(j)])'*([0;1;0;0]- beta'*[1;x1(i);x2(j)]) );
        prob(3) = (2*pi)^(-2)*exp( -0.5 * ([0;0;1;0]- beta'*[1;x1(i);x2(j)])'*([0;0;1;0]- beta'*[1;x1(i);x2(j)]) );
        prob(4) = (2*pi)^(-2)*exp( -0.5 * ([0;0;0;1]- beta'*[1;x1(i);x2(j)])'*([0;0;0;1]- beta'*[1;x1(i);x2(j)]) );
        [~,labels(i,j)]=max(prob);
    end
end

figure(10)
clf
hold on
plot(train(idx1,1), train(idx1,2), 'co');
plot(train(idx2,1), train(idx2,2), 'bo');
plot(train(idx3,1), train(idx3,2), 'go');
plot(train(idx4,1), train(idx4,2), 'ro');
contour(x1,x2,labels,[1 2 3 4]);
xlabel('v1')
ylabel('v2')
legend('y=1','y=2','y=3','y=4','boundaries')








