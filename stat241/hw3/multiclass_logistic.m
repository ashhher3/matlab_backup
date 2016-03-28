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
%% Newton-Raphson on training data
d=size(x_train,1);
k=size(y_train,1);
n=length(train);

beta = ones(d*(k-1),1); % initialize beta
%beta = [1;1;1;2;2;2;3;3;3];

% scramble the data so we don't get stuck
idx = randperm(n);

% iterate through data
for t=1:n
    x= x_train(:,idx(t));
    
    % while beta has these values, use (1+sum_i^k-1exp(-B_i^Tx) a lot
    denom=1;
    for i=1:(k-1)
        denom = denom + exp(-beta(((i-1)*d+1):i*d)'*x);
    end
    clear i
    
    % evaluate first deriv at current estimate
    grad=zeros(size(beta));
    for a=1:k-1 % for each Beta
        for y=1:k-1 % deriv = sum_y(-x1[y=a] + (xexp(B_a^Tx))/(1+sum_i^k-1exp(-B_i^Tx))
            if y==a % -x1[y=a]
                grad(((a-1)*d+1):a*d)=grad(((a-1)*d+1):a*d) - x;
            end
            % + (xexp(B_a^Tx))/(1+sum_i^k-1exp(-B_i^Tx)
            grad(((a-1)*d+1):a*d)=grad(((a-1)*d+1):a*d) + (x*exp(-1*beta(((a-1)*d+1):a*d)'*x))/(denom);
        end
        clear y
    end
    clear a
    
    % evaluate Hessian at current estimate
    Hess=nan((d*(k-1)),(d*(k-1)));
    for a=1:k-1
        betaA=beta(((a-1)*d+1):a*d);
        for b=1:k-1
            betaB=beta(((b-1)*d+1):b*d);
            if b==a
                Amtx=-1*(x*x')*exp(-1*betaA'*x);
                Hess(((a-1)*d+1):a*d,((a-1)*d+1):a*d)=(k-1)*(denom*Amtx - Amtx*exp(-1*betaA'*x))/(denom^2);
            else
                Amtx=x*exp(-1*betaA'*x);
                Bmtx=x*exp(-1*betaB'*x);
                Hess(((a-1)*d+1):a*d,((b-1)*d+1):b*d)=(k-1)*(Amtx*Bmtx')/(denom^2);
            end
        end
        clear b betaB Amtx Bmtx
    end
    clear a betaA
   
    beta = beta + 0.001*Hess\grad;
    
    clear A B

end


