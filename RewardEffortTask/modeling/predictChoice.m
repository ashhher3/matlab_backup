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