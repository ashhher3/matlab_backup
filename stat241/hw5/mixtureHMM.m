%% load data
clear all
train = importdata('hw5-2.data');
test = importdata('hw5-2.test');

%% initialize parameters

pi = [1/3;1/3;1/3]; % it says not to estimate initial state but it has to be something
A = ones(3,3)/3;
R = ones(3,2)/2;
lambda = [ 1, 5; 50, 100; 200, 300];

%% run EM
m=3; %if written as function, would pass in these, data, and initial params 
k=2;
T=70;
y=(train);

% make structures to hold alphas and betas
alpha = NaN(m,k,T);
aNormLog = NaN(T,1);
beta = NaN(m,k,T);

for iteration = 1:100
    % I can only do 6 iterations before things become NaN, and I don't know
    % whether that's because I'm updating wrong or need to do more tricks
    % to store things right?? I tried to follow every hint on piazza
    
    % 1 - calculate alphas
        % initialize apha(q_1,x_1) = p(y_1|x_1,q_1)r_q1x1 pi_q1
        for q1=1:m
            for x1=1:k
                alpha(q1,x1,1)=poisspdf(y(1),lambda(q1,x1))*R(q1,x1)*pi(q1);
                
            end
        end
        aNorm=sum(sum(alpha(:,:,1)));
        alpha(:,:,1) = alpha(:,:,1)/aNorm;
        aNormLog(1)=log(aNorm);
        clear q1 x1 aNorm
        
        % loop over rest of data points
        for tt=2:T
            for qtt=1:m
                for xtt=1:k
                    probNextY=poisspdf(y(tt),lambda(qtt,xtt));
                    % needs to be multiplied by a sum over previous qt,xt
                    s = 0;
                    for qt=1:m
                        for xt=1:k
                            %p=alpha(qt,xt,tt-1)*R(qtt,xtt)*A(qt,qtt);
                            p=alpha(qt,xt,tt-1)*exp(aNormLog(tt-1))*R(qtt,xtt)*A(qt,qtt);
                            s=s+p;
                        end
                    end
                    alpha(qtt,xtt,tt)=s*probNextY;
                end
            end
            aNorm=sum(sum(alpha(:,:,tt)));
            alpha(:,:,tt) = alpha(:,:,tt)/aNorm;
            aNormLog(tt)=log(aNorm);
        end
        clear tt qtt xtt probNextY s qt xt p aNorm
    
    % 2 - calculate betas
        % initialize beta(T) to be all ones
        for qT=1:m
            for xT=1:k
                beta(qT,xT,T)=1;
            end
        end
        clear qT xT
        
        % loop over rest of data points
        for t=T-1:-1:1
            for qt=1:m
                for xt=1:k
                    s=0;
                    for qtt=1:m
                        for xtt=1:k
                            probY=poisspdf(y(t+1),lambda(qtt,xtt));
                            p=beta(qtt,xtt,t+1)*probY*R(qtt,xtt)*A(qt,qtt);
                            s=s+p;
                        end
                    end
                    beta(qt,xt,t)=s;
                end
            end
        end
        clear t qt xt s qtt xtt probY p
    
    % 3 - update parameters (question says not to estimate initial state)
    
        % compute p(y) since it comes up a lot
        py = 0;
        t=1; % should be able to choose any t
        for qt=1:m
            for xt=1:k
                %py=py+alpha(qt,xt,t)*beta(qt,xt,t);
                py=py+alpha(qt,xt,t)*exp(aNormLog(t))*beta(qt,xt,t);
            end
        end
        clear t qt xt
        
        % new A matrix
        Anew=NaN(size(A));
        for qt=1:m
            for qtt=1:m
                num=0;
                for t=1:T-1
                    term1=0;
                    for xt=1:k
                        %term1=term1+alpha(qt,xt,t);
                        term1=term1+alpha(qt,xt,t)*exp(aNormLog(t));
                    end
                    
                    term2=0;
                    for xtt=1:k
                        prod=poisspdf(y(t+1),lambda(qtt,xtt))*R(qtt,xtt)^2*beta(qtt,xtt,t+1);
                        term2=term2+prod;
                    end
        
                    num=num+(term1*term2*A(qtt,qt)/py);
                end
                Anew(qt,qtt)=num;
            end
        end
        clear qt qtt t term1 xt term2 xtt prod num
        
        % new R matrix
        Rnew=NaN(size(R));
        for qt=1:m
            for xt=1:k
                num=0;
                for t=1:T
                    %num=num+(alpha(qt,xt,t)*beta(qt,xt,t)/py);
                    num=num+(alpha(qt,xt,t)*exp(aNormLog(t))*beta(qt,xt,t)/py);
                end
                Rnew(qt,xt)=num;
            end
        end
        clear qt xt num t
        
        % new lambda
        lambdaN=NaN(size(lambda));
        for qt=1:m
            for xt=1:k
                num=0;
                denom=0;
                for t=1:T
                    %num=num+(alpha(qt,xt,t)*beta(qt,xt,t)*y(t)/py);
                    num=num+(alpha(qt,xt,t)*exp(aNormLog(t))*beta(qt,xt,t)*poisspdf(y(t),lambda(qt,xt))/py);
                    denom=denom+(alpha(qt,xt,t)*exp(aNormLog(t))*beta(qt,xt,t));
                end
                lambdaN(qt,xt)=num/denom;
            end
        end
        clear qt xt num t
        
        % normalize and actually update
        for q=1:m
            A(q,:) = Anew(q,:) / sum(Anew(q,:));
            R(q,:) = Rnew(q,:) / sum(Rnew(q,:));
        end
        lambda = lambdaN;
        clear Anew Rnew lambdaN
        
        if isnan(sum(sum(A))) || isnan(sum(sum(R))) || isnan(sum(sum(lambda)))
            iteration
            break
        end
    
end



