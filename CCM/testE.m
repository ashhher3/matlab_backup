function [Evals, selfPredsX, selfPredsY] = testE(minE, maxE, Eskip, X, Y,tau, L, n)

if n*L > length(X) || n*L >length(Y)
    error('Not enough data given to test n=%d segments of length L=%L\n',n,L)
end

minL = (tau +1)*maxE + 2 - tau;
if L < minL
    error('minimum L for this tau and this maxE is %d, you gave %d\n',minL, L)
end
clear minL

startIdx = 1:L:length(X);
startIdx = startIdx(1:n);

Evals = minE:Eskip:maxE;

XxmapY = NaN(n, length(Evals));
YxmapX = NaN(n, length(Lvals));

for i = 1:length(Evals)
    for j = 1:length(startIdx)
        
        E = Evals(i);
        
        start = startIdx(j);
        
        fprintf('E = %d, iteration %d\n', L, j);
        
        Xtest = X(start:start+L-1);
        Ytest = Y(start:start+L-1);
        
        [ ~, Yest, Xest, ~, ~ ] = CCM( Xtest, Ytest, E, tau );
        
        R = corrcoef(Xtest(size(Xtest,1)-size(Xest,1)+1:end), Xest);
        YxmapX(j,i) = R(1,2);
        
        R = corrcoef(Ytest(size(Xtest,1)-size(Yest,1)+1:end), Yest);
        XxmapY(j,i) = R(1,2);
        
    end
end
    
end