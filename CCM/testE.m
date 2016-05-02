function [Evals, XxmapX, YxmapY, XxmapY] = testE(minE, Eskip, maxE, X, Y,tau, L, n)

% input checking
if n*L > length(X) || n*L >length(Y)
    error('Not enough data given to test n=%d segments of length L=%d\n',n,L)
end

if tau < maxE
    error('Best results when tau>= maxE\n')
end

minL = (tau +1)*maxE + 2 - tau;
if L < 2*minL
    error('minimum L for this tau and this maxE is %d, you gave %d\n',minL, L)
end
clear minL

% initialize variables
startIdx = 1:L:length(X);
startIdx = startIdx(1:n);

Evals = minE:Eskip:maxE;

XxmapX = NaN(n, length(Evals));
XxmapY = NaN(n, length(Evals));
YxmapY = NaN(n, length(Evals));

% do CCM
for i = 1:length(Evals)
    E = Evals(i);
    
    for j = 1:length(startIdx)
        
        start = startIdx(j);
        
        fprintf('E = %d, iteration %d\n', E, j);
        
        Xtest = X(start:start+L-1);
        Ytest = Y(start:start+L-1);
        
        [ ~, Yest, Xself, Yself, ~, ~ ] = CCM( Xtest, Ytest, E, tau );
        
        R = corrcoef(Xtest(size(Xtest,1)-size(Xself,1)+1:end), Xself);
        XxmapX(j,i) = R(1,2);
        
        R = corrcoef(Ytest(size(Xtest,1)-size(Xself,1)+1:end), Yest);
        XxmapY(j,i) = R(1,2);
        
        R = corrcoef(Ytest(size(Xtest,1)-size(Yself,1)+1:end), Yself);
        YxmapY(j,i) = R(1,2);
        
    end
end

% % make figure
% fig = figure();
% set (fig, 'Units', 'normalized', 'Position', [0,0,0.9,0.9], 'Color', [1,1,1]);
% hold on
% x=Evals;
% y1=XxmapX;
% y2=YxmapY;
% h1=shadedErrorBar(x,median(y1),[abs(prctile(y1,75)-median(y1)); abs(prctile(y1,25)-median(y1))],'g',1);
% h2=shadedErrorBar(x,median(y2),[abs(prctile(y2,75)-median(y2)); abs(prctile(y2,25)-median(y2))],'b',1);
% xlim([minE maxE]) 
% box on
% title({'Accuracy of self-predictions using time-lagged manifold',sprintf('(L=%d, \\tau = %d)',L, tau)})
% xlabel('E (embedding dimension)')
% ylabel('\rho')
% legend([h1.mainLine, h2.mainLine],'X xmap X','Y xmap Y', 'Location', 'best')
% hold off


end