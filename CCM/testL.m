function [Lvals, XxmapX ] = testE(minL, Lskip, maxL, X,tau, E, n)

% % input checking
% if n*E > length(X)
%     error('Not enough data given to test n=%d segments of length L=%d\n',n,E)
% end
% if tau < maxL
%     error('Best results when tau>= maxE\n')
% end
% minL = (tau +1)*maxL + 2 - tau;
% if E < minL
%     error('minimum L for this tau and this maxE is %d, you gave %d\n',minL, E)
% end
% clear minL

% initialize variables
startIdx = 1:maxL:length(X);
startIdx = startIdx(1:n);
Lvals = minL:Lskip:maxL;
XxmapX = NaN(n, length(Lvals));

% do CCM
for i = 1:length(Lvals)
    L = Lvals(i);
    for j = 1:length(startIdx)
        start = startIdx(j);
        fprintf('L = %d, iteration %d\n', L, j);
        
        Xtest = X(start:start+L-1);
        Mx = makeShadowManifold(Xtest, E, tau);
        [Xtrue, Xest] = testShadowManifold( Mx, E);
        
        R = corrcoef(Xtrue, Xest);
        XxmapX(j,i) = R(1,2);
    
    end
end

% make figure
fig = figure();
set (fig, 'Units', 'normalized', 'Position', [0,0,0.9,0.9], 'Color', [1,1,1]);
hold on
x=Lvals;
y1=XxmapX;
h1=shadedErrorBar(x,median(y1),[abs(prctile(y1,75)-median(y1)); abs(prctile(y1,25)-median(y1))],'g',1);
xlim([minL maxL]) 
box on
title({'Accuracy of self-predictions using time-lagged manifold',sprintf('(E=%d, \\tau = %d)',E, tau)})
xlabel('L')
ylabel('\rho')
legend([h1.mainLine],'X xmap X', 'Location', 'best')
hold off


end