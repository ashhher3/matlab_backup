function [Evals, XxmapX, XxmapY ] = testE(minE, Eskip, maxE, X, Y, tau, L, n)

% input checking
extra = (maxE-1)*tau;
if n*(L+extra) > length(X)
    error('Not enough data given to test n=%d segments of length L=%d\n',n,L)
end
% if tau < maxE
%     error('Best results when tau>= maxE\n')
% end
minL = (tau +1)*maxE + 2 - tau;
if L < minL
    error('minimum L for this tau and this maxE is %d, you gave %d\n',minL, L)
end
clear minL

% initialize variables
startIdx = 1:L+extra:length(X);
startIdx = startIdx(1:n);
Evals = minE:Eskip:maxE;
XxmapX = NaN(n, length(Evals));
XxmapY = NaN(n, length(Evals));

% do CCM
for i = 1:length(Evals)
    E = Evals(i);
    for j = 1:length(startIdx)
        start = startIdx(j);
        fprintf('E = %d, iteration %d\n', E, j);
        
        Xtest = X(start:start+L+extra-1);
        Mx = makeShadowManifold(Xtest, E, tau);
        Mx = Mx(end-L+1:end,:);
        [Xtrue, Xest] = testShadowManifold( Mx, E, tau);
        
        R = corrcoef(Xtrue, Xest);
        XxmapX(j,i) = R(1,2);
        
        Ytest = Y(start:start+L+extra-1);
        My = makeShadowManifold(Ytest, E, tau);
        My = My(end-L+1:end,:);
        Yest = crossMap(Mx, My, E);
        Ytrue = Ytest(end-length(Yest)+1:end);
        
        R = corrcoef(Ytrue, Yest);
        XxmapY(j,i) = R(1,2);
    
    end
end

% make figure
fig = figure();
set (fig, 'Units', 'normalized', 'Position', [0,0,0.9,0.9], 'Color', [1,1,1]);
hold on
x=Evals;
y1=XxmapX;
y2=XxmapY;
h1=shadedErrorBar(x,median(y1),[abs(prctile(y1,75)-median(y1)); abs(prctile(y1,25)-median(y1))],'g',1);
h2=shadedErrorBar(x,median(y2),[abs(prctile(y2,75)-median(y2)); abs(prctile(y2,25)-median(y2))],'b',1);
xlim([minE maxE]) 
box on
title({'Accuracy of predictions using time-lagged manifolds of different dimensions',sprintf('(L=%d, \\tau = %d)',L, tau)})
xlabel('E (embedding dimension)')
ylabel('\rho')
legend([h1.mainLine, h2.mainLine],'X xmap X', 'X xmap Y', 'Location', 'best')
hold off


end