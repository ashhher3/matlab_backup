%% make fake data from the coupled difference equation system

% set params of system
% X(t+1) = X(t)[rx - rx*X(t) - Bxy*Y(t)]
% Y(t+1) = Y(t)[ry - ry*Y(t) - Byx*X(t)]
rx = 3;
Bxy =0.3;
ry = 3.9;
Byx = 0.2;

% initialize
X = NaN(50000, 1);
Y = NaN(50000, 1);
X(1) = 0.2;
Y(1) = 0.4;

% fill in data
for t=2:length(X)
    X(t) = X(t-1)*[rx - rx*X(t-1) - Bxy*Y(t-1)];
    Y(t) = Y(t-1)*[ry - ry*Y(t-1) - Byx*X(t-1)];
end

clear rx Bxy ry Byx t

%% set params and call CCM

E = 2;

for tau = 2*E; %1:10
    [ Xest, Yest, Xself, Yself, Mx, My ] = CCM( X, Y, E, tau );
    figure()
    hold on
    scatter(Y(size(X,1)-size(Xself,1)+1:end),Yest);
    scatter(X(size(X,1)-size(Xself,1)+1:end),Xest);
    scatter(X(size(X,1)-size(Xself,1)+1:end),Xself);
    plot(X(size(X,1)-size(Xself,1)+1:end),X(size(X,1)-size(Xself,1)+1:end));
    title(sprintf('tau = %d',tau))
    xlabel('true')
    ylabel('cross-mapped estimate')
    hold off
end

%% look for convergence
E = 2;
tau = 2;

minL = 20;
step = 10;
maxL = 5000; %length(X);

rhoXxmapY = [];
rhoYxmapX = [];

for L = minL:step:maxL
    
    Xtest = X(1:L);
    Ytest = Y(1:L);
    [ Xest, Yest, ~, ~, ~, ~ ] = CCM( Xtest, Ytest, E, tau );
    
    R = corrcoef(Xtest(size(Xtest,1)-size(Xest,1)+1:end), Xest);
    rhoXxmapY = [rhoXxmapY, R(1,2)];
    
    R = corrcoef(Ytest(size(Xtest,1)-size(Yest,1)+1:end), Yest);
    rhoYxmapX = [rhoYxmapX, R(1,2)];
    
    fprintf('L = %d\n', L);
end

figure()
hold on
plot(minL:step:maxL,rhoXxmapY,minL:step:maxL, rhoYxmapX)
xlabel('L')
ylabel('corr coef')
title(sprintf('E = %d, tau = %d',E, tau))
hold off
pause()

%% (added 20160502) run over multiple data segments

E = 2;
tau = 4;
minL = (tau +1)*E + 2 - tau;
maxL =1000;
n=length(X)/maxL;

startIdx = 1:maxL:length(X);

Lvals = [minL:5:maxL, maxL];

XxmapY = NaN(n, length(Lvals));
YxmapX = NaN(n, length(Lvals));

for i = 1:length(Lvals)
    for j = 1:length(startIdx)
        
        L = Lvals(i);
        start = startIdx(j);
        
        fprintf('L = %d, iteration %d\n', L, j);
        
        Xtest = X(start:start+L-1);
        Ytest = Y(start:start+L-1);
        
        [ Xest, Yest, ~, ~, ~, ~ ] = CCM( Xtest, Ytest, E, tau );
        
        R = corrcoef(Xtest(size(Xtest,1)-size(Xest,1)+1:end), Xest);
        YxmapX(j,i) = R(1,2);
        
        R = corrcoef(Ytest(size(Xtest,1)-size(Yest,1)+1:end), Yest);
        XxmapY(j,i) = R(1,2);
        
    end
end


%% make the plot
fig = figure();
set (fig, 'Units', 'normalized', 'Position', [0,0,0.9,0.9], 'Color', [1,1,1]);
hold on
x=Lvals;
y1=XxmapY;
y2=YxmapX;
h1=shadedErrorBar(x,median(y1),[abs(prctile(y1,75)-median(y1)); abs(prctile(y1,25)-median(y1))],'g',1);
h2=shadedErrorBar(x,median(y2),[abs(prctile(y2,75)-median(y2)); abs(prctile(y2,25)-median(y2))],'b',1);
% xlim([minL maxL]) %maxL])
% %ylim([min([min(y1), min(y2)]), 1.05])
box on
title({'CCM for difference equations system'})
xlabel('L')
ylabel('\rho')
legend([h1.mainLine, h2.mainLine],'X xmap Y','Y xmap X', 'Location', 'best')
hold off




%%
clear all
close all