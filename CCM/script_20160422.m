%% make fake data from the coupled difference equation system

% set params of system
% X(t+1) = X(t)[rx - rx*X(t) - Bxy*Y(t)]
% Y(t+1) = Y(t)[ry - ry*Y(t) - Byx*X(t)]
rx = 3;
Bxy =0.3;
ry = 3.9;
Byx = 0.2;

% initialize
X = NaN(5000, 1);
Y = NaN(5000, 1);
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
    [ Xself, Yest, Xest, Mx, My ] = CCM( X, Y, E, tau );
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
tau = 4;

minL = 20;
step = 10;
maxL = length(X);

rhoXxmapY = [];
rhoYxmapX = [];

for L = minL:step:maxL
    
    Xtest = X(1:L);
    Ytest = Y(1:L);
    [ ~, Yest, Xest, ~, ~ ] = CCM( Xtest, Ytest, E, tau );
    
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

%%
clear all
close all