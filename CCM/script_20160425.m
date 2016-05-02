%%

% script for doing CCM with subnets data, original version 20160425 KHPD
% 20160502 adjusted for change in CCM output order - KD

%% load in some data
idx = [65, 96; 149, 158; 159, 168];
idx = idx-64;
ofc = []; am=[]; hp=[];

for f=1:6
load(sprintf('EC108_2fd7f1ac_all_%03d.mat',f))
ofc = [ofc, mean(data(idx(1,1):idx(1,2),:))];
am = [am, mean(data(idx(2,1):idx(2,2),:))];
hp = [hp, mean(data(idx(3,1):idx(3,2),:))];
end

%% first test out a couple values of E, make plots

E = 10;
tau = 20;

L = 10000;

X = ofc(1:L)';
Y = hp(1:L)';

[ Xest, Yest, Xself, Yself, Mx, My ] = CCM( X, Y, E, tau );


figure()
hold on
scatter(Y(size(X,1)-size(Xself,1)+1:end),Yest);
scatter(X(size(X,1)-size(Xself,1)+1:end),Xest);
scatter(X(size(X,1)-size(Xself,1)+1:end),Xself);
plot(X(size(X,1)-size(Xself,1)+1:end),X(size(X,1)-size(Xself,1)+1:end));
title(sprintf('E = %d, tau = %d',E, tau))
xlabel('true value')
ylabel('cross-mapped estimate')
hold off

%% look for convergence

minL = 1000;
minE = 19;
maxE = 30;

for E = minE:maxE
    for tau = [E+2, 2*E, 40];
        minL = (tau +1)*E -2 + tau;
        rhoOFCxmapHp = []; %nan(maxE - minE+1, length(ofc)-minL);
        rhoHpxmapOFC = []; % nan(maxE - minE+1, length(ofc)-minL);
        
        for L = minL:1000:10000; %length(ofc)
            % do xmap with ofc as x and hp as Y
            X = ofc(1:L)';
            Y = hp(1:L)';
            [ Xest, Yest, ~, ~, ~, ~ ] = CCM( X, Y, E, tau );
            
            R = corrcoef(X(size(X,1)-size(Xest,1)+1:end), Xest);
            rhoOFCxmapHp = [rhoOFCxmapHp, R(1,2)];
            
            R = corrcoef(Y(size(X,1)-size(Yest,1)+1:end), Yest);
            rhoHpxmapOFC = [rhoHpxmapOFC, R(1,2)];
            
            fprintf('L = %d\n', L);
        end
        
        figure()
        hold on
        plot(minL:1000:10000,rhoOFCxmapHp, minL:1000:10000, rhoHpxmapOFC)
        xlabel('L')
        ylabel('corr coef')
        title(sprintf('E = %d, tau = %d',E, tau))
        hold off
        pause()
    end
end

