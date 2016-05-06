%% dimension plot for ACC vs Am

fig = figure();
set (fig, 'Units', 'normalized', 'Position', [0,0,0.9,0.9], 'Color', [1,1,1]);
hold on
x=Evals;
y1=CxmapC;
y2=CxmapA;
h1=shadedErrorBar(x,median(y1),[abs(prctile(y1,75)-median(y1)); abs(prctile(y1,25)-median(y1))],'g',1);
h2=shadedErrorBar(x,median(y2),[abs(prctile(y2,75)-median(y2)); abs(prctile(y2,25)-median(y2))],'b',1);
xlim([min(Evals) max(Evals)]);
box on
title({'Embedding dimension vs. prediction accuracy for neural data',sprintf('(EC108, L=%d, \\tau = %d)',LN, tauN)})
xlabel('E')
ylabel('\rho')
legend([h1.mainLine, h2.mainLine],'ACC leave-one-out','ACC xmap Am', 'Location', 'best')
hold off

print('pngs/EC108_dimensionPlot','-dpng','-r600')
savefig('figs/EC108_dimensionPlot.fig')


%% dimension plot for difference equations system

fig = figure();
set (fig, 'Units', 'normalized', 'Position', [0,0,0.9,0.9], 'Color', [1,1,1]);
hold on
x=EvalsD;
y1=XxmapX;
y2=XxmapY;
h1=shadedErrorBar(x,median(y1),[abs(prctile(y1,75)-median(y1)); abs(prctile(y1,25)-median(y1))],'g',1);
h2=shadedErrorBar(x,median(y2),[abs(prctile(y2,75)-median(y2)); abs(prctile(y2,25)-median(y2))],'b',1);
box on
title({'Embedding dimension vs. prediction accuracy for difference equations',sprintf('(L=%d, \\tau = %d)',LD, tauD)})
xlabel('E')
ylabel('\rho')
legend([h1.mainLine, h2.mainLine],'X leave-one-out','X xmap Y', 'Location', 'best')
hold off

print('pngs/DiffEqs_dimensionPlot','-dpng','-r600')
savefig('figs/DiffEqs_dimensionPlot.fig')

%%

acc = acc_full;
am = am_full;
ofc = ofc_full;
hp = hp_full;

%% create filtered versions of the data

butterBeta = designfilt('bandpassiir','FilterOrder',4, ...
'HalfPowerFrequency1',13,'HalfPowerFrequency2',30, ...
'SampleRate',512);

butterGamma = designfilt('bandpassiir','FilterOrder',4, ...
'HalfPowerFrequency1',32,'HalfPowerFrequency2',126, ...
'SampleRate',512);

butterHighGamma = designfilt('bandpassiir','FilterOrder',4, ...
'HalfPowerFrequency1',128,'HalfPowerFrequency2',254, ...
'SampleRate',512);

acc_beta = filtfilt(butterBeta, acc);
acc_gamma = filtfilt(butterGamma, acc);
acc_highgamma = filtfilt(butterHighGamma, acc);

am_beta = filtfilt(butterBeta, am);
am_gamma = filtfilt(butterGamma, am);
am_highgamma = filtfilt(butterHighGamma, am);

ofc_beta = filtfilt(butterBeta, ofc);
ofc_gamma = filtfilt(butterGamma,ofc);
ofc_highgamma = filtfilt(butterHighGamma, ofc);

hp_beta = filtfilt(butterBeta, hp);
hp_gamma = filtfilt(butterGamma, hp);
hp_highgamma = filtfilt(butterHighGamma, hp);

%% run testE for each band

[Evals_b, CxmapC_beta, CxmapA_beta] = testE(2,2,12,acc_beta,am_beta,30,500,10);
[Evals_g, CxmapC_gamma, CxmapA_gamma] = testE(2,2,12,acc_gamma,am_gamma,30,500,10);
[Evals_hg, CxmapC_highgamma, CxmapA_highgamma] = testE(2,2,12,acc_highgamma,am_highgamma,30,500,10);

%%

LN = 500;
tauN=30;

%% dimension plot for beta band

fig = figure();
set (fig, 'Units', 'normalized', 'Position', [0,0,0.9,0.9], 'Color', [1,1,1]);
hold on
x=Evals_b;
y1=CxmapC_beta;
y2=CxmapA_beta;
h1=shadedErrorBar(x,median(y1),[abs(prctile(y1,75)-median(y1)); abs(prctile(y1,25)-median(y1))],'g',1);
h2=shadedErrorBar(x,median(y2),[abs(prctile(y2,75)-median(y2)); abs(prctile(y2,25)-median(y2))],'b',1);
xlim([min(Evals_b) max(Evals_b)]);
box on
title({'Embedding dimension vs. prediction accuracy for beta band (13-30 Hz)',sprintf('(EC108, L=%d, \\tau = %d)',LN, tauN)})
xlabel('E')
ylabel('\rho (median over 10 segments)')
legend([h1.mainLine, h2.mainLine],'ACC leave-one-out','ACC xmap Am', 'Location', 'best')
hold off

print('pngs/EC108_dimensionPlot_beta','-dpng','-r600')
savefig('figs/EC108_dimensionPlot_beta.fig')

%% dimension plot for gamma band

fig = figure();
set (fig, 'Units', 'normalized', 'Position', [0,0,0.9,0.9], 'Color', [1,1,1]);
hold on
x=Evals_g;
y1=CxmapC_gamma;
y2=CxmapA_gamma;
h1=shadedErrorBar(x,median(y1),[abs(prctile(y1,75)-median(y1)); abs(prctile(y1,25)-median(y1))],'g',1);
h2=shadedErrorBar(x,median(y2),[abs(prctile(y2,75)-median(y2)); abs(prctile(y2,25)-median(y2))],'b',1);
xlim([min(Evals_g) max(Evals_g)]);
box on
title({'Embedding dimension vs. prediction accuracy for gamma band (32-126 Hz)',sprintf('(EC108, L=%d, \\tau = %d)',LN, tauN)})
xlabel('E')
ylabel('\rho (median over 10 segments)')
legend([h1.mainLine, h2.mainLine],'ACC leave-one-out','ACC xmap Am', 'Location', 'best')
hold off

print('pngs/EC108_dimensionPlot_gamma','-dpng','-r600')
savefig('figs/EC108_dimensionPlot_gamma.fig')

%% dimension plot for high gamma band

fig = figure();
set (fig, 'Units', 'normalized', 'Position', [0,0,0.9,0.9], 'Color', [1,1,1]);
hold on
x=Evals_hg;
y1=CxmapC_highgamma;
y2=CxmapA_highgamma;
h1=shadedErrorBar(x,median(y1),[abs(prctile(y1,75)-median(y1)); abs(prctile(y1,25)-median(y1))],'g',1);
h2=shadedErrorBar(x,median(y2),[abs(prctile(y2,75)-median(y2)); abs(prctile(y2,25)-median(y2))],'b',1);
xlim([min(Evals_hg) max(Evals_hg)]);
box on
title({'Embedding dimension vs. prediction accuracy for high gamma band (128-256Hz)',sprintf('(EC108, L=%d, \\tau = %d)',LN, tauN)})
xlabel('E')
ylabel('\rho (median over 10 segments)')
legend([h1.mainLine, h2.mainLine],'ACC leave-one-out','ACC xmap Am', 'Location', 'best')
hold off

print('pngs/EC108_dimensionPlot_highgamma','-dpng','-r600')
savefig('figs/EC108_dimensionPlot_highgamma.fig')


%% do CCM on gamma
acc = acc_gamma;
ofc = ofc_gamma;
hp = hp_gamma;
am = am_gamma;

% set params and initialize variables

E = 10;
tau = 30;
minL = (tau +1)*E + 2 - tau;
n = 50;
maxL =5000;

startIdx = 1:maxL:length(ofc);
startIdx = startIdx(1:n);
maxL = 500;
Lvals = [minL:5:maxL, maxL];

OxmapA = NaN(n, length(Lvals));
AxmapO = NaN(n, length(Lvals));
OxmapH = NaN(n, length(Lvals));
HxmapO = NaN(n, length(Lvals));
AxmapH = NaN(n, length(Lvals));
HxmapA = NaN(n, length(Lvals));
CxmapO = NaN(n, length(Lvals));
OxmapC = NaN(n, length(Lvals));
CxmapA = NaN(n, length(Lvals));
AxmapC = NaN(n, length(Lvals));
CxmapH = NaN(n, length(Lvals));
HxmapC = NaN(n, length(Lvals));


% do CCM

for i = 1:length(Lvals)
    for j = 1:length(startIdx)
        
        L = Lvals(i);
        start = startIdx(j);
        
        fprintf('L = %d, iteration %d\n', L, j); 
        % xmap OFC & Am
        X = ofc(start:start+L-1)';
        Y = am(start:start+L-1)';
        [ Xest, Yest, ~, ~, ~, ~ ] = CCM( X, Y, E, tau );
        
        
        R = corrcoef(X(size(X,1)-size(Xest,1)+1:end), Xest);
        AxmapO(j,i) = R(1,2);
        
        R = corrcoef(Y(size(X,1)-size(Yest,1)+1:end), Yest);
        OxmapA(j,i) = R(1,2);
        
        % xmap OFC & Hp
        X = ofc(start:start+L-1)';
        Y = hp(start:start+L-1)';
        [ Xest, Yest, ~, ~, ~, ~ ] = CCM( X, Y, E, tau );
        
        R = corrcoef(X(size(X,1)-size(Xest,1)+1:end), Xest);
        HxmapO(j,i) = R(1,2);
        
        R = corrcoef(Y(size(X,1)-size(Yest,1)+1:end), Yest);
        OxmapH(j,i) = R(1,2);
        
        % xmap Hp & Am
        X = hp(start:start+L-1)';
        Y = am(start:start+L-1)';
        [ Xest, Yest, ~, ~, ~, ~ ] = CCM( X, Y, E, tau );
        
        R = corrcoef(X(size(X,1)-size(Xest,1)+1:end), Xest);
        AxmapH(j,i) = R(1,2);
        
        R = corrcoef(Y(size(X,1)-size(Yest,1)+1:end), Yest);
        HxmapA(j,i) = R(1,2);
        
        % xmap ACC & OFC
        X = acc(start:start+L-1)';
        Y = ofc(start:start+L-1)';
        [ Xest, Yest, ~, ~, ~, ~ ] = CCM( X, Y, E, tau );
        
        R = corrcoef(X(size(X,1)-size(Xest,1)+1:end), Xest);
        OxmapC(j,i) = R(1,2);
        
        R = corrcoef(Y(size(X,1)-size(Yest,1)+1:end), Yest);
        CxmapO(j,i) = R(1,2);
        
        % xmap ACC & Am
        X = acc(start:start+L-1)';
        Y = am(start:start+L-1)';
        [ Xest, Yest, ~, ~, ~, ~ ] = CCM( X, Y, E, tau );
        
        R = corrcoef(X(size(X,1)-size(Xest,1)+1:end), Xest);
        AxmapC(j,i) = R(1,2);
        
        R = corrcoef(Y(size(X,1)-size(Yest,1)+1:end), Yest);
        CxmapA(j,i) = R(1,2);
        
        % xmap ACC & Hp
        X = acc(start:start+L-1)';
        Y = hp(start:start+L-1)';
        [ Xest, Yest, ~, ~, ~, ~ ] = CCM( X, Y, E, tau );
        
        R = corrcoef(X(size(X,1)-size(Xest,1)+1:end), Xest);
        HxmapC(j,i) = R(1,2);
        
        R = corrcoef(Y(size(X,1)-size(Yest,1)+1:end), Yest);
        CxmapH(j,i) = R(1,2);
        
    end
    save('20160501_EC108_CCM_E10_tau10_maxL600_gamma.mat')
end

% OFC vs Am
fig = figure();
set (fig, 'Units', 'normalized', 'Position', [0,0,0.9,0.9], 'Color', [1,1,1]);
hold on
x=Lvals;
y1=AxmapO;
y2=OxmapA;
h1=shadedErrorBar(x,median(y1),[abs(prctile(y1,75)-median(y1)); abs(prctile(y1,25)-median(y1))],'g',1);
h2=shadedErrorBar(x,median(y2),[abs(prctile(y2,75)-median(y2)); abs(prctile(y2,25)-median(y2))],'b',1);
xlim([minL maxL]) %maxL])
%ylim([min([min(y1), min(y2)]), 1.05])
box on
title({'CCM for gamma band in OFC and Am','(EC108, E=10, \tau = 30)'})
xlabel('L')
ylabel('\rho')
legend([h1.mainLine, h2.mainLine],'Am xmap OFC','OFC xmap Am', 'Location', 'best')
hold off

print('pngs/EC108_CCM_OFCAm_E10_gamma','-dpng','-r600')
savefig('figs/EC108_CCM_OFCAm_E10_gamma.fig')

% OFC vs Hp
fig = figure();
set (fig, 'Units', 'normalized', 'Position', [0,0,0.9,0.9], 'Color', [1,1,1]);
hold on
x=Lvals;
y1=HxmapO;
y2=OxmapH;
h1=shadedErrorBar(x,median(y1),[abs(prctile(y1,75)-median(y1)); abs(prctile(y1,25)-median(y1))],'g',1);
h2=shadedErrorBar(x,median(y2),[abs(prctile(y2,75)-median(y2)); abs(prctile(y2,25)-median(y2))],'b',1);
xlim([minL maxL]) %maxL])
%ylim([min([min(y1), min(y2)]), 1.05])
box on
title({'CCM for gamma band in OFC and Hp','(EC108, E=10, \tau = 30)'})
xlabel('L')
ylabel('\rho')
legend([h1.mainLine, h2.mainLine],'Hp xmap OFC','OFC xmap Hp', 'Location', 'best')
hold off

print('pngs/EC108_CCM_OFCHp_E10_gamma','-dpng','-r600')
savefig('figs/EC108_CCM_OFCHp_E10_gamma.fig')

% OFC vs ACC
fig = figure();
set (fig, 'Units', 'normalized', 'Position', [0,0,0.9,0.9], 'Color', [1,1,1]);
hold on
x=Lvals;
y1=CxmapO;
y2=OxmapC;
h1=shadedErrorBar(x,median(y1),[abs(prctile(y1,75)-median(y1)); abs(prctile(y1,25)-median(y1))],'g',1);
h2=shadedErrorBar(x,median(y2),[abs(prctile(y2,75)-median(y2)); abs(prctile(y2,25)-median(y2))],'b',1);
xlim([minL maxL]) %maxL])
%ylim([min([min(y1), min(y2)]), 1.05])
box on
title({'CCM for gamma band in OFC and ACC','(EC108, E=10, \tau = 30)'})
xlabel('L')
ylabel('\rho')
legend([h1.mainLine, h2.mainLine],'ACC xmap OFC','OFC xmap ACC', 'Location', 'best')
hold off

print('pngs/EC108_CCM_OFCACC_E10_gamma','-dpng','-r600')
savefig('figs/EC108_CCM_OFCACC_E10_gamma.fig')

% ACC vs Am
fig = figure();
set (fig, 'Units', 'normalized', 'Position', [0,0,0.9,0.9], 'Color', [1,1,1]);
hold on
x=Lvals;
y1=AxmapC;
y2=CxmapA;
h1=shadedErrorBar(x,median(y1),[abs(prctile(y1,75)-median(y1)); abs(prctile(y1,25)-median(y1))],'g',1);
h2=shadedErrorBar(x,median(y2),[abs(prctile(y2,75)-median(y2)); abs(prctile(y2,25)-median(y2))],'b',1);
xlim([minL maxL]) %maxL])
%ylim([min([min(y1), min(y2)]), 1.05])
box on
title({'CCM for gamma band in ACC and Am','(EC108, E=10, \tau = 30)'})
xlabel('L')
ylabel('\rho')
legend([h1.mainLine, h2.mainLine],'Am xmap ACC','ACC xmap Am', 'Location', 'best')
hold off

print('pngs/EC108_CCM_ACCAm_E10_gamma','-dpng','-r600')
savefig('figs/EC108_CCM_ACCAm_E10_gamma.fig')


% ACC vs Hp
fig = figure();
set (fig, 'Units', 'normalized', 'Position', [0,0,0.9,0.9], 'Color', [1,1,1]);
hold on
x=Lvals;
y1=HxmapC;
y2=CxmapH;
h1=shadedErrorBar(x,median(y1),[abs(prctile(y1,75)-median(y1)); abs(prctile(y1,25)-median(y1))],'g',1);
h2=shadedErrorBar(x,median(y2),[abs(prctile(y2,75)-median(y2)); abs(prctile(y2,25)-median(y2))],'b',1);
xlim([minL maxL]) %maxL])
%ylim([min([min(y1), min(y2)]), 1.05])
box on
title({'CCM for gamma band in ACC and Hp','(EC108, E=10, \tau = 30)'})
xlabel('L')
ylabel('\rho')
legend([h1.mainLine, h2.mainLine],'Hp xmap ACC','ACC xmap Hp', 'Location', 'best')
hold off

print('pngs/EC108_CCM_ACCHp_E10_gammaa','-dpng','-r600')
savefig('figs/EC108_CCM_ACCHp_E10_gamma.fig')

% Am vs Hp
fig = figure();
set (fig, 'Units', 'normalized', 'Position', [0,0,0.9,0.9], 'Color', [1,1,1]);
hold on
x=Lvals;
y1=HxmapA;
y2=AxmapH;
h1=shadedErrorBar(x,median(y1),[abs(prctile(y1,75)-median(y1)); abs(prctile(y1,25)-median(y1))],'g',1);
h2=shadedErrorBar(x,median(y2),[abs(prctile(y2,75)-median(y2)); abs(prctile(y2,25)-median(y2))],'b',1);
xlim([minL maxL]) %maxL])
%ylim([min([min(y1), min(y2)]), 1.05])
box on
title({'CCM for gamma band in Am and Hp','(EC108, E=10, \tau = 30)'})
xlabel('L')
ylabel('\rho')
legend([h1.mainLine, h2.mainLine],'Hp xmap Am','Am xmap Hp', 'Location', 'best')
hold off

print('pngs/EC108_CCM_AmHp_E10_gamma','-dpng','-r600')
savefig('figs/EC108_CCM_AmHp_E10_gamma.fig')

%%
close all
%% okay now do high gamma

acc = acc_highgamma;
ofc = ofc_highgamma;
hp = hp_highgamma;
am = am_highgamma;

% set params and initialize variables

E = 10;
tau = 30;
minL = (tau +1)*E + 2 - tau;
n = 50;
maxL =5000;

startIdx = 1:maxL:length(ofc);
startIdx = startIdx(1:n);
maxL = 500;
Lvals = [minL:5:maxL, maxL];

OxmapA = NaN(n, length(Lvals));
AxmapO = NaN(n, length(Lvals));
OxmapH = NaN(n, length(Lvals));
HxmapO = NaN(n, length(Lvals));
AxmapH = NaN(n, length(Lvals));
HxmapA = NaN(n, length(Lvals));
CxmapO = NaN(n, length(Lvals));
OxmapC = NaN(n, length(Lvals));
CxmapA = NaN(n, length(Lvals));
AxmapC = NaN(n, length(Lvals));
CxmapH = NaN(n, length(Lvals));
HxmapC = NaN(n, length(Lvals));


% do CCM

for i = 1:length(Lvals)
    for j = 1:length(startIdx)
        
        L = Lvals(i);
        start = startIdx(j);
        
        fprintf('L = %d, iteration %d\n', L, j); 
        % xmap OFC & Am
        X = ofc(start:start+L-1)';
        Y = am(start:start+L-1)';
        [ Xest, Yest, ~, ~, ~, ~ ] = CCM( X, Y, E, tau );
        
        R = corrcoef(X(size(X,1)-size(Xest,1)+1:end), Xest);
        AxmapO(j,i) = R(1,2);
        
        R = corrcoef(Y(size(X,1)-size(Yest,1)+1:end), Yest);
        OxmapA(j,i) = R(1,2);
        
        % xmap OFC & Hp
        X = ofc(start:start+L-1)';
        Y = hp(start:start+L-1)';
        [ Xest, Yest, ~, ~, ~, ~ ] = CCM( X, Y, E, tau );
        
        R = corrcoef(X(size(X,1)-size(Xest,1)+1:end), Xest);
        HxmapO(j,i) = R(1,2);
        
        R = corrcoef(Y(size(X,1)-size(Yest,1)+1:end), Yest);
        OxmapH(j,i) = R(1,2);
        
        % xmap Hp & Am
        X = hp(start:start+L-1)';
        Y = am(start:start+L-1)';
        [ Xest, Yest, ~, ~, ~, ~ ] = CCM( X, Y, E, tau );
        
        R = corrcoef(X(size(X,1)-size(Xest,1)+1:end), Xest);
        AxmapH(j,i) = R(1,2);
        
        R = corrcoef(Y(size(X,1)-size(Yest,1)+1:end), Yest);
        HxmapA(j,i) = R(1,2);
        
        % xmap ACC & OFC
        X = acc(start:start+L-1)';
        Y = ofc(start:start+L-1)';
        [ Xest, Yest, ~, ~, ~, ~ ] = CCM( X, Y, E, tau );
        
        R = corrcoef(X(size(X,1)-size(Xest,1)+1:end), Xest);
        OxmapC(j,i) = R(1,2);
        
        R = corrcoef(Y(size(X,1)-size(Yest,1)+1:end), Yest);
        CxmapO(j,i) = R(1,2);
        
        % xmap ACC & Am
        X = acc(start:start+L-1)';
        Y = am(start:start+L-1)';
        [ Xest, Yest, ~, ~, ~, ~ ] = CCM( X, Y, E, tau );
        
        R = corrcoef(X(size(X,1)-size(Xest,1)+1:end), Xest);
        AxmapC(j,i) = R(1,2);
        
        R = corrcoef(Y(size(X,1)-size(Yest,1)+1:end), Yest);
        CxmapA(j,i) = R(1,2);
        
        % xmap ACC & Hp
        X = acc(start:start+L-1)';
        Y = hp(start:start+L-1)';
        [ Xest, Yest, ~, ~, ~, ~ ] = CCM( X, Y, E, tau );
        
        R = corrcoef(X(size(X,1)-size(Xest,1)+1:end), Xest);
        HxmapC(j,i) = R(1,2);
        
        R = corrcoef(Y(size(X,1)-size(Yest,1)+1:end), Yest);
        CxmapH(j,i) = R(1,2);
        
    end
    save('20160501_EC108_CCM_E10_tau10_maxL600_highgamma.mat')
end

% OFC vs Am
fig = figure();
set (fig, 'Units', 'normalized', 'Position', [0,0,0.9,0.9], 'Color', [1,1,1]);
hold on
x=Lvals;
y1=AxmapO;
y2=OxmapA;
h1=shadedErrorBar(x,median(y1),[abs(prctile(y1,75)-median(y1)); abs(prctile(y1,25)-median(y1))],'g',1);
h2=shadedErrorBar(x,median(y2),[abs(prctile(y2,75)-median(y2)); abs(prctile(y2,25)-median(y2))],'b',1);
xlim([minL maxL]) %maxL])
%ylim([min([min(y1), min(y2)]), 1.05])
box on
title({'CCM for high gamma band in OFC and Am','(EC108, E=10, \tau = 30)'})
xlabel('L')
ylabel('\rho')
legend([h1.mainLine, h2.mainLine],'Am xmap OFC','OFC xmap Am', 'Location', 'best')
hold off

print('pngs/EC108_CCM_OFCAm_E10_highgamma','-dpng','-r600')
savefig('figs/EC108_CCM_OFCAm_E10_highgamma.fig')

% OFC vs Hp
fig = figure();
set (fig, 'Units', 'normalized', 'Position', [0,0,0.9,0.9], 'Color', [1,1,1]);
hold on
x=Lvals;
y1=HxmapO;
y2=OxmapH;
h1=shadedErrorBar(x,median(y1),[abs(prctile(y1,75)-median(y1)); abs(prctile(y1,25)-median(y1))],'g',1);
h2=shadedErrorBar(x,median(y2),[abs(prctile(y2,75)-median(y2)); abs(prctile(y2,25)-median(y2))],'b',1);
xlim([minL maxL]) %maxL])
%ylim([min([min(y1), min(y2)]), 1.05])
box on
title({'CCM for high gamma band in OFC and Hp','(EC108, E=10, \tau = 30)'})
xlabel('L')
ylabel('\rho')
legend([h1.mainLine, h2.mainLine],'Hp xmap OFC','OFC xmap Hp', 'Location', 'best')
hold off

print('pngs/EC108_CCM_OFCHp_E10_highgamma','-dpng','-r600')
savefig('figs/EC108_CCM_OFCHp_E10_highgamma.fig')

% OFC vs ACC
fig = figure();
set (fig, 'Units', 'normalized', 'Position', [0,0,0.9,0.9], 'Color', [1,1,1]);
hold on
x=Lvals;
y1=CxmapO;
y2=OxmapC;
h1=shadedErrorBar(x,median(y1),[abs(prctile(y1,75)-median(y1)); abs(prctile(y1,25)-median(y1))],'g',1);
h2=shadedErrorBar(x,median(y2),[abs(prctile(y2,75)-median(y2)); abs(prctile(y2,25)-median(y2))],'b',1);
xlim([minL maxL]) %maxL])
%ylim([min([min(y1), min(y2)]), 1.05])
box on
title({'CCM for high gamma band in OFC and ACC','(EC108, E=10, \tau = 30)'})
xlabel('L')
ylabel('\rho')
legend([h1.mainLine, h2.mainLine],'ACC xmap OFC','OFC xmap ACC', 'Location', 'best')
hold off

print('pngs/EC108_CCM_OFCACC_E10_highgamma','-dpng','-r600')
savefig('figs/EC108_CCM_OFCACC_E10_highgamma.fig')

% ACC vs Am
fig = figure();
set (fig, 'Units', 'normalized', 'Position', [0,0,0.9,0.9], 'Color', [1,1,1]);
hold on
x=Lvals;
y1=AxmapC;
y2=CxmapA;
h1=shadedErrorBar(x,median(y1),[abs(prctile(y1,75)-median(y1)); abs(prctile(y1,25)-median(y1))],'g',1);
h2=shadedErrorBar(x,median(y2),[abs(prctile(y2,75)-median(y2)); abs(prctile(y2,25)-median(y2))],'b',1);
xlim([minL maxL]) %maxL])
%ylim([min([min(y1), min(y2)]), 1.05])
box on
title({'CCM for high gamma band in ACC and Am','(EC108, E=10, \tau = 30)'})
xlabel('L')
ylabel('\rho')
legend([h1.mainLine, h2.mainLine],'Am xmap ACC','ACC xmap Am', 'Location', 'best')
hold off

print('pngs/EC108_CCM_ACCAm_E10_highgamma','-dpng','-r600')
savefig('figs/EC108_CCM_ACCAm_E10_highgamma.fig')


% ACC vs Hp
fig = figure();
set (fig, 'Units', 'normalized', 'Position', [0,0,0.9,0.9], 'Color', [1,1,1]);
hold on
x=Lvals;
y1=HxmapC;
y2=CxmapH;
h1=shadedErrorBar(x,median(y1),[abs(prctile(y1,75)-median(y1)); abs(prctile(y1,25)-median(y1))],'g',1);
h2=shadedErrorBar(x,median(y2),[abs(prctile(y2,75)-median(y2)); abs(prctile(y2,25)-median(y2))],'b',1);
xlim([minL maxL]) %maxL])
%ylim([min([min(y1), min(y2)]), 1.05])
box on
title({'CCM for high gamma band in ACC and Hp','(EC108, E=10, \tau = 30)'})
xlabel('L')
ylabel('\rho')
legend([h1.mainLine, h2.mainLine],'Hp xmap ACC','ACC xmap Hp', 'Location', 'best')
hold off

print('pngs/EC108_CCM_ACCHp_E10_highgammaa','-dpng','-r600')
savefig('figs/EC108_CCM_ACCHp_E10_highgamma.fig')

% Am vs Hp
fig = figure();
set (fig, 'Units', 'normalized', 'Position', [0,0,0.9,0.9], 'Color', [1,1,1]);
hold on
x=Lvals;
y1=HxmapA;
y2=AxmapH;
h1=shadedErrorBar(x,median(y1),[abs(prctile(y1,75)-median(y1)); abs(prctile(y1,25)-median(y1))],'g',1);
h2=shadedErrorBar(x,median(y2),[abs(prctile(y2,75)-median(y2)); abs(prctile(y2,25)-median(y2))],'b',1);
xlim([minL maxL]) %maxL])
%ylim([min([min(y1), min(y2)]), 1.05])
box on
title({'CCM for high gamma band in Am and Hp','(EC108, E=10, \tau = 30)'})
xlabel('L')
ylabel('\rho')
legend([h1.mainLine, h2.mainLine],'Hp xmap Am','Am xmap Hp', 'Location', 'best')
hold off

print('pngs/EC108_CCM_AmHp_E10_highgamma','-dpng','-r600')
savefig('figs/EC108_CCM_AmHp_E10_highgamma.fig')
































