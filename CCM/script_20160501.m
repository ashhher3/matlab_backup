%% load data and trim

idx = [65, 96; 149, 158; 159, 168; 179, 188];
idx = idx-64;
ofc = []; am=[]; hp=[]; acc=[];

for f=1:4
load(sprintf('EC108_2fd7f1ac_all_%03d.mat',f))
ofc = [ofc, mean(data(idx(1,1):idx(1,2),:))];
am = [am, mean(data(idx(2,1):idx(2,2),:))];
hp = [hp, mean(data(idx(3,1):idx(3,2),:))];
acc = [acc, mean(data(idx(4,1):idx(4,2),:))];
end

last = 5000*100;
ofc = ofc(1:last);
am = am(1:last);
hp = hp(1:last);
acc = acc(1:last);

clear idx last f data

%% filter out beta

butterBeta = designfilt('bandpassiir','FilterOrder',4, ...
'HalfPowerFrequency1',13,'HalfPowerFrequency2',30, ...
'SampleRate',512);

acc = filtfilt(butterBeta, acc);
ofc = filtfilt(butterBeta, ofc);
am = filtfilt(butterBeta, am);
hp = filtfilt(butterBeta, hp);


%% set params and initialize variables

E = 10;
tau = 30;
minL = (tau +1)*E + 2 - tau;
n = 100;
maxL =5000;

startIdx = 1:maxL:length(ofc);
maxL = 600;
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


%% do CCM

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
    save('20160501_EC108_CCM_E10_tau10_maxL600_beta.mat')
end

matlabmail('kderosier6@gmail.com', 'CCM is done', 'CCM is done, beta version')

%% make plots from CCM data

%% OFC vs Am
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
title({'CCM for beta band in OFC and Am','(EC108, E=10, \tau = 30)'})
xlabel('L')
ylabel('\rho')
legend([h1.mainLine, h2.mainLine],'Am xmap OFC','OFC xmap Am', 'Location', 'best')
hold off

print('pngs/EC108_CCM_OFCAm_E10_beta','-dpng','-r600')
savefig('figs/EC108_CCM_OFCAm_E10_beta.fig')

%% OFC vs Hp
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
title({'CCM for beta band in OFC and Hp','(EC108, E=10, \tau = 30)'})
xlabel('L')
ylabel('\rho')
legend([h1.mainLine, h2.mainLine],'Hp xmap OFC','OFC xmap Hp', 'Location', 'best')
hold off

print('pngs/EC108_CCM_OFCHp_E10_beta','-dpng','-r600')
savefig('figs/EC108_CCM_OFCHp_E10_beta.fig')

%% OFC vs ACC
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
title({'CCM for beta band in OFC and ACC','(EC108, E=10, \tau = 30)'})
xlabel('L')
ylabel('\rho')
legend([h1.mainLine, h2.mainLine],'ACC xmap OFC','OFC xmap ACC', 'Location', 'best')
hold off

print('pngs/EC108_CCM_OFCACC_E10_beta','-dpng','-r600')
savefig('figs/EC108_CCM_OFCACC_E10_beta.fig')

%% ACC vs Am
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
title({'CCM for beta band in ACC and Am','(EC108, E=10, \tau = 30)'})
xlabel('L')
ylabel('\rho')
legend([h1.mainLine, h2.mainLine],'Am xmap ACC','ACC xmap Am', 'Location', 'best')
hold off

print('pngs/EC108_CCM_ACCAm_E10_beta','-dpng','-r600')
savefig('figs/EC108_CCM_ACCAm_E10_beta.fig')


%% ACC vs Hp
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
title({'CCM for beta band in ACC and Hp','(EC108, E=10, \tau = 30)'})
xlabel('L')
ylabel('\rho')
legend([h1.mainLine, h2.mainLine],'Hp xmap ACC','ACC xmap Hp', 'Location', 'best')
hold off

print('pngs/EC108_CCM_ACCHp_E10_beta','-dpng','-r600')
savefig('figs/EC108_CCM_ACCHp_E10_beta.fig')

%% Am vs Hp
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
title({'CCM for beta band in Am and Hp','(EC108, E=10, \tau = 30)'})
xlabel('L')
ylabel('\rho')
legend([h1.mainLine, h2.mainLine],'Hp xmap Am','Am xmap Hp', 'Location', 'best')
hold off

print('pngs/EC108_CCM_AmHp_E10_beta','-dpng','-r600')
savefig('figs/EC108_CCM_AmHp_E10_beta.fig')
