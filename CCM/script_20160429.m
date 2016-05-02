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
title({'CCM for full signal in OFC and Am',sprintf('(EC108, E=%d, \\tau = %d)',E, tau)})
xlabel('L')
ylabel('\rho')
legend([h1.mainLine, h2.mainLine],'Am xmap OFC','OFC xmap Am', 'Location', 'best')
hold off

% print('pngs/EC108_CCM_OFCAm_E10','-dpng','-r600')
% savefig('figs/EC108_CCM_OFCAm_E10.fig')

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
title({'CCM for full signal in OFC and Hp',sprintf('(EC108, E=%d, \\tau = %d)',E, tau)})
xlabel('L')
ylabel('\rho')
legend([h1.mainLine, h2.mainLine],'Hp xmap OFC','OFC xmap Hp', 'Location', 'best')
hold off

% print('pngs/EC108_CCM_OFCHp_E10','-dpng','-r600')
% savefig('figs/EC108_CCM_OFCHp_E10.fig')

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
title({'CCM for full signal in OFC and ACC',sprintf('(EC108, E=%d, \\tau = %d)',E, tau)})
xlabel('L')
ylabel('\rho')
legend([h1.mainLine, h2.mainLine],'ACC xmap OFC','OFC xmap ACC', 'Location', 'best')
hold off

% print('pngs/EC108_CCM_OFCACC_E10','-dpng','-r600')
% savefig('figs/EC108_CCM_OFCACC_E10.fig')

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
title({'CCM for full signal in ACC and Am',sprintf('(EC108, E=%d, \\tau = %d)',E, tau)})
xlabel('L')
ylabel('\rho')
legend([h1.mainLine, h2.mainLine],'Am xmap ACC','ACC xmap Am', 'Location', 'best')
hold off

% print('pngs/EC108_CCM_ACCAm_E10','-dpng','-r600')
% savefig('figs/EC108_CCM_ACCAm_E10.fig')


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
title({'CCM for full signal in ACC and Hp',sprintf('(EC108, E=%d, \\tau = %d)',E, tau)})
xlabel('L')
ylabel('\rho')
legend([h1.mainLine, h2.mainLine],'Hp xmap ACC','ACC xmap Hp', 'Location', 'best')
hold off
% 
% print('pngs/EC108_CCM_ACCHp_E10','-dpng','-r600')
% savefig('figs/EC108_CCM_ACCHp_E10.fig')

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
title({'CCM for full signal in Am and Hp',sprintf('(EC108, E=%d, \\tau = %d)',E, tau)})
xlabel('L')
ylabel('\rho')
legend([h1.mainLine, h2.mainLine],'Hp xmap Am','Am xmap Hp', 'Location', 'best')
hold off

% print('pngs/EC108_CCM_AmHp_E10','-dpng','-r600')
% savefig('figs/EC108_CCM_AmHp_E10.fig')










































