function fig = plotCCM(Lvals,XxmapY,YxmapX,E, tau,dataName,xName,yName,saveName)
%% plotCCM.m Makes convergent cross-mapping figure using shadedErrorBar.m
%  

%% Changelog
% initial version written 20160519 by KHPD, based on Sugihara et al 2012
%

%%

fig = figure();
set (fig, 'Units', 'normalized', 'Position', [0,0,0.9,0.9], 'Color', [1,1,1]);
hold on
x=Lvals;
y1=XxmapY;
y2=YxmapX;
h1=shadedErrorBar(x,median(y1),[abs(prctile(y1,75)-median(y1)); abs(prctile(y1,25)-median(y1))],'g',1);
h2=shadedErrorBar(x,median(y2),[abs(prctile(y2,75)-median(y2)); abs(prctile(y2,25)-median(y2))],'b',1);
xlim([min(Lvals) max(Lvals)])
box on
title({sprintf('CCM for %s and %s',xName,yName),sprintf('(%s, E=%d, \tau = %d)',dataName,E,tau)})
xlabel('L')
ylabel('\rho')
legend([h1.mainLine, h2.mainLine],sprintf('%s xmap %s',xName,yName),sprintf('%s xmap %s',yName,xName), 'Location', 'best')
hold off

if ~isempty(saveName)
    print(sprintf('pngs/%s',saveName),'-dpng','-r600')
    savefig(sprintf('figs/%s.fig',saveName))
end

end