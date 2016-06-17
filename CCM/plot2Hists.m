function [fig, h1, h2] = plot2Hists( data1, data2, binWidth, data1Name, data2Name,xName,yName, figTitle,saveName )

fig = figure();
set (fig, 'Units', 'normalized', 'Position', [0,0,0.9,0.9], 'Color', [1,1,1]);
hold on

h1=histogram(data1, 'Normalization','probability','BinWidth',binWidth);
h2=histogram(data2, 'Normalization','probability','BinWidth',binWidth);

box on
title(figTitle)
xlabel(xName)
ylabel(yName)
legend({data1Name, data2Name})
hold off


if ~isempty(saveName)
    print(sprintf('pngs/%s',saveName),'-dpng','-r600')
    savefig(sprintf('figs/%s.fig',saveName))
end

end

