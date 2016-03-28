function [fig,p,stats] = makeCountVCohPlot( count, coh)

idx=find([coh]>=0.455);
count=count(idx);
coh=coh(idx);

x=[ count', ones(size(count'))] ;
y=coh';

% p=polyfit(x,y,1);
% yfit=polyval(p,x);
% yresid=y-yfit;
% SSresid=sum(yresid.^2);
% SStotal=(length(y)-1)*var(y);
% rsq=1-SSresid/SStotal;

[p,~,~,~,stats]=regress(y,x);

xplot=min(count):(max(count)-min(count))/100000:max(count);
yplot=polyval(p,xplot);

fig=figure();
hold on
scatter(count,y,'bo')
plot(xplot,yplot,'r')
mTextBox = uicontrol('style','text');
set(mTextBox,'String',sprintf('y=%.2dx+%.2d\nR^2=%.2d, p=%.2d',p(1),p(2),stats(1),stats(3)))
set(mTextBox,'Units','characters')
set(mTextBox,'BackgroundColor',[1 1 1])
set(fig,'Color',[1 1 1])


end