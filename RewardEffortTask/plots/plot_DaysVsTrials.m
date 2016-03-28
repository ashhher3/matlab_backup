figure(61)
hold on
scatter(22:63, [TrialsPerDay(22:end).Trials],'bo','filled')
ylabel('Number of successes')
ax=gca;
ticks=22:5:64;
set(ax,'XTick',ticks);
set(ax,'XTickLabel',{'Sept. 01','Sept. 08','Sept. 15','Sept. 22','Sept. 29','Oct. 07','Oct. 14','Oct.21','Oct. 28'})
linex=20:0.1:65;
liney=ones(size(linex))*1223.9;
plot(linex,liney,'r--');
xlim([20 65]);
title('Successful Trials per day in September & October 2015')