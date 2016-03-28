%% add column for max forces

n=numel(success);

for i=1:n
    xTrace = success(i).ForceTrace(:,4);
    yTrace = success(i).ForceTrace(:,5);
    
   % transform to same units as target
    xTrace = 300*xTrace*(50/10) + 0;
    yTrace = 300*yTrace*(50/10) + 0;
    
    switch success(i).TrialChoiceID
        case 1 %chose up
            ymax = max(yTrace);
            xmax = max(abs(xTrace));
            
        case 0 %chose down
            ymax = abs(min(yTrace));
            xmax = max(abs(xTrace));
    end
    
    success(i).MaxX = xmax;
    success(i).MaxY = ymax;
end

clear i n xTrace yTrace xmax ymax

%%
maxY = NaN(max([idx(:).end] -[idx(:).start]),numel(idx));
maxX = maxY;

for d = 1:numel(idx) % for each day
    % pull out maxes and store
    mx = [success(idx(d).start:idx(d).end).MaxX];
    my = [success(idx(d).start:idx(d).end).MaxY];
    
    maxX(1:numel(mx),d) = mx;
    maxY(1:numel(my),d) = my;
end

clear mx my d

%%
EU = zeros(size(maxY));
ED = EU;
UU = EU;
UD = EU;
DU = EU;
DD = EU;

for d = 1:numel(idx) % for each day
    s=idx(d).start;
    e=idx(d).end;
    
    EU(1:e-s+1,d)=[[success(s:e).TrialType]==3] & [[success(s:e).TrialChoiceID]==1];
    ED(1:e-s+1,d)=[success(s:e).TrialType]==3 & [success(s:e).TrialChoiceID]==0;
    
    UU(1:e-s+1,d)=[[success(s:e).UpEffort]<[success(s:e).DownEffort]] & [[success(s:e).TrialChoiceID]==1];
    UD(1:e-s+1,d)=[[success(s:e).UpEffort]<[success(s:e).DownEffort]] & [[success(s:e).TrialChoiceID]==0];
    DU(1:e-s+1,d)=[[success(s:e).UpEffort]>[success(s:e).DownEffort]] & [[success(s:e).TrialChoiceID]==1];
    DD(1:e-s+1,d)=[[success(s:e).UpEffort]>[success(s:e).DownEffort]] & [[success(s:e).TrialChoiceID]==0];
end

EU([EU==0])=NaN;
ED([ED==0])=NaN;
UU([UU==0])=NaN;
UD([UD==0])=NaN;
DU([DU==0])=NaN;
DD([DD==0])=NaN;

clear s e d

%%
figure(41)

subplot(2,3,1)
hold on
boxplot(maxY.*EU/30)
xlabel('Days')
ylabel('Peak force')
title('Vertical force; Equal effort; Chose Up')
ylim([0 9])

subplot(2,3,4)
boxplot(maxY.*ED/30)
xlabel('Days')
ylabel('Peak force')
title('Vertical force; Equal effort; Chose Down')
ylim([0 9])

subplot(2,3,2)
boxplot(maxY.*UU/30)
xlabel('Days')
ylabel('Peak force')
title('Vertical force; Up Hard; Chose Up')
ylim([0 9])

subplot(2,3,5)
boxplot(maxY.*UD/30)
xlabel('Days')
ylabel('Peak force')
title('Vertical force; Up Hard; Chose Down')
ylim([0 9])

subplot(2,3,3)
boxplot(maxY.*DU/30)
xlabel('Days')
ylabel('Peak force')
title('Vertical force; Down Hard; Chose Up')
ylim([0 9])

subplot(2,3,6)
boxplot(maxY.*DD/30)
xlabel('Days')
ylabel('Peak force')
title('Vertical force; Down Hard; Chose Down')
ylim([0 9])

%%
figure(51)

subplot(2,3,1)
boxplot(maxX.*EU/30)
xlabel('Days')
ylabel('Peak force')
title('Horizontal force; Equal effort; Chose Up')
ylim([0 9])

subplot(2,3,4)
boxplot(maxX.*ED/30)
xlabel('Days')
ylabel('Peak force')
title('Horizontal force; Equal effort; Chose Down')
ylim([0 9])

subplot(2,3,2)
boxplot(maxX.*UU/30)
xlabel('Days')
ylabel('Peak force')
title('Horizontal force; Up Hard; Chose Up')
ylim([0 9])

subplot(2,3,5)
boxplot(maxX.*UD/30)
xlabel('Days')
ylabel('Peak force')
title('Horizontal force; Up Hard; Chose Down')
ylim([0 9])

subplot(2,3,3)
boxplot(maxX.*DU/30)
xlabel('Days')
ylabel('Peak force')
title('Horizontal force; Down Hard; Chose Up')
ylim([0 9])

subplot(2,3,6)
boxplot(maxX.*DD/30)
xlabel('Days')
ylabel('Peak force')
title('Horizontal force; Down Hard; Chose Down')
ylim([0 9])
