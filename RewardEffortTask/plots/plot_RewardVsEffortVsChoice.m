%% blocks to points

Date = '20151210';
Session = '1003';

totalBlocks = DATA(end).BlockNum;
startIdx = numel(blockAvg); %0 for first session
trial=1;
for currBlock=1:totalBlocks
    choseUp=0;
    total=0;
    while DATA(trial).BlockNum<=currBlock
        if (DATA(trial).OutcomeID==0) || (DATA(trial).OutcomeID==6)
            total = 1+total;
            choseUp = choseUp + DATA(trial).TrialChoiceID;
        end
        UpR = DATA(trial).UpReward;
        DownR = DATA(trial).DownReward;
        UpE = DATA(trial).UpEffort;
        DownE = DATA(trial).DownEffort;
        trial = 1+trial;
    end
    idx = startIdx + currBlock;
    blockAvg(idx).Date = Date;
    blockAvg(idx).Session = Session;
    blockAvg(idx).BlockNum = currBlock;
    blockAvg(idx).UpR = UpR;
    blockAvg(idx).DownR = DownR;
    blockAvg(idx).UpE = UpE;
    blockAvg(idx).DownE = DownE;
    blockAvg(idx).Choice = choseUp/total;
end

%%
clear choseUp currBlock Date DownE DownR idx Session startIdx total totalBlocks trial UpE UpR

%% make new columns

for i=1:numel(blockAvg)
    blockAvg(i).RUpToTot = blockAvg(i).UpR/(blockAvg(i).UpR + blockAvg(i).DownR);
    blockAvg(i).EUpToTot = (1/blockAvg(i).UpE)/((1/blockAvg(i).UpE) + (1/blockAvg(i).DownE));
end
clear i

%% linear discriminator for R(up/total) vs E(up/total)
X=[[blockAvg(:).RUpToTot]',[blockAvg(:).EUpToTot]'];

Y=round([blockAvg(:).Choice]);

class = fitcdiscr(X,Y);

clear X Y

%% R(up/total) vs E(up/total)
figure(31)
hold on
scatter([blockAvg(:).RUpToTot],[blockAvg(:).EUpToTot],[],[blockAvg(:).Choice],'filled');

K = class.Coeffs(1,2).Const;
L = class.Coeffs(1,2).Linear;

% Plot the curve K + [x1,x2]*L  = 0:
f = @(x1,x2) K + L(1)*x1 + L(2)*x2;
h3 = ezplot(f,[0 1 0.25 0.75]);
h3.Color = 'k';
h3.LineWidth = 2;

xlabel( 'Reward index offered (up/total)' )
ylabel( 'Effort index offered (up/total)' )
title({'NHP choices for month of October' '(Averaged over each block)'})

%% R(up/total) vs E(up/down)
figure(32)
scatter([blockAvg(:).RUpToTot],[blockAvg(:).EUpToDown],[],[blockAvg(:).Choice],'filled');
xlabel( 'Up Reward / (Up Reward + Down Reward)' )
ylabel( 'Up Effort / Down Effort' )
title('#2 - effort is Up/Down')

%% R(up/total) vs E(up-down)
figure(33)
UminusD = (1./[blockAvg(:).UpE])-(1./[blockAvg(:).DownE]);
%l=log(UminusD);
scatter([blockAvg(:).RUpToTot],UminusD,[],[blockAvg(:).Choice],'filled');
xlabel( 'Up Reward / (Up Reward + Down Reward)' )
ylabel( 'Up Effort - Down Effort' )
title('#3 - effort is Up - Down')


%% R(up/total) vs E(log(up/down))
figure(34)
scatter([blockAvg(:).RUpToTot],[blockAvg(:).LogEUpToDown],[],[blockAvg(:).Choice],'filled');
xlabel( 'Up Reward / (Up Reward + Down Reward)' )
ylabel( 'log(Up Effort / Down Effort)' )
title('#4 - effort is log(Up/Down)')


%% binned version of figure 31

rmin = 0;
rmax = 1;
rstep = (rmax-rmin)/50;

emin = 0;
emax = 0.75;
estep = (emax - emin)/50;

raxis = rmin:rstep:rmax;
eaxis = emin:estep:emax;

num = zeros(length(raxis),length(eaxis));
denom = num;
e=nan(size(num));
r=e;

for i=1:length(blockAvg)
    ridx = fix((blockAvg(i).RUpToTot - rmin)/rstep)+1;
    eidx = fix((blockAvg(i).EUpToTot - emin)/estep)+1;
    num(ridx,eidx) = blockAvg(i).Choice + num(ridx,eidx);
    denom(ridx,eidx) = 1 + denom(ridx,eidx);
    r(ridx,eidx) = ridx*rstep;
    e(ridx,eidx) = eidx*estep;
end

denom(denom==0)=NaN;

grid = num ./ denom;

figure(35)
hold on
pcolor(raxis,eaxis,grid');
shading flat

K = class.Coeffs(1,2).Const;
L = class.Coeffs(1,2).Linear;

% Plot the curve K + [x1,x2]*L  = 0:
f = @(x1,x2) K + L(1)*x1 + L(2)*x2;
h3 = ezplot(f,[0 1 0.25 0.75]);
h3.Color = 'k';
h3.LineWidth = 2;


ylim([0.25 0.75])
xlabel( 'Up Reward / (Up Reward + Down Reward)' )
ylabel( 'Up Effort / (Up Effort + Down Effort)' )
title('#5 - binned version of #1')
