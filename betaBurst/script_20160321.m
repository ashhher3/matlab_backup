%% beta power averages

%acc
acc.betaPwrAvg=[];
for j=1:size(acc.signal,1)
    acc.betaPwrAvg(j,:)=getPwr10Sec(acc.betaPwr(j,:),fs);
end

%am
am.betaPwrAvg=[];
for j=1:size(am.signal,1)
    am.betaPwrAvg(j,:)=getPwr10Sec(am.betaPwr(j,:),fs);
end

%hp
hp.betaPwrAvg=[];
for j=1:size(hp.signal,1)
    hp.betaPwrAvg(j,:)=getPwr10Sec(hp.betaPwr(j,:),fs);
end

%ofc
ofc.betaPwrAvg=[];
for j=1:size(ofc.signal,1)
    ofc.betaPwrAvg(j,:)=getPwr10Sec(ofc.betaPwr(j,:),fs);
end

%% MEAN COHERENCE
% ACC-AM, ACC-HP
for i=1:10
    for j=1:10
        coh.acc(i).am(j)=mean(getCoh10Sec(acc.beta(i,:), am.beta(j,:), fs));
        coh.acc(i).hp(j)=mean(getCoh10Sec(acc.beta(i,:), hp.beta(j,:), fs));
    end
end

% ACC - OFC

for i=1:10
    for j=1:32
        coh.acc(i).ofc(j)=mean(getCoh10Sec(acc.beta(i,:), ofc.beta(j,:), fs));
    end
end

% AM-HP
for i=1:10
    for j=1:10
        coh.am(i).hp(j)=mean(getCoh10Sec(am.beta(i,:), hp.beta(j,:), fs));
    end
end


% AM- OFC
for i=1:10
    for j=1:32
        coh.am(i).ofc(j)=mean(getCoh10Sec(am.beta(i,:), ofc.beta(j,:), fs));
    end
end

% HP - OFC
for i=1:10
    for j=1:32
        coh.hp(i).ofc(j)=mean(getCoh10Sec(hp.beta(i,:), ofc.beta(j,:), fs));
    end
end

%% COHERENCE FOR SPECIFIC CHANNELS

coh.acc_chan = 4;
coh.am_chan = 10;
coh.hp_chan = 10;
coh.ofc_chan = 19;

coh.AccAm = getCoh10Sec(acc.beta(coh.acc_chan,:),am.beta(coh.am_chan,:),fs);
coh.AccHp = getCoh10Sec(acc.beta(coh.acc_chan,:),hp.beta(coh.hp_chan,:),fs);
coh.AccOFC = getCoh10Sec(acc.beta(coh.acc_chan,:),ofc.beta(coh.ofc_chan,:),fs);
coh.AmHp = getCoh10Sec(am.beta(coh.am_chan,:),hp.beta(coh.hp_chan,:),fs);
coh.AmOFC = getCoh10Sec(am.beta(coh.am_chan,:),ofc.beta(coh.ofc_chan,:),fs);
coh.HpOFC = getCoh10Sec(hp.beta(coh.hp_chan,:),ofc.beta(coh.ofc_chan,:),fs);

%% put all the things into one matrix
vars = {'ACC beta power', 'Am beta power', 'Hp beta power', ...
    'OFC beta power', 'Acc-Am coherence', 'Acc-Hp coherence', ...
    'Acc-OFC coherence', 'Am-Hp coherence', 'Am-OFC coherence', ...
    'Hp-OFC coherence'};

x = NaN(size(vars,2), size(coh.AccAm,2));

x(1,:) = acc.betaPwrAvg(coh.acc_chan,:);
x(2,:) = am.betaPwrAvg(coh.am_chan,:);
x(3,:) = hp.betaPwrAvg(coh.hp_chan,:);
x(4,:) = ofc.betaPwrAvg(coh.ofc_chan,:);

x(5,:) = coh.AccAm;
x(6,:) = coh.AccHp;
x(7,:) = coh.AccOFC;
x(8,:) = coh.AmHp;
x(9,:) = coh.AmOFC;
x(10,:) = coh.HpOFC;

%% regress all the things

n = size(x,1);
m = size(x,2);
R = zeros(n,n);
P = NaN(n,n);
L = NaN(n,n,2);

for i=1:n
    for j=1:n
        if j~=i
            [l,~,~,~,stats]=regress(x(j,:)',[x(i,:)', ones(m,1)]);
            
            R(i,j)=stats(1);
            P(i,j)=stats(3);
            L(i,j,1:2)=l;
            
        end
    end
end

clear n m l stats i j

%% make some plots

Rpad = [ flipud(R), zeros( size(R,1), 1)];
Rpad = [ Rpad ;
    zeros(1, size(Rpad, 2))];


Ppad = [ flipud(P), zeros( size(P,1), 1)];
Ppad = [ Ppad ;
    zeros(1, size(Ppad, 2))];

f1 = figure();
pcolor(Rpad)
ax1=gca;
shading flat
caxis([0 1])
colorbar
f1.Color=[1 1 1];
title(sprintf('R-squared values for linear regression\n (EC108, values averaged over 10s windows)'))
ax1.YTick = 1.5:1:10.5;
ax1.YTickLabel = fliplr(vars);
ax1.XTickLabel = {};

f2 = figure();
pcolor(-1*log10(Ppad))
ax2=gca;
shading flat
colorbar
f2.Color=[1 1 1];
title(sprintf('p values for linear regression\n (EC108, values averaged over 10s windows)'))
% c = get(colorbar,'YTick');
% colorbar('Ytick',c,'YTickLabel',-1*10.^c);
c = [1e-60 1e-50 1e-40 1e-30 1e-20 1e-10 0.005 1];
c=fliplr(c);
caxis(-1*log10([c(1) c(end)]));
colorbar('FontSize',11,'YTick',-1*log10(c),'YTickLabel',c);
ax2.YTick = 1.5:1:10.5;
ax2.YTickLabel = fliplr(vars);
ax2.XTickLabel = {};

%% clean up workspace
clear ax ax1 ax2 f f1 f2 L2 P2 R2 Rpad Ppad c

%% make some linear regression plots

%% ACC beta vs Hp beta
a = 1;
b = 3;

indep = x(a,:)';
dep = x(b,:)';
l=L(a,b,:);
l=l(:);
xplot=min(indep):(max(indep)-min(indep))/100000:max(indep);
yplot=polyval(l,xplot);

fig=figure();
hold on
c=linspace(0,10*length(indep),length(indep));
scatter(indep,dep,[],c)
plot(xplot,yplot,'b')
text((max(indep)-min(indep))/2,(max(dep)-min(dep))/4, sprintf('y=%.2dx+%.2d\nR^2=%.2d, p=%.2d',l(1),l(2),R(a,b),P(a,b)))
set(fig,'Color',[1 1 1])
xlabel(vars(a))
ylabel(vars(b))

clear a b dep indep l xplot yplot

%% ACC-Am coh vs ACC-Hp coh
a = 5;
b = 6;

indep = x(a,:)';
dep = x(b,:)';
l=L(a,b,:);
l=l(:);
xplot=min(indep):(max(indep)-min(indep))/100000:max(indep);
yplot=polyval(l,xplot);

fig=figure();
hold on
c=linspace(0,10*length(indep),length(indep));
scatter(indep,dep,[],c)
plot(xplot,yplot,'b')
text((max(indep)-min(indep))/2,(max(dep)-min(dep))/4, sprintf('y=%.2dx+%.2d\nR^2=%.2d, p=%.2d',l(1),l(2),R(a,b),P(a,b)))
set(fig,'Color',[1 1 1])
xlabel(vars(a))
ylabel(vars(b))

clear a b dep indep l xplot yplot

%% OFC-Am coh vs OFC-Hp coh
a = 9;
b = 10;

indep = x(a,:)';
dep = x(b,:)';
l=L(a,b,:);
l=l(:);
xplot=min(indep):(max(indep)-min(indep))/100000:max(indep);
yplot=polyval(l,xplot);

fig=figure();
hold on
c=linspace(0,10*length(indep),length(indep));
scatter(indep,dep,[],c)
plot(xplot,yplot,'b')
text((max(indep)-min(indep))/2,(max(dep)-min(dep))/4, sprintf('y=%.2dx+%.2d\nR^2=%.2d, p=%.2d',l(1),l(2),R(a,b),P(a,b)))
set(fig,'Color',[1 1 1])
xlabel(vars(a))
ylabel(vars(b))

clear a b dep indep l xplot yplot

%% comparing coh to ACC and coh to OFC at high or low Am-Hp coh times
a = 7; 

% b = 5; % acc - am
% c = 9; % ofc - am

b = 6; % acc - hp
c = 10; % ofc - hp

indep1 = [];
dep1 = [];
indep2 = [];
dep2 = [];

thresh = median(x(a,:));

for i=1:length(x(a,:))
    if x(a,i)>=thresh
        indep1 = [indep1; x(b,i)];
        dep1 = [dep1; x(c,i)];
    else
        indep2 = [indep2; x(b,i)];
        dep2 = [dep2; x(c,i)];
    end
end

fig = figure();
hold on

scatter(indep1, dep1, 'co')
scatter(indep2, dep2, 'mo')

[l,~,~,~,stats]=regress(dep1,[indep1, ones(size(indep1))]);
xplot=min(indep1):(max(indep1)-min(indep1))/100000:max(indep1);
yplot=polyval(l,xplot);
plot(xplot,yplot,'b')
text(0.405,0.6 , sprintf('High %s\ny=%.2dx+%.2d\nR^2=%.2d, p=%.2d',char(vars(a)),l(1),l(2),stats(1),stats(3)), 'Color', 'blue')

[l,~,~,~,stats]=regress(dep2,[indep2, ones(size(indep2))]);
xplot=min(indep2):(max(indep2)-min(indep2))/100000:max(indep2);
yplot=polyval(l,xplot);
plot(xplot,yplot,'r')
text(0.6, 0.44, sprintf('Low %s\ny=%.2dx+%.2d\nR^2=%.2d, p=%.2d',char(vars(a)),l(1),l(2),stats(1),stats(3)), 'Color', 'red')

xlabel(vars(b))
ylabel(vars(c))
fig.Color=[1 1 1];


clear a b dep indep l xplot yplot


%% clean workspace

clear ans fig c












