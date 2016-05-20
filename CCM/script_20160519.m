%% KHPD 20160519

%% load EC100 data

clear all

load('/sabes_data3/kderosier/CCM/EC100_forCCM/excerpt_from_oneStruct_EC100_1f1a51b1_001.mat')

cd('/newhome/kderosier/Documents/MATLAB/CCM/')

%% take mean of each area to be CCM input data

ofc = mean(data.signal(1:8,:));
am = mean(data.signal(9:12,:));
hp = mean(data.signal(13:22,:));
acc = mean(data.signal(23:31,:));

%% check that they data looks like what i think it does

f=figure();
hold on
plot(ofc(1:512))
plot(am(1:512))
plot(hp(1:512))
plot(acc(1:512))
hold off
pause()
close(f)
clear f

%% set params for basic CCM

E = 10;
tau = 10;
n=100;
maxL = 400;
Lskip = 5;

%% run CCM for all pairs

% note Lvals, startIdx depend on length & params, which are same for all

[ OxmapA, AxmapO, ~, ~ ] = CCM( ofc, am, E, tau, n, maxL, Lskip );

[ OxmapH, HxmapO, ~, ~ ] = CCM( ofc, hp, E, tau, n, maxL, Lskip );

[ OxmapC, CxmapO, ~, ~ ] = CCM( ofc, acc, E, tau, n, maxL, Lskip );

[ CxmapA, AxmapC, ~, ~ ] = CCM( acc, am, E, tau, n, maxL, Lskip );

[ CxmapH, HxmapC, ~, ~ ] = CCM( acc, hp, E, tau, n, maxL, Lskip );

[ AxmapH, HxmapA, Lvals, startIdx ] = CCM( am, hp, E, tau, n, maxL, Lskip );


%% save

save('data/EC100_100contig_E10_tau10_maxL400.mat')


%%

dataName='EC100';

f1 = plotCCM(Lvals,OxmapA,AxmapO,E,tau,dataName,'OFC','Am','EC100_100contig_E10_tau10_OFCAm');

f2 = plotCCM(Lvals,OxmapH,HxmapO,E,tau,dataName,'OFC','Hp','EC100_100contig_E10_tau10_OFCHp');

f3 = plotCCM(Lvals,OxmapC,CxmapO,E,tau,dataName,'OFC','ACC','EC100_100contig_E10_tau10_OFCACC');

f4 = plotCCM(Lvals,CxmapA,AxmapC,E,tau,dataName,'ACC','Am','EC100_100contig_E10_tau10_ACCAm');

f5 = plotCCM(Lvals,CxmapH,HxmapC,E,tau,dataName,'ACC','Hp','EC100_100contig_E10_tau10_ACCHp');

f6 = plotCCM(Lvals,AxmapH,HxmapA,E,tau,dataName,'Am','Hp','EC100_100contig_E10_tau10_AmHp');

clear dataName

%% just for fun test against random noise

rdm = 500*rand(size(acc));

[ OxmapR, RxmapO, ~, ~ ] = CCM( ofc, rdm, E, tau, n, maxL, Lskip );


dataName='EC100';

f1 = plotCCM(Lvals,OxmapR,RxmapO,E,tau,dataName,'OFC','Random noise','');

clear dataName

%% clear stuff from above
clear OxmapR RxmapO rdm


