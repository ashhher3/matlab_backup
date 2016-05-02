%%
% script for doing CCM on subnets data, original version 20160426 KHPD
% 20160502 adjusted code for change in order of CCM.m outputs - KHPD

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

%% set params and initialize variables

E = 10;
tau = 30;
minL = (tau +1)*E + 2 - tau;
n = 100;
maxL = length(ofc)/n;

startIdx = 1:maxL:length(ofc);
Lvals = [minL:100:maxL, maxL];

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
    save('20160428_EC108_CCM_E10.mat')
end

matlabmail('kderosier6@gmail.com', 'CCM is done', 'CCM is done')