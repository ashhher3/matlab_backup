%% load data and trim

idx = [65, 96; 149, 158; 159, 168];
idx = idx-64;
ofc = []; am=[]; hp=[];

for f=1:4
load(sprintf('EC108_2fd7f1ac_all_%03d.mat',f))
ofc = [ofc, mean(data(idx(1,1):idx(1,2),:))];
am = [am, mean(data(idx(2,1):idx(2,2),:))];
hp = [hp, mean(data(idx(3,1):idx(3,2),:))];
end

last = 1000000;
ofc = ofc(1:last);
am = am(1:last);
hp = hp(1:last);

clear idx last

%% set params and initialize variables

E = 30;
tau = 30;
minL = (tau +1)*E + 2 - tau;
n = 100;
maxL = length(ofc)/n;

startIdx = 1:maxL:length(ofc);
Lvals = minL:2:maxL;

OxmapA = NaN(n, length(Lvals));
AxmapO = NaN(n, length(Lvals));
OxmapH = NaN(n, length(Lvals));
HxmapO = NaN(n, length(Lvals));
AxmapH = NaN(n, length(Lvals));
HxmapA = NaN(n, length(Lvals));

%% do CCM

for i = 1:length(Lvals)
    for j = 1:length(startIdx)
        
        L = Lvals(i);
        start = startIdx(j);
        
        fprintf('L = %d, iteration %d\n', L, j);
        
        % xmap OFC & Am
        X = ofc(start:start+L-1)';
        Y = am(start:start+L-1)';
        [ ~, Yest, Xest, ~, ~ ] = CCM( X, Y, E, tau );
        
        R = corrcoef(X(size(X,1)-size(Xest,1)+1:end), Xest);
        AxmapO(j,i) = R(1,2);
        
        R = corrcoef(Y(size(X,1)-size(Yest,1)+1:end), Yest);
        OxmapA(j,i) = R(1,2);
        
        % xmap OFC & Hp
        X = ofc(start:start+L-1)';
        Y = hp(start:start+L-1)';
        [ ~, Yest, Xest, ~, ~ ] = CCM( X, Y, E, tau );
        
        R = corrcoef(X(size(X,1)-size(Xest,1)+1:end), Xest);
        HxmapO(j,i) = R(1,2);
        
        R = corrcoef(Y(size(X,1)-size(Yest,1)+1:end), Yest);
        OxmapH(j,i) = R(1,2);
        
        % xmap Am & Hp
        X = hp(start:start+L-1)';
        Y = am(start:start+L-1)';
        [ ~, Yest, Xest, ~, ~ ] = CCM( X, Y, E, tau );
        
        R = corrcoef(X(size(X,1)-size(Xest,1)+1:end), Xest);
        AxmapH(j,i) = R(1,2);
        
        R = corrcoef(Y(size(X,1)-size(Yest,1)+1:end), Yest);
        HxmapA(j,i) = R(1,2);
         
        
    end
end

save('20160426_EC108_CCM.mat')