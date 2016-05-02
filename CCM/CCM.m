function [ Xest, Yest, Xself, Yself, Mx, My ] = CCM( X, Y, E, tau )
% CCM Summary of this function goes here
%   Detailed explanation goes here

% 20160502 - CHANGED ORDER OF OUTPUTS!!

% make the time-lagged manifolds
Mx = makeShadowManifold(X, E, tau);
My = makeShadowManifold(Y, E, tau);

% check how well data estimates itself with this E and this tau
Xself = crossMap(Mx, Mx, E);
Yself = crossMap(My, My, E);

% do cross-mapping
Yest = crossMap(Mx, My, E);
Xest = crossMap(My, Mx, E);

end

function Yest = crossMap(Mx, My, E)
Yest = NaN( size(My, 1), 1);

% get E+1 nearest neighbors for all points in Mx
[N, D] = knnsearch(Mx, Mx,'k',E+2);
N = N(:,2:end); %nearest nbr is the point itself, so throw out 1st clm
D = D(:,2:end);

% compute weights
W=NaN(1, size(D,2));
% seems like should be possible to do this without a loop but I can't
% figure out how to do the indices.....
for i = 1:size(Mx, 1)
    W = exp( -1*(D(i,:)) / (D(i,1)));
    W = W / sum(W);
    Yest(i) = W*My(N(i,:),1);
end


end

function Mx = makeShadowManifold( X , E, tau)
L = length(X);
tstart = 1+(E-1)*tau;
Mx = NaN( (L-tstart), E);

for t=tstart:L
    for j=0:E-1
        Mx(t-tstart+1,j+1) = X(t-j*tau);
    end
end

end
