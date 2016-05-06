function [ Xest, Yest, Mx, My ] = CCM_old( X, Y, E, tau )
% CCM Summary of this function goes here
%   Detailed explanation goes here

% 20160502 - CHANGED ORDER OF OUTPUTS!!
% 20160505 - removed self-mapping part (redundant with
%               testShadowManifold.m)
%            removed makeShadowManifold as subfunction since it is a
%               separate .m file now
% 20160506 - LEGACY VERSION OF CCM.M FOR SCRIPTS FROM BEFORE THIS DATE

% make the time-lagged manifolds
Mx = makeShadowManifold(X, E, tau);
My = makeShadowManifold(Y, E, tau);

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
