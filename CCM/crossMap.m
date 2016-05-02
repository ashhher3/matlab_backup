function Yest = crossMap(Mx, My, E)
%% crossMap.m - part of Matlab CCM code
% Takes two E-dimensional "shadow manifolds" (as built by 
% makeShadowManifold.m) and uses cross-mapping to estimate points on My
% from points on Mx

%% Changelog
% based on Sugihara et al 2012
% pulled out from CCM.m and turned into own file 20160502 KHPD

%%

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
