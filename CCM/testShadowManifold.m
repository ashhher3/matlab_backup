function [ Xtrue, Xest ] = testShadowManifold( Mx, E)
%% testShadowManifold.m - part of Matlab CCM code
% Takes E-dimensional "shadow manifold" (as built by 
% makeShadowManifold.m), splits it in two pieces, and uses each to estimate 
% the other

%% Changelog
% based on multispatial CCM R package by Adam Clark
% https://cran.r-project.org/web/packages/multispatialCCM/index.html
%
% initial version written 20160502 by KHPD

%%

start = E;
stop = size(Mx,1)-E+1;

Xtrue = Mx(start:stop,1);
Xest = NaN(size(Xtrue)); % can't calculate for actually every point

% for each point in Mx
for i = start:stop
    % remove nearest TIME neighbors
    
    
    % find nearest manifold neighbors
    [N1, D1] = knnsearch(testmfld, Mx(i,:),'k',E+2);
    N1 = N1(:,2:end); %nearest nbr is the point itself, so throw out 1st clm
    D1 = D1(:,2:end);
    
    % calculate and store estimate
    W = exp( -1*(D1(i,:)) / (D1(i,1)));
    W = W / sum(W);
    Xest(i) = W*chunkA(N1(i,:),1);
    

end



end

