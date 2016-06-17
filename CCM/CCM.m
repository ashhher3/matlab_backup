function [ XxmapY, YxmapX, Lvals, startIdx ] = CCM( X, Y, E, tau, n, maxL, Lskip )
%% CCM.m Tests for causality between data series X and Y using convergent cross-mapping
%   Inputs: X, Y - data series of same length that may be causally related.
%                  Should be column vectors, but I will flip them for you.
%           E    - embedding dimension, ie, dimension of underlying
%                  dynamics, use testE.m to help choose
%           tau  - time lag to use when making shadow manifolds, measured
%                  in number of samples, NOT SECONDS
%           n    - number of segments to test and average over
%           maxL - maximum length of data to test, OPTIONAL, default is
%                  length(X)/n
%           Lskip- step between values of L to test, OPTIONAL, default is 1
%
%   Outputs: XxmapY - matrix of corrcoef values for predicting Y from X
%                     each row is a diff starting point in orginal data
%                     each column is a different value of L
%           YxmapX  - same as XxmapY but predicting X from Y
%           Lvals   - vector of L values used, measured in samples
%           startIdx - indices of original data X and Y used as starting
%                    points

%% Changelog
% initial version written 20160420 by KHPD, based on Sugihara et al 2012
%
% 20160502 - CHANGED ORDER OF OUTPUTS!!
% 20160505 - removed self-mapping part (redundant with
%               testShadowManifold.m)
%            removed makeShadowManifold as subfunction since it is a
%               separate .m file now
% 20160506 - reworked to actually do convergence test 
%            changed inputs/outputs, incompatible w earlier scripts
%            removed crossMap as subfunction
% 20160519 - moved print statement to outer loop
%            

%%

% do some input checking
if isrow(X)
    X = X';
end
if isrow(Y)
    Y = Y';
end
if ~(iscolumn(X) && iscolumn(Y))
    error('X and Y must be 1-dimensional data series.')
end
if ~(length(X) == length(Y))
    warning('X and Y are different lengths, are you sure about this?\n')
end
if E<=1
    error('E must be >=2')
end
if tau<1
    error('tau is measured in samples and must be >=1')
end
if isempty(n) && isempty(maxL)
    error('You must provide either n or maxL')
end
if isempty(maxL)
    maxL = floor(length(X)/n);
end
if isempty(n)
    n = floor(length(X)/maxL);
end
if n*maxL > length(X)
    error('To test n=%d segments of length maxL=%d, need %d samples, you gave %d.',n, maxL, n*maxL, length(X))
end
if isempty(Lskip)
    Lskip = 1;
end

% initialize variables
minL = (tau +1)*E + 2 - tau;
Lvals = [minL:Lskip:maxL, maxL]; % test maxL even if Lskip doesn't hit it
startIdx = 1:maxL:length(X);
startIdx = startIdx(1:n);
XxmapY = NaN(n, length(Lvals));
YxmapX = NaN(n, length(Lvals));

% loop over start points
for j = 1:length(startIdx)
    % set start idx
    start = startIdx(j);
    
    fprintf('starting point %d\n', j); 
    
    % loop over Lvals backwards
    for i = length(Lvals):-1:1
        % set L
        L = Lvals(i);
      
        % if L is max L, make biggest shadow manifold
        if L==maxL
            Mx = makeShadowManifold( X(start:start+L-1), E, tau);
            My = makeShadowManifold( Y(start:start+L-1), E, tau);
        else % trim off Lskip many points from the end
            Mx = Mx(1:end-Lskip,:);
            My = My(1:end-Lskip,:);
        end
        
        % do cross mapping
        Xest = crossMap(My, Mx, E);
        Yest = crossMap(Mx, My, E);
        
        % compute corrcoef for Xtrue & Xest, and store in YxmapX
        R = corrcoef(Mx(1:end,1), Xest);
        YxmapX(j,i) = R(1,2);
        
        % compute corrcoef for Xest & Xtrue, and store in XxmapY
        R = corrcoef(My(1:end,1), Yest);
        XxmapY(j,i) = R(1,2);
        
    end
end
    


end