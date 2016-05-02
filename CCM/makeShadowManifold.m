function Mx = makeShadowManifold( X , E, tau)
%% makeShadowManifold.m - part of Matlab CCM code
% Takes vector X and constructs "shadow manifold" Mx using time-lagged
% versions of X. E is dimension and tau is time lag, measured in samples.
% Mx(i,:) = [ X(i), X(i-tau), ... , X(t-(E-1)*tau) ]

%% Changelog
% initial version written 20160420 by KHPD, based on Sugihara et al 2012
% added description and changelog 20160502 KHPD

%%

L = length(X);
tstart = 1+(E-1)*tau;
Mx = NaN( (L-tstart), E);

for t=tstart:L % is it possible to eliminate a loop by being clever?
    for j=0:E-1
        Mx(t-tstart+1,j+1) = X(t-j*tau);
    end
end

end
