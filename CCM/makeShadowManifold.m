function Mx = makeShadowManifold( X , E, tau)
L = length(X);

tstart = 1+(E-1)*tau;

Mx = NaN( (L-tstart), 3);

for t=tstart:L
    Mx(t-tstart+1,:) = [ X(t) , X(t-tau) , X(t-2*tau) ];
end

end

