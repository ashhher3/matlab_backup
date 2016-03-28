function ac = autocorrelationForCircularVars(S,params)
% NB: this fxn assumes that Ndims = 1

% params
[Ntraj,T] = size(S);
N = params.N;
trajmin = params.smin;
trajmax = N/(N-1)*(params.smax - params.smin) + params.smin;
tauMax = T-1;

% autocorrelation of the data
acx = 0;
acy = 0;
for iTraj = 1:Ntraj
    s = scalefxn(squeeze(S(iTraj,:)),trajmin,trajmax,0,2*pi);
    acx = acx + xcorr(cos(s),cos(s),tauMax,'none');
    acy = acy + xcorr(sin(s),sin(s),tauMax,'none');
    
    %%%%
%     s = squeeze(S(iTraj,:));
%     acx = acx + xcorr(s,s,tauMax,'coeff');
%     acy = acy + xcorr(s,s,tauMax,'coeff');
    %%%%
end
ac = [acx; acy]/Ntraj;




end













