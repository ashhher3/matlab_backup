function mixtureoferlangs
% mixtureoferlangs
%
% Try to reproduce the two distribution of amygdala-hippocampus coherences
% on slide 6 of Lowry Kirkby's presentation (labmeetingLong_13nov15.pdf) as
% mixtures of two Erlang distributions.
%
% As it turns out, you can match both quite nicely using a single set of
% mixing proportions, and a single "low-coherence state."  Only the high-
% coherence states differ: the "sad state" has a smaller shape parameter
% (waiting for fewer arrivals) but a higher slightly larger mean (k*mu;
% i.e., arrivals come less frequently).
%
% To double check, you also compute analytically the variance of the
% mixtures.  They match closely to the reported results in the figures on
% slide 9.
%
% Parameters were found by hand.

%-------------------------------------------------------------------------%
% Created: 12/09/15 (happy b'day, L.N.)
%   by JGM
%-------------------------------------------------------------------------%


alp = 1;
%%% multiplying the scale by 2 "doubles" the x axis and "halves" the y axis
dx = 0.001;
xf = 650*dx*alp;
x = 0:dx:xf;
figure(2); clf; hold on;


% low-cohererence state
p = 5/8;                                % fraction of time in this state
k1 = 1;                                 % shape parameter
mu1 = 0.05*alp;                         % scale parameter

% high-coherence state A
k2 = 16;
mu2 = 0.01*alp;
y = p*gampdf(x,k1,mu1) + (1-p)*gampdf(x,k2,mu2);
plot(x,y,'k');

vrnc = getErlangMixtureVariance([k1;k2],[mu1;mu2],[p;1-p]);
fprintf('first Erlang-mixture variance is %2.4f\n',vrnc);

% high-coherence state B
k3 = 7;
mu3 = 0.03*alp;
y = p*gampdf(x,k1,mu1) + (1-p)*gampdf(x,k3,mu3);
plot(x,y,'r');

axis([0,xf,0,14/alp]);

vrnc = getErlangMixtureVariance([k1;k3],[mu1;mu3],[p;1-p]);
fprintf('second Erlang-mixture variance is %2.4f\n',vrnc);


Nsamples = 10000;
Y = generateErlangMixtures([k1 k1; k2 k3],[mu1 mu1; mu2 mu3],...
    [p, p; 1-p, 1-p],Nsamples);


% y = (10/16)*gampdf(x,1,0.05) + (3/16)*gampdf(x,16,0.01) + (3/16)*gampdf(x,40,0.01);
% plot(x,y,'r'); axis([0,xf,0,14/alp]);


end
%-------------------------------------------------------------------------%
%-------------------------------------------------------------------------%



%-------------------------------------------------------------------------%
function [xpct,m2_1] = getErlangMoments(ks,mus)
% 1st and 2nd moment of Erlang distribution

xpct = ks.*mus;
vrnc = ks.*mus.^2;
m2_1 = xpct.^2 + vrnc;

end
%-------------------------------------------------------------------------%


%-------------------------------------------------------------------------%
function vrnc = getErlangMixtureVariance(ks,mus,ps)

[m1s,m2s] = getErlangMoments(ks,mus); % 1st and 2nd moment of Erlang
vrnc = ps'*m2s - (ps'*m1s)^2;

end
%-------------------------------------------------------------------------%


%-------------------------------------------------------------------------%
function Y = generateErlangMixtures(K,Mu,P,Nsamples)

% Ns
Ndstrbs = size(P,2);
Nbins = 20;

% prepare colors
clrNames(:,1) = mat2cell(repmat('gray  ,opacity=0.5',Nbins,1),ones(Nbins,1),18);
clrNames(:,2) = mat2cell(repmat('red   ,opacity=0.5',Nbins,1),ones(Nbins,1),18);
clrNames(:,3) = mat2cell(repmat('yellow,opacity=0.5',Nbins,1),ones(Nbins,1),18);
clrNames(:,4) = mat2cell(repmat('blue  ,opacity=0.5',Nbins,1),ones(Nbins,1),18);


% malloc
Ns = zeros(Ndstrbs,Nbins);
Y = zeros(Ndstrbs,Nsamples);

% draw samples
X = linspace(0,0.6,Nbins);
for iDstrb = 1:Ndstrbs
    Y(iDstrb,:) = ErlangMixtureRnd(K(:,iDstrb),Mu(:,iDstrb),P(:,iDstrb),Nsamples);
    Ns(iDstrb,:) = hist(Y(iDstrb,:),X);
end
Ns = bsxfun(@rdivide,Ns,sum(Ns,2))';

tikzBarGraph(Ns,zeros(Nbins,2,Ndstrbs),4,0.32,...
    mat2cell(num2str(X','%2.2f'),ones(Nbins,1),4),...
    'correlation coefficient','fraction','',clrNames(:,1:Ndstrbs),...
    ['erlangsHist',date]);

end
%-------------------------------------------------------------------------%

















