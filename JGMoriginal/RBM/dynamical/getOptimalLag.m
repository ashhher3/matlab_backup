function [taus,maxMI,Ptau] = getOptimalLag(V,S,params)
% getOptimalLag     Using mutual information
% USAGE:
%   [taus,maxMI,Ptau] = getOptimalLag(V,S,params)
%
% Find the "optimal" lag of the stimuli in S for Bernoulli probabilities V.
% That is, V is a tensor (Ntraj x Nunits x T) of conditional probabilities,
% p(V=1|r), for Bernoulli random variables, and S is a tensor (Ntraj x T) 
% of the corresponding stimuli.  This function sample-averages out R under 
% p(r|s) to get p(V=1|s), but it only reports the result, Ptau, at the lag,
% taus, that maximizes the mutual information between S and V.  This MI is
% also reported as maxMI.
%
% The usual params structure is also required.

%-------------------------------------------------------------------------%
% Created: 07/20/14
%   by JGM
%-------------------------------------------------------------------------%


% params
[Ntraj,Nunits,T] = size(V);
N = params.N;
trajmin = params.smin;
trajmax = N/(N-1)*(params.smax - params.smin) + params.smin;
Nsbins = 30;
% tauMax = 40;
tauMax = 80;


% malloc
P = zeros(tauMax,Nsbins);
MI = zeros(tauMax,1);
maxMI = zeros(Nunits,1);
taus = zeros(Nunits,1);
Ptau = zeros(Nunits,Nsbins);

tic;
fprintf('looping through hidden units...\n');
for iUnit = 1:Nunits;
    
    for t0 = 1:tauMax    
        Nsamples = (T-t0+1)*Ntraj;
        
    
        % marginal distribution over V
        v = squeeze(V(:,iUnit,t0:T));
        pofV = [1-sum(v(:))/Nsamples; sum(v(:))/Nsamples];

        % marginal distribution over S
        s = S(:,1:(T-t0+1));
        [Ns,sbins] = histc(s(:),linspace(trajmin,trajmax,Nsbins+1));
        pofS = Ns(1:Nsbins)/Nsamples;
        
        % joint distribution
        Vincidence = [1-v(:), v(:)];
        Sincidence = sparse(1:Nsamples,sbins,1,Nsamples,Nsbins);
        pofSandV = Vincidence'*Sincidence/Nsamples;
        
        % mutual information
        pofSandV(pofSandV==0) = pofSandV(pofSandV==0) + eps;
        MI(t0) = sum(sum(pofSandV.*(log2(pofSandV) - log2(pofV*pofS'))));
        
        % Pr(v|s) for this lag
        pofVgivenS = bsxfun(@mrdivide,pofSandV,pofS');
        P(t0,:) = pofVgivenS(2,:); %%% just V=1
    end
    fprintf('.');
    
    % keep the maximum MI, pofVgivenS, and the lag
    [maxMI(iUnit),indMax] = max(MI);
    taus(iUnit) = indMax-1;
    Ptau(iUnit,:) = P(indMax,:);
    
    
end
fprintf('\n');
toc

% find lag that maximizes mutual information
% [maxMI,indMax] = max(MI,[],2);
% taus = indMax-1;


end