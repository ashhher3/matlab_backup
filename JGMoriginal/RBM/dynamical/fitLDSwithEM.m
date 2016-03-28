tic

% params
Niter = 20;
Nruns = 12;
params = setParams;
params.dynamics.meta = 'RandInit';
%%% params.dynamics.meta = 'RandInitWithEC';
[~,machine] = system('hostname');
params.machine = strtrim(machine);
params.dynamics.T = 1000;

% malloc
allLLs = zeros(Niter,Nruns);

for iRun = 1:Nruns
    LLbest = -Inf;
    for iter = 1:Niter
        
        fprintf('iteration: %i\n',iter);
        %%% LDSparamsEM = EM4LDS(size(params.dynamics.A,2),params);
        %%% LDSparamsEM = EM4LDS(size(params.dynamics.G,1),params);
        LDSparamsEM = EM4LDS(3,params);
        allLLs(iter,iRun) = LDSparamsEM.LL;
        
        if LDSparamsEM.LL > LLbest
            LLbest = LDSparamsEM.LL;
            LDSparamsEMthebest = LDSparamsEM;
        end
    end

    Allparams(iRun) = LDSparamsEMthebest;
end

toc


%%% save(['LDSparamsEM',params.MODEL,params.dynamics.meta],'LDSparamsEM');