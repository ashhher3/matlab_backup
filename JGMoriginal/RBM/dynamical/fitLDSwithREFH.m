% fitLDSwithREFH
%   This script trains a series of "EFH filters" on data that encode a 
%   dynamical system.  The dynamics themselves are set in setDynamics.m.  
%   Then it tests the trained network, with the classic test (error 
%   covariances across all trials).
%
%   The point is to generate errors bars across tests of the model.  Cf.
%   the m-files that trains with EM, fitLDSwithEM.m.

%-------------------------------------------------------------------------%
% Cribbed: 01/22/15
%   -from LearnDynamics
%   by JGM
%-------------------------------------------------------------------------%


% train an RBM-filter on dynamical data
clear; clc; 


% initialize
params = setParams;
[~,machine] = system('hostname');
params.machine = strtrim(machine);
params.dynamics.T = 1000;
testData = getLDSdata(params);
TESTDECODING = 1;


% Ns
Ntest = 5;
hid2propRatio = 6;
Nxprmt = 12;
Nruns = 20;
Nsensory = params.Nmods*params.N^params.Ndims;
params.numsUnits = [Nsensory + hid2propRatio*Nsensory, hid2propRatio*Nsensory];
params.t = hid2propRatio*Nsensory;
params.DBNmaxepoch = 90;
numsUnits = params.numsUnits;
numRBMs = length(numsUnits)-1;
Nbatches = params.dynamics.T;
Ncases = params.Ncases;
%%%%%%%%%%%
params.DBNmaxepoch = 1200;
yrsigmoid = @(iEp)(1000/(1 + exp(-iEp/8 + 7.5)) + 0.5);
%%%%%%%%%%%


% loop through experiments (to gather stats over) and run (find one best)
datagenargs = {};
for iXprmt = 1:Nxprmt

    % malloc/reint
    errorStatVec = nan(Nruns,1,'like',testData.Z);
    for iRun = 1:Nruns

        % start over: display and reset wts
        paramDisplay(params);
        wts = cell(numRBMs*2,1);
        
        % train an EFH/DBN
        for i_rbm = 1:numRBMs
            Nvis = numsUnits(i_rbm);
            Nhid = numsUnits(i_rbm+1);
            fprintf(1,'Pretraining Layer %i w/RBM: %d-%d \n',i_rbm,Nvis,Nhid);
            restart = 1;
            
            % train
            tic; EFH; toc;
            
            % pack together weights for saving (hid => recog., vis => gener.)
            wts{i_rbm} = [vishid; hidbiases];
            wts{numRBMs*2-i_rbm+1} = [vishid'; visbiases'];
            
            % for next time through
            batchdata = batchposhidmeans; 
        end
        
        % store the error variance
        [~,~,pEFH] = EFHfilter(testData,wts,params);
        pEFH.name = 'rEFH';
        EFHstats = testDynamics(testData,params,0,pEFH);
        close all;
        errorStatVec(iRun) = det(EFHstats(strcmp(params.mods,params.NS)).Cvrn);
    
        % store these weights if they're the best
        if errorStatVec(iRun) == min(errorStatVec)
            bestwts = wts;
        end
       
    end
   
    Allwts{iXprmt} = bestwts;

end

%%%%%
save('results\finalwts\wts1DrEFHManyXprmts1502XX','params','Allwts');
%%%%%

