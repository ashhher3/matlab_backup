function pNAIVE = getNaiveFilter(pSENSORY,params)

% params
C = params.dynamics.C;
Ncases = params.Ncases;

% convert posterior covariance (back) into sufficient statistics
T0.Shat = shortdata(Ncases,4,(C*longdata(pSENSORY.Xhat)')');
T0.CvrnMat = pSENSORY.CvrnMat;

% pretend that you don't know anything about the repelling walls
params.dynamics.walls = 'other';

% just give it all the parameters (but these are naive)
T = size(T0.Shat,4);
KFparams.T = T;
KFparams.A = params.dynamics.A;
KFparams.C = params.dynamics.C;
KFparams.SigmaX = params.dynamics.SigmaX;
KFparams.mu0 = [params.dynamics.muX0; params.dynamics.muV0];
KFparams.Info0 = setInfoMatrix(params.dynamics.SigmaX0,...
    params.dynamics.SigmaV0,params.Ndims);

% run the Kalman filter
%%%%%%%%%%%%%%%%%%%%%% this is now teh wrong syntax
pNAIVE = KF4PPC(T0,KFparams,params);
%%%%%%%%%%%%%%%%%%%%%%

% label
setColors;
pNAIVE.color = EYEcolor;        pNAIVE.name = 'naive';


end