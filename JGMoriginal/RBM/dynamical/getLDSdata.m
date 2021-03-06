function LDSdata = getLDSdata(params,varargin)
% getLDSdata    Get data from an LTI system
%
%   USAGE:
%       LDSdata = getLDSdata(params,varargin)
%
% Generate some new training data for the KF: 
%
%   LDSdata.Z:              stats (not wrapped)
%   LDSdata.U:              control??
%   LDSdata.Y:              state observation
%   LDSdata.SigmaY          Cvrn[Y|X]
%   LDSdata.R:              PPCs
%   LDSdata.S:              stimuli (wrapped!!)
%   LDSdata.restarts{1:T}:  restart vector

%-------------------------------------------------------------------------%
% Revised: 08/26/14
%   -eliminated "noisily observed control"
% Revised: 01/10/14
%   -restored readout of Cvrn as well as Shat from PPCs
% Revised: 01/06/14
%   -changed output to structure to accommodate a controlled LDS
% Revised: 12/19/13
%   -added varargout for Udata
% Revised: 12/18/13
%   -put X into "canonical structure" (a dimension for Nmods)
% Cribbed: 10/30/13
%   -from getLDSparams
% Created: 10/21/13
%   by JGM
%-------------------------------------------------------------------------%

% Ns
t = params.t;
Ncases = params.Ncases;
Ndims = params.Ndims;
mods = params.mods;
T = params.dynamics.T;

% get the data in PPC form
if isempty(varargin)
    [D0,S,~,Q] = DATAGENPP(T,params);
else
    S = varargin{1};
    [D0,S,~,Q] = DATAGENPP(size(S,1)/params.Ncases,params,'stimuli',S);
end
filterdata = Q.filterdata; clear Q;
switch params.typeUnits{1}
    case {'BP','GP','BG','Bernoulli'}
        PPCs = D0(:,(t+1):end,:);
    case 'PB'
        PPCs = D0(:,1:(end-t),:);
    otherwise
        error('strange type of units for this function -- jgm');
end

% collect the states and observations into useful form
Z = cat(3,filterdata(:).states);
if strcmp(params.dynamics.walls,'controlled')
    EC = S(:,:,strcmp(mods,'Efference-Copy'),:);
    Z = cat(2,X,reshape(EC,[Ncases,Ndims*size(EC,3),T]));
end
[Y,SigmaY] = PPC2FltrCmlnts(PPCs,S,Z,params);


% collect outputs
LDSdata.Z = Z;                          % state and control, *not* wrapped
LDSdata.Y = Y;                          % sufficient stats for 
LDSdata.SigmaY = SigmaY;                %   the state: T_z(R^z)
LDSdata.R = PPCs;                       % population responses
LDSdata.S = S;                          % the encoded stimuli, wrapped!

% restart data
[LDSdata.restarts{1:T}] = filterdata(:).RESTART;
%-------------------------------------------------------------------------%

end
