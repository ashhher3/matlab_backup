function params = setParams
% SETPARAMS     Sets parameters for sensory integration learning
%   SETPARAMS sets the (reused) parameters of DBN, DATAGEN, FK2link, etc.
%-------------------------------------------------------------------------%
% Revised: 12/09/13
%   -rationalized: put in case statements for all the different models
% Revised: 09/07/12
%   -fixed max and min calculations
% Revised: 12/20/10
%   -changed to accomodate 3-modality scheme and non-2D-stimuli.
% Revised: 7/14/10
%   -changed the grid size---while keeping the relative tuning acuity
%       constant.
% Revised: 7/7/10
%   -fixed posmax (should be measured from origin, not from left edge)
% Created: 7/5/10
%   by JGM
%-------------------------------------------------------------------------%


% model
% params.MODEL = '1DrEFHbern';
params.MODEL = '1Dinteg';


% PATH
path(path,'../ee125/fxns');
path('../utils',path);      % contains tex.m
% path('../utils/matlab2tikz',path);
path(path,'../tools')
% path(path,'retired');
% path(path,'scratch');
% path(path,'tuningcurves');
path(path,'results');
path(path,'parallel');
path(path,'dynamical');
% path(path,'../PPC models');



% network architecture
params.nettype = 'ENCODER';                 % /'CLASSIFIER'/??
% numsUnits = [500 500 2000];                 % GEH CLASSIFIER

% training params                           % [GEH params??]
params.Ncases = 40;                         % [10 cases]
params.DBNmaxepoch = 90;                    % [50 epochs]
params.mw = 500;  % 500;                    % [50 newtons]
params.mvb = 120; % 120;                    % [50 newtons]
%params.mw2 = 500;  % 500;                   % [50 newtons]
params.mhb = 120; % 120;                    % [50 newtons]
params.b = params.mw/2;                     % [250? N-s/m]
params.k = 0.001;                           % [0.02 N/m]
params.Ts = 1;                              % "sampling interval"
params.amass = 1.10;                        % A "matrix" for masses
params.numCDsteps = 1;

% backprop params [unused]
params.BPmaxepoch = 5;
params.max_iter = 3;                        % number of linesearches
params.numMinibatches = 10;                 % to be combined in a big batch



% which model?
switch params.MODEL
    
    case '2Dinteg'
        
        % data generation
        params.Ndims = 2;                   % encoded vars/neuron
        params.Nmods = 2;                   % number of modalities
        params.N = 30;                      % number of neurons *** (2) ***
        params.g = 15;                      % gain
        params.swing = 0.2;                 % swing in gain (max=1)
        params.NS = 'Joint-Angle';          % "neutral space"
       
        % arm properties
        params.thmin = [-pi/2; pi/4];   % [-pi/4; 0];
        params.thmax = [pi/4; 3*pi/4];  % [pi/4; pi/2];
        params.L1 = 12;                     % upper arm
        params.L2 = 20;                     % forearm
        
        params = getLimits(params);
        params.smin = [params.posmin params.thmin];
        params.smax = [params.posmax params.thmax];
        
        % modality names
        params.mods{1} = 'Hand-Position';
        params.mods{2} = 'Joint-Angle';
        
        % RBM units
        Nvis = params.Nmods*params.N^params.Ndims;
        params.numsUnits = [Nvis, Nvis/2];
        params.typeUnits = {'Poisson','Bernoulli'};
        
        % learning rates/momentum
        params.mw =  500;                   % [50 newtons]
        params.mvb = 120;                   % [50 newtons]
        params.mhb = 120;                   % [50 newtons]
        params.b = params.mw/2;             % [250? N-s/m]


    case '2Dtwoarms'
        
        % data generation
        params.Ndims = 2;                   % encoded vars/neuron
        params.Nmods = 2;                   % number of modalities
        params.N = 30;                      % number of neurons *** (2) ***
        params.g = 15;                      % gain
        params.swing = 0.2;                 % swing in gain (max=1)
        params.NS = 'Joint-Angle-Left';    	% "neutral space"
       
        % arm properties
        params.smin = [-pi/2 -pi/2; pi/4 pi/4];
        params.smax = [pi/4 pi/4; 3*pi/4 3*pi/4];  % [pi/4; pi/2];
                
        % modality names
        params.mods{1} = 'Joint-Angle-Left';
        params.mods{2} = 'Joint-Angle-Right';
        
        % RBM units
        Nvis = params.Nmods*params.N^params.Ndims;
        params.numsUnits = [Nvis, Nvis/2];
        params.typeUnits = {'Poisson','Bernoulli'};
        
        % learning rates/momentum
        params.mw =  500;                   % [50 newtons]
        params.mvb = 120;                   % [50 newtons]
        params.mhb = 120;                   % [50 newtons]
        params.b = params.mw/2;             % [250? N-s/m]

        
    case '1Daddition'
        
        % data generation
        params.Ndims = 1;                   % encoded vars/neuron
        params.Nmods = 3;                   % number of modalities
        params.N = 60;                      % number of neurons *** (2) ***
        params.g = 15;                      % gain
        params.gains = [5 20 5]; % [5 15 5];        % overwrite gains %%%%%
        params.swing = 0.2;                 % swing in gain (max=1)
        params.NS = 'Hand-Position';        % "neutral space"
        
        % arm properties
        params.thmin = pi/6;    % pi/8; % pi/4;
        params.thmax = 5*pi/6;  % 7*pi/8; % 3*pi/4;
        params.L1 = 12;
        params = getLimits(params);
        params.smin = [params.posmin params.thmin params.eyemin];
        params.smax = [params.posmax params.thmax params.eyemax];
        
        % modality names
        params.mods{1} = 'Hand-Position';
        params.mods{2} = 'Joint-Angle';
        params.mods{3} = 'Gaze-Angle';
        
        % RBM units
        Nvis = params.Nmods*params.N^params.Ndims;
        params.numsUnits = [Nvis, 160];
        params.typeUnits = {'Poisson','Bernoulli'};
        
        
    case '1Dinteg'
        
        % data generation
        params.Ndims = 1;                   % encoded vars/neuron
        params.Nmods = 2;                   % number of modalities
        params.N = 60;                      % number of neurons *** (2) ***
        %%% params.N = 120;                      % number of neurons *** (2) ***
        params.g = 15;                      % gain
        % params.g = 0.5; % 0.08;
        params.swing = 0.2;                 % swing in gain (max=1)
        % params.swing = 0;
        params.NS = 'Joint-Angle';
        
        % arm properties
        params.thmin = pi/6;    % pi/8; % pi/4;
        params.thmax = 5*pi/6;  % 7*pi/8; % 3*pi/4;
        params.L1 = 12;
        params = getLimits(params);
        params.smin = [params.posmin params.thmin];
        params.smax = [params.posmax params.thmax];
        
        % modality names
        params.mods{1} = 'Hand-Position';
        params.mods{2} = 'Joint-Angle';
        
        % RBM units
        Nvis = params.Nmods*params.N^params.Ndims;
        params.numsUnits = [Nvis, Nvis/2];
        params.typeUnits = {'Poisson','Bernoulli'};
        %params.typeUnits = {'Bernoulli','Bernoulli'};
        %params.typeUnits = {'Binomial','Binomial'};
        
        %params.nexperiments = 200; %for binomial
                
        
    case '2Daddition'
        
        % data generation
        params.Ndims = 2;                   % encoded vars/neuron
        params.Nmods = 3;                   % number of modalities
        params.N = 30;                      % number of neurons *** (2) ***
        params.g = 15;                      % gain
        params.swing = 0.2;                 % swing in gain (max=1)
        params.NS = 'Hand-Position';        % "neutral space"
        
        % arm properties
        params.thmin = [-pi/2; pi/4];   % [-pi/4; 0];
        params.thmax = [pi/4; 3*pi/4];  % [pi/4; pi/2];
        params.L1 = 12;                     % upper arm
        params.L2 = 20;                     % forearm
        params = getLimits(params);
        params.smin = [params.posmin params.thmin params.eyemin];
        params.smax = [params.posmax params.thmax params.eyemax];
        
        % modality names
        params.mods{1} = 'Hand-Position';
        params.mods{2} = 'Joint-Angle';
        params.mods{3} = 'Gaze-Angle';
        
        % RBM units
        Nvis = params.Nmods*params.N^params.Ndims;
        params.numsUnits = [Nvis, Nvis];
        params.typeUnits = {'Poisson','Bernoulli'};
        
        
        
        
    case '2DrEFH'
        
        % data generation
        params.Ndims = 2;                   % encoded vars/neuron
        params.Nmods = 1;                   % number of modalities
        params.N = 15;                      % number of neurons *** (2) ***
        params.g = 10; % 6;                       % gain
        params.swing = 0.2;                 % swing in gain (max=1)
        params.NS = 'Joint-Angle';
        
        % arm properties
        params.thmin = [-3*pi/8; -pi/4]; % [-pi/2; pi/4]; % shoulder, elbow
        params.thmax = [+3*pi/8; +pi/4]; % [pi/4; 3*pi/4];% shoulder, elbow
        params.L1 = 12;                     % upper arm
        params.L2 = 20;                     % forearm
        params.smin = params.thmin;
        params.smax = params.thmax;
        
        % modality names
        params.mods = {'Joint-Angle'};
        
        % dynamical params
        params.dynamics = setDynamics(params);
        
        % RBM units
        Nvis = params.Nmods*params.N^params.Ndims;
        params.numsUnits = [Nvis + 10*Nvis, 10*Nvis];
        params.typeUnits = {'BP','Bernoulli'};
        params.t = params.numsUnits(2);
        
        
    case '2DrEFHwithEC'
        
        % data generation
        params.Ndims = 2;                   % encoded vars/neuron
        params.Nmods = 2;                   % number of modalities
        params.N = 15;                      % number of neurons *** (2) ***
        params.g = 10; % 6;                       % gain
        params.swing = 0.2;                 % swing in gain (max=1)
        params.NS = 'Joint-Angle';
        
        % arm properties
        params.thmin = [-3*pi/8; -pi/4]; % [-pi/2; pi/4];       % shoulder, elbow
        params.thmax = [+3*pi/8; +pi/4]; % [pi/4; 3*pi/4];      % shoulder, elbow
        params.L1 = 12;                     % upper arm
        params.L2 = 20;                     % forearm
        
        % EC properties
        params.ctrlmin = -[0.2; 0.2];       % -[0.06; 0.06];
        params.ctrlmax = [0.2; 0.2];        % [0.06; 0.06]; % *** (5) ***
        params.smin = [params.thmin, params.ctrlmin];
        params.smax = [params.thmax, params.ctrlmax];
        
        
        % modality names
        params.mods = {'Joint-Angle','Efference-Copy'};
        
        
        % dynamical params
        params.dynamics = setDynamics(params);
        
        % RBM units
        Nvis = params.Nmods*params.N^params.Ndims;
        params.numsUnits = [Nvis + 12*Nvis, 12*Nvis];
        params.typeUnits = {'BP','Bernoulli'};
        params.t = params.numsUnits(2);
        
        
    case '1DrEFH'
        
        % data generation
        params.Ndims = 1;                   % encoded vars/neuron
        params.Nmods = 1;                   % number of modalities
        params.N = 15;                      % number of neurons *** (2) ***
        params.g = 8;                       % gain
        params.swing = 0.2;                 % swing in gain (max=1)
        params.NS = 'Joint-Angle';
        
        % arm properties
        params.thmin = -pi/3; % pi/6;    % 0;
        params.thmax = +pi/3; % 5*pi/6;  % 2*pi;
        params.L1 = 12;
        params.smin = params.thmin;
        params.smax = params.thmax;
        
        % modality names
        params.mods = {'Joint-Angle'};
        
        % dynamical params
        params.dynamics = setDynamics(params);
        
        % RBM units
        Nvis = params.Nmods*params.N^params.Ndims;
        params.numsUnits = [Nvis + 16*Nvis, 16*Nvis];
        params.typeUnits = {'BP','Bernoulli'};
        params.t = params.numsUnits(2);
        
        
        
    case 'HVNdamped'
        
        % data generation
        params.Ndims = 1;                   % encoded vars/neuron
        params.Nmods = 1;                   % number of modalities
        params.N = 15;                      % number of neurons *** (2) ***
        params.g = 10;                      % gain
        params.swing = 0.2;                 % swing in gain (max=1)
        params.NS = 'Joint-Angle';
        
        % arm properties
        params.thmin = -pi/3;
        params.thmax = +pi/3;
        params.L1 = 12;
        params.smin = params.thmin;
        params.smax = params.thmax;
        
        % modality names
        params.mods = {'Joint-Angle'};
        
        % dynamical params 
        k=0; b=0.25; m=5; dt=0.05;
        %%%%%% slower!
        b = b*4;
        %%%%%%
        params.dynamics.A = [1.0000, dt; -k/m*dt, -(b/m*dt-1)];
        params.dynamics.C = [1,0];
        params.dynamics.SigmaX = 40*[5e-7, 0; 0, 5e-5];
        params.dynamics.muX0 = 0;
        params.dynamics.SigmaX0 = Inf;
        params.dynamics.SigmaV0 = 1e-4;
        params.dynamics.muV0 = 0;
        params.dynamics.walls = 'wrapping';
        
        % RBM units
        Nvis = params.Nmods*params.N^params.Ndims;
        params.numsUnits = [Nvis + 15*Nvis, 15*Nvis];
        params.typeUnits = {'BP','Bernoulli'};
        params.t = params.numsUnits(2);
        
        
    case '1DrEFHwithEC'
        
        % data generation
        params.Ndims = 1;                   % encoded vars/neuron
        params.Nmods = 2;                   % number of modalities
        params.N = 15;                      % number of neurons *** (2) ***
        params.g = 8;                       % gain
        params.swing = 0.2;                 % swing in gain (max=1)
        params.NS = 'Joint-Angle';
        
        % arm properties
        params.thmin = -pi/3; % pi/6;
        params.thmax = +pi/3; % 5*pi/6;
        params.L1 = 12;
        params.smin = params.thmin;
        params.smax = params.thmax;
        
        % EC properties
        params.ctrlmin = -1.25;
        params.ctrlmax = 1.25;
        %%% you can solve for these in closed form by computing the
        %%% variance of a drift/diffusion process, -> 3*stddev
        params.smin = [params.thmin, params.ctrlmin];
        params.smax = [params.thmax, params.ctrlmax];
        
        % modality names
        params.mods = {'Joint-Angle','Efference-Copy'};
        
        % dynamical params
        params.dynamics = setDynamics(params);
        
        % RBM units
        Nvis = params.Nmods*params.N^params.Ndims;
        params.numsUnits = [Nvis + 8*Nvis, 8*Nvis];
        params.typeUnits = {'BP','Bernoulli'};
        params.t = params.numsUnits(2);
        
        % learning rates
        params.mw = 500;
        params.mvb = 120;
        params.mhb = 120;
        params.b = params.mw/2;
        
        
    case '3DrEFH'
        
        % data generation
        params.Ndims = 3;                   % encoded vars/neuron
        params.Nmods = 1;                   % number of modalities
        params.N = 15;                      % number of neurons *** (2) ***
        params.g = 10; % 6;                       % gain
        params.swing = 0.2;                 % swing in gain (max=1)
        params.NS = 'Joint-Angle';
        
        % arm properties
        params.thmin = [0; 0; 0];
        params.thmax = [2*pi; 2*pi; 2*pi];
        params.smin = params.thmin;
        params.smax = params.thmax;
        
        % modality names
        params.mods = {'Joint-Angle'};
        
        % dynamical params
        params.dynamics = setDynamics(params);
        
        % RBM units
        Nvis = params.Nmods*params.N^params.Ndims;
        params.numsUnits = [Nvis + 1*Nvis, 1*Nvis];
        params.typeUnits = {'BP','Bernoulli'};
        params.t = params.numsUnits(2);
        
        
        
        
        
    case 'Binomial' % or whatever
        %%%%%%%%%%%%%%%%%
        params.nexperiments = 200;          % for binomial neurons
        %%%%%%%%%%%%%%%%%
        % case
        
        %%%%%%%%%%%%%%%%%%%%%%%
        % params.gains = [22 8 8];
        % params.gains = [10 30 10];
        % params.gains = [2 5 2];
        %%%%%%%%%%%%%%%%%%%%%%%
        
    case 'MCDexperiment'
        
        % just for this case
        dmin = 40;      % mm
        dmax = 115;     % mm
        Fmin = 5; % 100;     % Hz
        Fmax = 15; % 300;     % Hz
        Nelectrodes = 8;
        PDs = 2*pi*(0:Nelectrodes-1)/Nelectrodes - pi;
        
        params.experiment.Fmin = Fmin;
        params.experiment.Fmax = Fmax;
        params.experiment.dmin = dmin;
        params.experiment.dmax = dmax;
        params.experiment.Nelectrodes = Nelectrodes;
        params.experiment.PDs = PDs;
        
        electrodeFunc = @(th,d)(...
            Fmax/dmax*bsxfun(@times,d,(1+cos(bsxfun(@minus,th,PDs)))/2));
        
        params.experiment.electrodeFunc = electrodeFunc;
        params.experiment.neuronsPerElectrode = 30;
        
        
        % data generation
        params.Ndims = 2;                   % encoded vars/neuron
        params.Nmods = 2;                   % number of modalities
        params.N = 15;                      % 24^2 = 8*72
        params.g = 15;                      % gain
        params.swing = 0;                   % swing in gain (max=1)
        params.NS = 'Motion-Dots';          % "neutral space"
        
        % arm properties
        params.smin = [-dmax; -dmax];
        params.smax = [dmax; dmax];
        
        % modality names
        params.mods = {'Motion-Dots','ICMS'};
        
        % RBM units
        Nvis = params.N^params.Ndims +...
            Nelectrodes*params.experiment.neuronsPerElectrode;
        params.numsUnits = [Nvis, floor(Nvis/2)];
        params.typeUnits = {'Poisson','Bernoulli'};
        
        
        
    case 'MCDdarpa'
        
        % distance min and max
        dmin = 40;      % mm
        dmax = 115;     % mm
        params.experiment.dmin = dmin;
        params.experiment.dmax = dmax;
        
        
        % data generation
        params.Ndims = 2;                   % encoded vars/neuron
        params.Nmods = 2;                   % number of modalities
        params.N = 15;                      % number of neurons *** (2) ***
        params.g = 15;                      % gain
        params.swing = 1;                   % swing in gain (max=1)
        params.NS = 'Motion-Dots';
        
        % arm properties
        params.smin = [-dmax -dmax; -dmax -dmax];
        params.smax = [dmax dmax; dmax dmax];
        
        % modality names
        params.mods = {'Motion-Dots','ICMS'};
        
        % RBM units
        Nvis = params.Nmods*params.N^params.Ndims;
        params.numsUnits = [Nvis, Nvis/2];
        params.typeUnits = {'Poisson','Bernoulli'};
        
        
        
    case 'BMMdata'
        path(path,'../SUNY');
        
        % data generation
        params.Ndims = 2;                   % encoded vars/neuron
        params.Nmods = 2;                   % number of modalities
        params.N = 30;                      % number of neurons *** (2) ***
        params.g = 15;                      % gain
        params.swing = 0.2;                 % swing in gain (max=1)
        params.NS = 'Joint-Angle';          % "neutral space"
        
        % modality names
        params.mods{1} = 'Hand-Position';
        params.mods{2} = 'Joint-Angle';
        
        % RBM units
        Nvis = params.Nmods*params.N^params.Ndims;
        params.numsUnits = [Nvis, Nvis/2];
        params.typeUnits = {'Poisson','Bernoulli'};
        
        
        
        
        
    case '1DrEFHbern'
        
        % data generation
        params.Ndims = 1;                   % encoded vars/neuron
        params.Nmods = 1;                   % number of modalities
        params.N = 60;                      % number of neurons *** (2) ***
        params.g = 0.08;                    % gain
        %%%params.swing = 0.2;                 % swing in gain (max=1)
        params.swing = 0;
        params.NS = 'Joint-Angle';
        
        % arm properties
        params.thmin = -pi/3; % pi/6;    % 0;
        params.thmax = +pi/3; % 5*pi/6;  % 2*pi;
        params.L1 = 12;
        params.smin = params.thmin;
        params.smax = params.thmax;
        
        % modality names
        params.mods = {'Joint-Angle'};
        
        % dynamical params
        Nwindows = 100;
        k=3; b=0.25; m=5; dt=0.05/Nwindows;
        params.dynamics.A = [1.0000, dt; -k/m*dt, -(b/m*dt-1)];
        params.dynamics.C = [1,0];
        params.dynamics.SigmaX = [5e-7, 0; 0, 5e-5]/Nwindows;
        params.dynamics.muX0 = 0;
        params.dynamics.SigmaX0 = Inf;
        params.dynamics.SigmaV0 = 5e-10;
        params.dynamics.muV0 = 0;
        params.dynamics.walls = 'wrapping';
        
        % RBM units
        Nvis = params.Nmods*params.N^params.Ndims;
        params.numsUnits = [Nvis + 10*Nvis, 10*Nvis];
        params.typeUnits = {'Bernoulli','Bernoulli'};
        params.t = params.numsUnits(2);
        
        
    otherwise
        
        error('unrecognized model -- jgm');
        
end


% universal PPC params
FWHM = 1/6;                                 % deg (Deneve01: "45-75deg")
c = FWHM/(2*sqrt(2*log(2)));                % std. dev. in degrees
params.C = eye(params.Ndims)*c^2;
%%% essentially, one std of the tuning curve (in one direction only) is
%%% about 7% of the total grid; so the 1-std (in both directions) coverage
%%% of each tuning curve is 14% of the total area.

% grid parameters
respLength = 1;                             % normalized to 1  *** (3) ***
margin = 4*c;                               % 99.99%           *** (4) ***
if isfield(params,'dynamics')
    if strcmp(params.dynamics.walls,'wrapping'), margin = 0; end
end
params.respLength = respLength;
params.margin = margin;
params.gridsize = margin + respLength + margin;
params.granularity = params.N/params.gridsize;


% setColors
params.VIScolor = [1 0 1];                  % magenta
params.PROPcolor = [1 0.5 0];               % orange
params.EYEcolor = [0 0 1];                  % blue
params.EFHcolor = [0 1 0];                  % green
params.OPTcolor = [0 0 0];                  % black


% say
fprintf('Using %f units...\n',params.granularity);

% check
if (length(params.numsUnits) ~= length(params.typeUnits))
    error('mismatch b/n number and types of layers --- jgm');
end


params.hier = 0;


end




% *** (1) ***
% the reachable workspace is hard to state in general, so I've hard-coded
% in the relevant params---see drawing in your notes.  Ought to be changed

% *** (2) ***
% cf. Ma et al.: "p(s|r) converges to a Gaussian as the number of neurons
%  increases."

% *** (3) ***
% Deneve01 has her neurons cover 45-75 degrees---so use 60 degrees.  Now,
% I think the space she's working with spans 360 degrees, so the FWHM for
% each of her neurons is 1/6 (60/360) the total space.  You want to keep
% the same *relative* precision, so your neurons will have FWHM, not 60,
% but 1/6 the total space---which you have normalized to unity.
% ---But in fact [06/06/11], in the 2D case this doesn't even matter! since
% the tuning-curve widths are irrelevant there.

% *** (4) ***
% The margin should include at least 99.99% of the Gaussian bubbles; so for
%  bubbles at the edge of the respLength, there needs to be a 4 stddev
%  margin.

% *** (5) ***
% The control limits were were determined empirically: 0.06 is
%   approximately 3 stddevs out from the mean control (0) that was used to
%   generate Lissajous curves *with no limit on the magnitude of the
%   control*.  The limits were subsequently put into
%   generateLissajousControls.m.

% *** (6) ***
% Rationale for the IC priors:
% Total number of trajectories (and therefore ICs):
%   DBNmaxepoch*Ncases = 90*40 = 3600
% Want E[number of rejections] < 1
% From lab notes:
%   (1) v0_rare = sqrt(umax/dt/A*(A^2 - (x0_rare - b)^2)).
% If x0_rare  = x0_{>2SD} and v0_rare = v0_{>3SD} then
%   (2) E[#rej] = (1-0.997)*(1-0.954)*3600 = 0.5 < 1 (ok)
% And x0 needs to be kept far away from the edges--say 5 SDs.
% Then:
%   (3) x0minusb_rare = 4*(xmax - xmin)/10; % (2SD each way)
%   (4) v0_rare = sqrt(0.06/dt./A.*(A.^2 - x0minusb_rare.^2))
% And 1 SD of v0 is:
%   (5) v0_{1SD} = v0_rare/3 = [0.5447, 0.4448]
% See assembleLissajousTensor for more details




% params.posmin = [-(cos(pi/4)*L1 + L2); sin(pi/4)*(L1 - L2)];
% params.posmax = [cos(pi/4)*(L1+L2); L1+L2]; % *** see (1) ***





%-------------------------------------------------------------------------%
% to update old params structures:
%
% params.smax = params.xmax;
% params.smin = params.xmin;
% params.Ndims = params.m;
% params.Nmods = params.r;
%
% params = rmfield(params,'xmax');
% params = rmfield(params,'xmin');
% params = rmfield(params,'r');
% params = rmfield(params,'m');
%
% params.MODEL = ...
% params.NS = ...
%-------------------------------------------------------------------------%



%-------------------------------------------------------------------------%
% code cleanup, July 2014:
% (1) Check/rehabilitate the following (not a complete list):
%       standardCond2.m
%       gainmods.m
%       shiftTest.m
%       [yr adaptation files]
%
% (2) Try to eliminate the following:
%       ntrlJacobian.m
%       estError.m
%       estGather.m
%       estStatsCorePP.m
%       decode.m
%       decoder.m
%       marginalErrors.m
%       marginalErrorsPP.m
%       PPCinfo.m
%       dispErrCovs.m
%       covInCalc.m
%       PPCinputStats.m
%       SICE.m
%       tuner.m [appears not to be in use]
%
% (3) See if the following are still useful:
%       testStim.m
%       condInf.m
%       condInfPP.m
%       condInfCondErrs.m
%       confabulate.m
%       getNoiseFloorPP.m
%       optMeanCheck.m
%       rbmPP.m
%
% (N) back up all the "final wts" to zamfir and mushroom
%-------------------------------------------------------------------------%



