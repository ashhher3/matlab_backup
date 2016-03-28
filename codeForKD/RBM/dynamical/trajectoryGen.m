function [src,filterdata] = trajectoryGen(T,smin,smax,params)
% [states,RESTART] = TRAJECTORYGEN(Nbatches,smin,smax,params)
% Generate "trajectories" through source space.  Batch 1 contains the first
% sample of Ncases trajectories; batch 2 contains the second sample of
% those same Ncases trajectories; and so on.

%-------------------------------------------------------------------------%
% Revised: 09/03/14 (JGM)
%   -added cases for wrapping + control
% Revised: 12/18/13 (JGM)
%   -added case "controlled"
% Revised: 11/20/13 (BKD)
%   -changed "resetting" condition to use wall_reset.m
%   -fixed bug that filterdata.states wasn't getting updated after resets
% Revised: 07/01/13 (JGM)
%   -renamed "flag" to "RESTART," since "flag" is protected
%   -changed "RESTART" from a huge tensor to a list of restart indices
%   -put the "new" outputs, "RESTART" and and "states," into a structure,
%   "filterdata"
% Revised: 06/27/13 (BKD)
%   -WARNING: OUTPUTS WERE CHANGED HERE
%   -added 'restart' condition
%   -now outputs states, which contains velocity
% Revised: 05/11/13
%   -wrote a function to generate ICs from either a normal or a uniform
%   prior, according to your conventions
% Created: 05/07/13
%   by JGM
%-------------------------------------------------------------------------%



% params
NSIND = strcmp(params.mods,params.NS);
trajmin = smin(:,NSIND);
trajmax = smax(:,NSIND);
Ncases = params.Ncases;
Ndims = params.Ndims;
Nstates = size(params.dynamics.A,1);
range = trajmax - trajmin;
walls = params.dynamics.walls;


% dynamics
A = params.dynamics.A;
SigmaX = params.dynamics.SigmaX;
muX0 = params.dynamics.muX0;
SigmaX0 = params.dynamics.SigmaX0;
muV0 = params.dynamics.muV0;
SigmaV0 = params.dynamics.SigmaV0;
if isfield(params.dynamics,'muX'),
    muX = params.dynamics.muX; 
else
    muX = zeros(Nstates,1); 
end

% use a GPU?
if strcmp(params.machine,'domestica')
    A = gpuArray(A);        
    muX = gpuArray(muX);            SigmaX = gpuArray(SigmaX);  
    muX0 = gpuArray(muX0);          SigmaX0 = gpuArray(SigmaX0);
    muV0 = gpuArray(muV0);          SigmaV0 = gpuArray(SigmaV0);
    trajmin = gpuArray(trajmin);    trajmax = gpuArray(trajmax);
end


% initial conditions
mrgn = 0.05;
x0 = sampleStatePrior(muX0,SigmaX0,Ncases,trajmin,trajmax,mrgn);
v0 = sampleStatePrior(muV0,SigmaV0,Ncases,trajmin,trajmax,mrgn);

% init
Z = zeros(Ncases,size(A,1),T,'like',A);
Z(:,:,1) = [x0 v0];
filterdata(1).states = Z(:,:,1);
filterdata(1).RESTART = 1:Ncases;


% what happens at the walls?
switch walls
    
    case 'wrapping'
        
        switch params.MODEL
            case {'1DrEFH','2DrEFH','1DtRBM','3DrEFH','HVNdamped'}
                
                % loop through time
                for t = 1:(T-1)
                    
                    % propagate all trajectories forward one step
                    Z(:,:,t+1) = (A*Z(:,:,t)' +...
                        bsxfun(@plus,mvnrnd(zeros(Ncases,2*Ndims),SigmaX)',muX))';
                    filterdata(t+1).states = Z(:,:,t+1);
                end
                
                % wrap the position elements
                src = state2stim(params.dynamics.C,Z,trajmin,trajmax,...
                    params.N,Ncases);
                
            case {'1DrEFHwithEC','2DrEFHwithEC'}
                
                % control params
                F = params.dynamics.F;
                G = params.dynamics.G;
                H = params.dynamics.H;
                muU0 = params.dynamics.muU0;
                SigmaU0 = params.dynamics.SigmaU0;
                SigmaU = params.dynamics.SigmaU;
                umin = smin(:,strcmp(params.mods,'Efference-Copy'));
                umax = smax(:,strcmp(params.mods,'Efference-Copy'));
                
                % initialize controls
                U = zeros(Ncases,size(G,2),T);
                U(:,:,1) = sampleStatePrior(muU0,SigmaU0,Ncases,...
                    umin,umax,mrgn);
                filterdata(1).states = cat(2,Z(:,:,1),U(:,:,1));
                
                % forward simulate through time
                for t = 1:(T-1)
                    Z(:,:,t+1) = (A*Z(:,:,t)' + G*U(:,:,t)' + bsxfun(...
                        @plus,mvnrnd(zeros(Ncases,2*Ndims),SigmaX)',muX))';
                    U(:,:,t+1) = (F*U(:,:,t)' + bsxfun(@plus,...
                        mvnrnd(zeros(Ncases,size(G,2)),SigmaU)',0))';
                    filterdata(t+1).states = cat(2,Z(:,:,t+1),U(:,:,t+1));
                end
                
                % store output and input for encoding in PPCs
                src(:,:,strcmp(params.mods,params.NS),:) = state2stim(...
                    params.dynamics.C,Z,trajmin,trajmax,params.N,Ncases);
                src(:,:,strcmp(params.mods,'Efference-Copy'),:) =...
                    state2stim(H,U,umin,umax,params.N,Ncases);
                
                % for debugging purposes
                TOPLOT = 0;
                if TOPLOT
                    th0 = get2DOutline(params.thmin,params.thmax,42);
                    for iCase = 1:Ncases
                        
                        % plot
                        figure(1); clf; hold on;
                        
                        plot(squeeze(R(iCase,1,:)),squeeze(R(iCase,2,:)),'r');
                        plot(squeeze(Ybar(iCase,1,:)),squeeze(Ybar(iCase,2,:)),'m');
                        plot(th0(:,1),th0(:,2),'k');
                        
                        axis equal;
                        title(num2str(iCase));
                        hold off;
                        pause()
                    end
                end
        
                
                
            case {'1DrEFHbern'}
                
                % re-malloc for one very long trajectory
                Z = zeros(size(A,1),T*Ncases,'like',A);
                Z(:,1) = [x0(1) v0(1)];
                
                % generate all the "inputs" first; then loop thru time
                bias = bsxfun(@plus,mvnrnd(zeros(T*Ncases,2*Ndims),SigmaX)',muX);
                for t = 1:(T*Ncases-1)
                    Z(:,t+1) = A*Z(:,t) + bias(:,t);
                end
                clear bias;
                Z = shiftdim(reshape(Z,[2*Ndims,T,Ncases]),2);
                C = mat2cell(Z,40,2,ones(T,1));
                [filterdata(1:T).states] = C{:};
                clear C;
               
                % wrap the position elements
                src = state2stim(params.dynamics.C,Z,trajmin,trajmax,...
                    params.N,Ncases);
                
                
            otherwise
                error('no such EFH model -- jgm\n\n');
        end
        
        
        
        
        
        
        
        
        
        
        
    case 'controlled'
        % deprecated!  But this case is still in here to support the
        % results in 'results\finalwts\wts1dControlled140428.mat'
        
        %%% NB that this setup assumes that the Ndims=Ninputs!
        umin = params.smin(:,strcmp(params.mods,'Efference-Copy'));
        umax = params.smax(:,strcmp(params.mods,'Efference-Copy'));
        eta = 0.8;
        dt = params.dynamics.A(1,params.Ndims+1);
        
        % make sure the ICs won't result in Lissajous curves that
        % demand excessive controls
        Amp = eta*(trajmax-trajmin)/2;
        dc = (trajmax+trajmin)/2;
        ustarmax = dt*bsxfun(@times,Amp,v0'.^2)./...
            bsxfun(@minus,Amp.^2,bsxfun(@minus,x0',dc).^2);
        
        checkCtrls = @(uTensor)(~any(...
            sum(bsxfun(@gt, uTensor, umax)) &...
            sum(bsxfun(@lt, -uTensor, umin))));
        OKCTRLS = checkCtrls(ustarmax);
        
        while OKCTRLS == 0
            fprintf('excessive controls requested; resampling...\n')
            x0 = sampleStatePrior(muX0,SigmaX0,Ncases,...
                trajmin,trajmax,mrgn);
            v0 = sampleStatePrior(muV0,SigmaV0,Ncases,...
                trajmin,trajmax,mrgn);
            ustarmax = Amp.*v0'.^2./(Amp.^2 - (x0'-dc).^2)*dt;
            OKCTRLS = checkCtrls(ustarmax);
        end
        
        % the reference: tensor of Lissajous curves (Ncases x
        % Ndims x T)
        R = assembleLissajousTensor(x0,v0,params.thmin,...
            params.thmax,T,dt,eta);
        %%% R = assembleLissajousTensor2(4,params.thmin,params.thmax,Ncases,T,dt);
        
        % get controls that will track a noisy version of these
        [Z,U] = generateTrackingControls(R,umin,umax,...
            params.dynamics);
        
        % store output and input for encoding in PPCs
        C = params.dynamics.C;
        Ybar = shortdata(Ncases,3,(C*longdata(Z)')');
        src(:,:,strcmp(params.mods,params.NS),:) = Ybar;
        
        H = params.dynamics.H;
        Vbar = shortdata(Ncases,3,(H*longdata(U)')');
        src(:,:,strcmp(params.mods,'Efference-Copy'),:) = Vbar;
        
        % for debugging purposes
        TOPLOT = 0;
        if TOPLOT
            th0 = get2DOutline(params.thmin,params.thmax,42);
            for iCase = 1:Ncases
                
                % plot
                figure(1); clf; hold on;
                
                plot(squeeze(R(iCase,1,:)),squeeze(R(iCase,2,:)),'r');
                plot(squeeze(Ybar(iCase,1,:)),squeeze(Ybar(iCase,2,:)),'m');
                plot(th0(:,1),th0(:,2),'k');
                
                axis equal;
                title(num2str(iCase));
                hold off;
                pause()
            end
        end
        
        % store the state in the filterdata structure
        for t = 1:T
            filterdata(t).states = Z(:,:,t);
        end
        
        
        
        
        
        
    case 'bouncing'
        
        % loop through time
        for t = 1:(T-1)
            
            % propagate all trajectories forward one step
            Z(:,:,t+1) = (A*Z(:,:,t)' + mvnrnd(zeros(Ncases,2*Ndims),SigmaX)')';
            filterdata(t+1).states = Z(:,:,t+1);
            
            % loop through trajectories
            for traj = 1:Ncases
                Z(traj,:,t+1) = wallbounce(Z(traj,:,t+1),Z(traj,:,t),trajmin,trajmax);
            end
        end
        src = Z(:,1:Ndims,:);
        
        
        
        
    case 'repelling'
        %%% assumes second-order dynamics
        
        % precalculate
        aa = 0.01;
        wmin = trajmin*ones(1,Ncases,'like',trajmin);
        wmax = trajmax*ones(1,Ncases,'like',trajmax);
        %%%%%%%%% hack
        wmin = wmin + 0.2;
        wmax = wmax - 0.2;
        %%%%%%%%%
        F = zeros(Nstates,Ncases);
        FLAG = 1;
        
        while FLAG == 1
            FLAG = 0;
            
            % loop through time
            for t = 1:(T-1)
                
                % create a wall force
                F((Ndims+1):2*Ndims,:) = aa./(Z(:,1:Ndims,t)'-wmin) +...
                    aa./(Z(:,1:Ndims,t)'-wmax) -...
                    0.001*Z(:,(Ndims+1):2*Ndims,t)';
                %%% note this is extra damping!
                
                % propagate all trajectories forward one step
                Z(:,:,t+1) = (A*Z(:,:,t)' +...
                    mvnrnd(zeros(Ncases,Nstates),SigmaX)' + F)';
                filterdata(t+1).states = Z(:,:,t+1);
                
                % check for bounary crossings anyway (it can happen)
                if any(sum(bsxfun(@gt, Z(:,1:Ndims,t+1), trajmax'),2))||...
                        any(sum(bsxfun(@lt, Z(:,1:Ndims,t+1), trajmin'),2))
                    fprintf('went out of bounds...regenerating -- jgm\n');
                    % this is dubious, of course...
                    x0 = sampleStatePrior(muX0,SigmaX0,Ncases,trajmin,trajmax,mrgn);
                    v0 = sampleStatePrior(muV0,SigmaV0,Ncases,trajmin,trajmax,mrgn);
                    Z(:,:,1) = [x0 v0];
                    filterdata(1).states = Z(:,:,1);
                    FLAG = 1;
                    break
                end
            end
        end
        src = Z(:,1:Ndims,:);
        
        
    case 'quadraticbowl'
        
        qbA = 0.6;
        FLAG = 1;
        cntrs = repmat((range/2 + trajmin)',Ncases,1);
        while FLAG == 1
            FLAG = 0;
            
            % loop through time
            for t = 1:(T-1)
                
               
                % propagate all trajectories forward one step
                noise = mvnrnd(zeros(Ncases,2*Ndims),SigmaX)';
                u = Z(:,1:2,t) - cntrs;  
                Z(:,:,t+1) = (A*Z(:,:,t)' + noise)';
                Z(:,3:4,t+1) = Z(:,3:4,t+1) - qbA*sign(u).*u.^2;                 
                % z(:,3:4,t+1) = z(:,3:4,t+1) - qbA.*u.^5;   
                
                
                filterdata(t+1).states = Z(:,:,t+1);
                
                % check for bounary crossings anyway (it can happen)
                if any(sum(Z(:,1:2,t+1) > repmat(trajmax',Ncases,1),2))||...
                        any(sum(Z(:,1:2,t+1) < repmat(trajmin',Ncases,1),2))
                    fprintf('went out of bounds...regenerating -- jgm\n');
                    % this is dubious, of course...
                    x0 = sampleStatePrior(muX0,SigmaX0,Ncases,trajmin,trajmax,mrgn);
                    v0 = sampleStatePrior(muV0,SigmaV0,Ncases,trajmin,trajmax,mrgn);
                    Z(:,:,1) = [x0 v0];
                    filterdata(1).states = Z(:,:,1);
                    FLAG = 1;
                    break
                end
            end
        end
        src = Z(:,1:Ndims,:);
        
        
        
    case 'resetting'
        
        Z(:,:,1) = wall_reset(trajmin,trajmax,Ncases,muV0,SigmaV0);
        filterdata(1).states = Z(:,:,1);
        % loop through time
        for t = 1:(T-1)
            
            % propagate all trajectories forward one step
            Z(:,:,t+1) = (A*Z(:,:,t)' + mvnrnd(zeros(Ncases,2*Ndims),SigmaX)')';
            %filterdata(t+1).states = z(:,:,t+1);
            
            % get logical indices of trajectories that exceed the bounds
            sfuture = Z(:,1:Ndims,t+1);
            resetinds = sum((sfuture<repmat(trajmin(:)',Ncases,1)),2)|...
                sum((sfuture>repmat(trajmax(:)',Ncases,1)),2);
            
            % reset these trajectories
            if sum(resetinds)
                %x0 = sampleStatePrior(muX0,SigmaX0,sum(resetinds),varmin,varmax);
                %v0 = sampleStatePrior(muV0,SigmaV0,sum(resetinds),varmin,varmax);
                %z(resetinds,:,t+1) = [x0 v0];
                Z(resetinds,:,t+1) = wall_reset(trajmin,trajmax,...
                    sum(resetinds),muV0,SigmaV0,Z(resetinds,:,t));
            end
            
            %moved this to below the reset - BKD
            filterdata(t+1).states = Z(:,:,t+1);
            
            % store indices
            filterdata(t+1).RESTART = find(resetinds)';
            
        end
        src = Z(:,1:Ndims,:);
        
        
        
        
        
    otherwise
        fprintf('warning: trajectories are allowed ');
        fprintf('to leave the workspace -- jgm\n\n');
        
end

% tensor -> matrix
src = longdata(src);
%%% you may want to store higher derivatives of the state only....

end
%-------------------------------------------------------------------------%
%-------------------------------------------------------------------------%



%-------------------------------------------------------------------------%
function x = sampleStatePrior(mu,Sigma,Ncases,xmin,xmax,mrgn)

Ndims = length(xmin);

if any(isinf(Sigma(:)))
    x = scalefxn(rand(Ndims,Ncases,'like',mu),...
        zeros(size(xmin),'like',mu),ones(size(xmin),'like',mu),...
        xmin+mrgn,xmax-mrgn)';
else
    x = bsxfun(@plus,mu',randn(Ncases,Ndims)*chol(Sigma)');
end

end
%-------------------------------------------------------------------------%


%-------------------------------------------------------------------------%
function Swrapped = state2stim(M,X,smin,smax,N,Ncases)
% multiply by the emission matrix and then wrap

Sunwrapped = shortdata(Ncases,3,(M*longdata(X)')');
srange = N/(N-1)*(smax - smin);
Swrapped = bsxfun(@plus,bsxfun(@mod,bsxfun(@minus,Sunwrapped,smin'),...
    srange'),smin');

end
%-------------------------------------------------------------------------%

%%% check
% [min(min(src(:,1,:))), max(max(src(:,1,:))); min(min(src(:,2,:))), max(max(src(:,2,:)))]
% [xmin xmax]
