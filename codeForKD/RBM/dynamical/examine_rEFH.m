function examine_rEFH(wts,params)
% examine_rEFH  Exploratory analysis of rEFH weights, tuning curves
% USAGE:
%
%   load('results/finalwts/1DwrappingHVNdamped.mat')
%   examine_rEFH(wts,params)
%
%   load('results/finalwts/wts1DrEFH141009.mat')
%   load('results/nonfinalwts/wts1DrEFHoverdamped141030.mat')

%-------------------------------------------------------------------------%
% Revised: 11/12/14 (JGM)
%   -rewrote analyses in terms of lag and PD rather than b, m
%   -computed autocorrelation of the inputs
%   -functionized, rationalized etc.
% Created: ??/??/??
%   by BKD
%-------------------------------------------------------------------------%

%-------------------------------------------------------------------------%
% TO DO:
% (2) find peak of cross correlation...
%-------------------------------------------------------------------------%


% init for this computer
[~,machine] = system('hostname');
params.machine = strtrim(machine);
T = 1000;
params.dynamics.T = T;
tauMax = 30; %%% hard-coded

% the predictions of the NN
% LDSdataTrain = getLDSdata(T,params);
% ...



% filter
LDSdataTest = getLDSdata(params);
[V0,pSENSORYtest,pEFHtest] = EFHfilter(LDSdataTest,wts,params);

% for each hidden unit, get MI-maximizing tau, max-MI, and p(V=1|th=th*)
Th = squeeze(LDSdataTest.S(:,:,strcmp(params.mods,params.NS),:));
[taus,maxMI,Ptau] = getOptimalLag(V0,Th,params);

% get the indices of units that are really tuned for past positions
thr = 0.2;  % minimum mutual information (in bits) required
hidInds = maxMI>thr;
visInds = [hidInds; true(params.N,1)];

% plot the receptive fields, p(V=1|x_{t-tau*}), where tau* maximizes M.I.
[iTh,maxInds] = plotDelayedRFs(Ptau(hidInds,:),params);

% histogram of taus, with autocorrelation
ac2D = autocorrelationForCircularVars(Th,params);
histLagsWithAutocorrelation(taus(hidInds),ac2D,tauMax,8845);

% plot weight matrix, sorted by PD and lag
[~,it] = sort(taus(hidInds)');
plotSortedWeightMatrix(iTh,it,wts{1}([visInds;false],hidInds),params,109);


% look at RF paramsHVN
% histJointDstrbOfDelaysAndPos(maxInds)
% figure(651);
% PDs = scalefxn(maxInds,1,30,trajmin,trajmax);
% scatterhist(PDs,taus,'Nbins',10,'Kernel','overlay')
% axis tight


% "receptive fields" of input units
% sPatch = linspace(params.margin,params.margin + params.respLength,size(Ptau,2))';
% F = tiledTuning(sPatch,1,params)';
% figure(3535);
% imagesc([Ptau(iTh,:);F]);






% Now redo using lines fit to the RFs in th-thdot space

% check if MI and thdot-to-th-linear-fit give the same results
[m,b,Rsq]  = getrfparams(V0,LDSdataTest,params);

% get the indices of units that are really tuned as stripes in th-thdot
thr = 0.3;  % minimum Rsq required
hidIndsMB = (Rsq>thr)';
visIndsMB = [hidIndsMB; true(params.N,1)];

plotTausAgainstSlopes(taus,m,params,1301);
plotTausAgainstSlopes(taus(hidInds),m(hidInds),params,1302);
dt = plotTausAgainstSlopes(taus(hidIndsMB),m(hidIndsMB),params,1303);


%%% foo = 1:225;
%%% caca = foo(hidInds);
%%% caca(abs(taus(hidInds) - m(hidInds)'/dt)>1.5)


% sort etc
[~,ib] = sort(b(hidIndsMB));
[~,im] = sort(m(hidIndsMB));
%%% [~,ir] = sort(Rsq);
histLagsWithAutocorrelation(m(hidIndsMB)/dt,ac2D,tauMax,8846);
plotSortedWeightMatrix(ib,im,wts{1}([visIndsMB;false],hidIndsMB),params,110);



% look at individual RFs
Thdot = squeeze(LDSdataTest.Z(:,2,:));
Nthbins = 30;
Nthdotbins = 30;
[ThIncidence,ThdotIncidence,ThdotBins,thdotmin,thdotmax] =...
    getIncidences(Nthbins,Nthdotbins,Th,Thdot,params);



% % using mutual information 
% plotThThdotRFs(ThIncidence,ThdotIncidence,ThdotBins,thdotmin,thdotmax,...
%     V0(:,hidInds,:),m(hidInds),b(hidInds),it,params,989);
% pause()
% % using line-fitting
% plotThThdotRFs(ThIncidence,ThdotIncidence,ThdotBins,thdotmin,thdotmax,...
%     V0(:,hidIndsMB,:),m(hidIndsMB),b(hidIndsMB),im,params,989);




% sorting the hidden units by type
thinAroundZero = [9,16,31,33,55,69,82,84,99,119,130,141,154,190,...
    167,173,179,195,205,218,224];
smallOnes = [22,39,42,44,83,91,96,101,110,131,137,139,199,225];
multimodal = [27,93,126,140,181,220];
weirdLowMI = [20,62,48,79,94,103,122,123,145,155,156,161,169,187,191,...
    202,209,216,223];
weirdHighMI = [133,197,148];

nonWeird = ones(size(taus))';
nonWeird([weirdLowMI,weirdHighMI]) = 0;



keyboard


% Compare responses for *ideal* lagged-position tuning curves
X0 = getIdealResponses(ThIncidence,Ptau,taus,T,params);
fakah(ThIncidence,ThdotIncidence,ThdotBins,thdotmin,thdotmax,...
    V0,X0,m,b,nonWeird,params,1,989);
fakah(ThIncidence,ThdotIncidence,ThdotBins,thdotmin,thdotmax,...
    V0,X0,m,b,[weirdHighMI weirdLowMI],params,1,989);




% put the weird ones at the end, and plot all RFs
adjustedTaus = taus;
adjustedTaus([weirdLowMI,weirdHighMI]) = max(taus)+1;
[~,orderedInds] = sort(adjustedTaus);
[allRFsV,allRFsX] = fakah(ThIncidence,ThdotIncidence,ThdotBins,...
    thdotmin,thdotmax,V0,X0,m,b,orderedInds',params,0,989);
figure; showrfs(allRFsV);
figure; showrfs(allRFsX);



%%%%%%%%%%
X0 = getIdealResponses(ThIncidence,Ptau,taus,T,params);
[~,orderedInds] = sort(taus);
[allRFsV,allRFsX] = fakah(ThIncidence,ThdotIncidence,ThdotBins,...
    thdotmin,thdotmax,V0,X0,m,b,orderedInds',params,0,989);
figure; showrfs(allRFsV);
figure; showrfs(allRFsX);
%%%%%%%%%%

end













% % Steady-state analysis
% A = params.dynamics.A;
% C = params.dynamics.C;
% SigmaX = params.dynamics.SigmaX;
% P = dareJGM(A',C',LDSparamsEM2ndOrd.SigmaY(1),SigmaX);
% K = P*C'*inv(LDSparamsEM2ndOrd.SigmaY(1) + C*P*C');
% D = eye(2) - K*C;
% 
% clear caca;
% caca(:,1) = K;
% for iFoo = 1:25, caca = cat(2,K,D*A*caca); end
% 
% NNN = histc(taus,linspace(0,24,10)); bar(linspace(0,25,10),NNN,1); 
% axis tight; hold on;
% plot(caca(1,:)/caca(1,1)*NNN(1)); %%%
% hold off;







% % analytical derivation of autocorrelation, in continuous time
% % The autocorrelation of the output (y) of an LTI system to a white-noise
% % input is: 
% %
% %   R_y(tau) = L^{-1}{H(s)*H(-s)}
% %
% % For the harmonic oscillator *with k = 0*: 
% %  
% %   H(s) = 1/(s^2 + b/m*2)
% %
% % Hence (via partial-fraction expansion), the PSD (in Laplace var) is:
% %
% %   S_y(s) = -1/b^2*(1/s^2) - m/(2*b^3)*(1/(s+b/m)) + m/(2*b^3)*(1/(s-b/m))
% %
% % Taking an inverse Laplace transform,
% %
% %   R_y(t) = -1/b^2*t - m/(2*b^3)*exp{-b(t/m)}* + m/(2*b^3)*exp{b(t/m)}
% %
% % (times the Heaviside function....).
% T = params.dynamics.T;
% t = linspace(0,T*dt,T);
% plot(t,-t - 1/2*exp(-t));






% % autocorrelation, analytical, discrete-time
% sx2 = params.dynamics.SigmaX(1,1);
% sv2 = params.dynamics.SigmaX(2,2);
% alp = params.dynamics.A(2,2);
% dt = params.dynamics.A(1,2);
% gamma = sv2*dt^2;
% 
% clear Ryy;
% A = -1/(1-alp)^2;
% C = -alp/((1-alp)^2*(1-alp^2));
% 
% 
% n = -999:1:999;
% Ryy = -sx2*abs(n) + gamma*A*abs(n) + gamma*C*alp.^abs(n);
% 
% 
% figure(542); clf; hold on;
% plot(n*dt,Ryy);
% hold off;







% % check
% clc; 
% 
% z = [2 5 7];
% invz = z.^-1;
% foo = 1./((z-1).*(invz - 1).*(z - alp).*(invz - alp));
% foo2 = -invz./((1-invz).^2 .*(1-alp*invz).*(1-alp*z));
% foo3 = z.^2./(alp*z.^4 + (- alp^2 - 2*alp - 1)*z.^3 +...
%     (2*alp^2 + 2*alp + 2)*z.^2 + (- alp^2 - 2*alp - 1)*z + alp);
% 
% bb = [0 0 1 0 0];
% aa = [alp, (-alp^2 - 2*alp - 1), (2*alp^2 + 2*alp + 2),...
%     (-alp^2 - 2*alp - 1) alp];
% num = bb*[z.^4; z.^3; z.^2; z.^1; z.^0];
% den = aa*[z.^4; z.^3; z.^2; z.^1; z.^0];
% foo4 = num./den;
% 
% 
% [r,p,k] = residue(bb,aa);
% foo5 = sum(bsxfun(@rdivide,r,bsxfun(@minus,z,p)),1);
% 
% 
% fooInf = A*invz./(1-invz).^2 + C./(1-alp*invz) + C./(1-alp*z);
% 
% foo
% foo4
% foo5
% fooInf
% 
% % sum(abs(foo-poop))





% % cf. autocorrelationForCircularVars
% T = params.dynamics.T;
% 
% ac = 0;
% acWrapped = 0;
% for iTraj = 1:size(LDSdataTest.Z,1)
%     z1 = squeeze(LDSdataTest.Z(iTraj,1,:));
%     s1 = squeeze(LDSdataTest.S(iTraj,1,1,:));
%     
%     ac = ac + xcorr(z1,z1,T-1,'none');
%     acWrapped = acWrapped + xcorr(s1,s1,T-1,'none');
% end
% %%% ac = ac/size(LDSdataTest.Z,1); %%% get rid of?
% % acWrapped = acWrapped/size(LDSdataTest.Z,1);
% 
% figure(431); clf; hold on;
% plot(-(length(ac)-1)/2:(length(ac)-1)/2,ac); %%% ac/max(ac)) %%% change to ac???
% % plot(-(length(ac)-1)/2:(length(ac)-1)/2,acWrapped/max(acWrapped),'r')
% hold off;






% %% for saving to illustrator
% figure
% clims = [-.5 .5];
% imagesc([],[],out(:,1:params.N),clims);
% axis equal
% axis tight
% ylabel('Hidden Units_t')
% xlabel('Encoding Units')
% set(gca,'XTick',[1 15],'XTickLabel',[1 30],'TickLength',[.005 0],'TickDir','out')
% colormap(flipud(cbrewer('div','RdBu',256)))
% 
% figure
% imagesc([],[],out(:,params.N+1:end),clims);
% axis equal
% axis tight
% %ylabel('Hidden Units_t')
% xlabel('Hidden Units_{t-1}')
% set(gca,'TickLength',[.005 0],'TickDir','out','YTick',[])
% colormap(flipud(cbrewer('div','RdBu',256)))
% colorbar



%-------------------------------------------------------------------------%
%-------------------------------------------------------------------------%


%-------------------------------------------------------------------------%
function [iTh,maxInds] = plotDelayedRFs(Ptau,params)

% params
N = params.N;
trajmin = params.smin;
trajmax = N/(N-1)*(params.smax - params.smin) + params.smin;

% do things
[~,maxInds] = max(Ptau');
%%% adjustedPtau = Ptau; 
%%% adjustedPtau([weirdLowMI,weirdHighMI],:) = repmat(1:1:30,22,1);
[~,iTh] = sort(maxInds);
figure(5313);
imagesc(Ptau(iTh,:));
xlabel('$\prop$','Interpreter','none')
ylabel('neuron $\#$','Interpreter','none')
xticklabels = [trajmin, pi/2, trajmax]; %%% hard-coded
xticks = linspace(1, size(Ptau, 2), numel(xticklabels));
set(gca, 'XTick', xticks, 'XTickLabel', xticklabels)
colormap(cbrewer('seq','BuGn',256))
hold on;
plot(maxInds(iTh),1:length(maxInds),'k','Linewidth',2.0)
hold off;

end
%-------------------------------------------------------------------------%

%-------------------------------------------------------------------------%
function histLagsWithAutocorrelation(taus,ac2D,tauMax,fignum)

% params
Nedges = 15;

% do things
edges = linspace(0,tauMax,Nedges);
figure(fignum);
Ns = histc(taus,edges); bar(edges,Ns,1); axis tight; hold on;
truncatedAC = ac2D(:,(size(ac2D,2)+1)/2:(size(ac2D,2)+1)/2 + tauMax);
histArea = tauMax/Nedges*sum(Ns);
acArea = sum(truncatedAC,2);
normalizedAC = bsxfun(@rdivide,truncatedAC,acArea)*histArea;
%%% plot(0:edges(end),normalizedAC(1,:),'r')
plot(0:edges(end),mean(normalizedAC),'k','Linewidth',2.0); 
%%% avg across cos(th) and sin(th)
hold off;

end
%-------------------------------------------------------------------------%

%-------------------------------------------------------------------------%
function plotSortedWeightMatrix(iTh,it,W,params,fignum)

figure(fignum); clf;
subplot(2,1,1)
imagesc(W([iTh,(end-params.N+1):end],iTh)',[-0.5 .5]);
%%% hard-codes BP rather than PB
%%% hard-codes a color range to make the picture look interesting
axis xy
title('sorted by preferred (delayed) position')
colorbar('location','West')

subplot(2,1,2)
imagesc(W([it,(end-params.N+1):end],it)',[-0.5 .5]); 
%%% hard-codes BP rather than PB
%%% hard-codes a color range to make the picture look interesting
axis xy
title('sorted by optimal lag')
% colormap(flipud(cbrewer('div','RdBu',256)))

end
%-------------------------------------------------------------------------%

%-------------------------------------------------------------------------%
function dt = plotTausAgainstSlopes(taus,m,params,fignum)
if isfield(params.dynamics,'dt')
    dt = params.dynamics.dt;
else 
    dt = params.dynamics.A(1,2); % hard-coded for 1D
end
figure(fignum); clf; hold on;
scatter(m/dt,taus);
xlabel('steps via linfit');
ylabel('steps via MI');
plot(m/dt,m/dt,'r');
hold off

end
%-------------------------------------------------------------------------%

%-------------------------------------------------------------------------%
function [ThIncidence,ThdotIncidence,ThdotBins,thdotmin,thdotmax] =...
    getIncidences(Nthbins,Nthdotbins,Th,Thdot,params)

% present positions, past velocities
Th = Th(:,2:end);
Thdot = Thdot(:,1:end-1);

% Th edges are easy
N = params.N;
trajmin = params.smin;
trajmax = N/(N-1)*(params.smax - params.smin) + params.smin;
ThEdges = linspace(trajmin,trajmax,Nthbins+1);

% Thdot edges are more complicated
thdotmin = -3*sqrt(var(Thdot(:)));
thdotmax = 3*sqrt(var(Thdot(:)));
if min(Thdot(:)) > thdotmin, thdotmin = min(Thdot(:)) + eps; end
if max(Thdot(:)) < thdotmax, thdotmax = max(Thdot(:)) - eps; end
ThdotEdges = linspace(thdotmin,thdotmax,Nthdotbins+1-2);
ThdotEdges = [-inf ThdotEdges inf];

% Ns
Nsamples = length(Th(:));
[~,ThBins] = histc(Th(:),ThEdges);
[~,ThdotBins] = histc(Thdot(:),ThdotEdges);

% p(th,thdot)
ThIncidence = sparse(1:Nsamples,ThBins,1,Nsamples,length(ThEdges)-1);
ThdotIncidence = sparse(1:Nsamples,ThdotBins,1,Nsamples,length(ThdotEdges)-1);


end
%-------------------------------------------------------------------------%

%-------------------------------------------------------------------------%
function plotThThdotRFs(ThIncidence,ThdotIncidence,ThdotBins,...
    thdotmin,thdotmax,V0,m,b,unitVec,params,fignum)


% params
Ndots = 100;
N = params.N;
trajmin = params.smin;
trajmax = N/(N-1)*(params.smax - params.smin) + params.smin;

% p(th,thdot)
pThThdot = full(ThIncidence'*ThdotIncidence);	% omitting normalizer

% set the color map!
figure(fignum); 
% colormap(cbrewer('seq','PuBuGn',256));

% loop across hidden units
for iUnit = unitVec
            
    % compute p(V=1,th,thdot) and p(V=1|th,thdot)
    v0 = squeeze(V0(:,iUnit,2:end));
    Nsamples = length(v0(:));
    ThDotIncidenceTimesPofVis1 = sparse(1:Nsamples,ThdotBins,v0(:),...
        length(v0(:)),size(ThdotIncidence,2));
    PrVis1ThThdot = full(ThDotIncidenceTimesPofVis1'*ThIncidence)'; 
    %%% omitting normalizer
    PrVis1givenThThdot = PrVis1ThThdot./(pThThdot + eps);

    
    % plot p(V=1|th,thdot)
    clf; hold on
    imagesc([thdotmin,thdotmax],[trajmin,trajmax],PrVis1givenThThdot);
    axis tight
    
    % plot line fitting th to thdot
    thdotPts = linspace(thdotmin,thdotmax,Ndots)';
    thhat = m(iUnit)*thdotPts + b(iUnit);
    plot(thdotPts,thhat,'g.')
    
    % label
    title(iUnit)
    ylabel('position')
    xlabel('velocity')
    hold off
    
    pause(0.1)
    
end

end
%-------------------------------------------------------------------------%


%-------------------------------------------------------------------------%
function X0 = getIdealResponses(ThIncidence,Ptau,taus,T,params)

Nneurons = size(Ptau,1);
Ncases = params.Ncases;

PrXis1 = reshape(ThIncidence*Ptau',[Ncases,T-1,Nneurons]);
X0 = zeros(Ncases,Nneurons,T);
for iNeuron = 1:Nneurons
    thisTau = taus(iNeuron);
    X0(:,iNeuron,(2+thisTau):end) = PrXis1(:,1:(end-thisTau),iNeuron);
end

end
%-------------------------------------------------------------------------%

%-------------------------------------------------------------------------%
function plotSortedWeightMatrix2(iTh,it,W,params,fignum)


figure(fignum); clf; 
colormap(flipud(cbrewer('div','RdBu',256)))
imagesc(W(iTh,iTh)',[-0.5 .5]);
axis xy
colorbar('location','West')
keyboard
fhandle = figure(fignum); 
matlab2tikzWrapper('WfbSortedByPD',fhandle);
fignum = fignum + 10;

figure(fignum); clf; 
colormap(flipud(cbrewer('div','RdBu',256)))
imagesc(W((end-params.N+1):end,iTh)',[-0.5 .5]);
axis xy
fhandle = figure(fignum); 
matlab2tikzWrapper('WpropSortedByPD',fhandle);
fignum = fignum + 10;


figure(fignum); clf; 
colormap(flipud(cbrewer('div','RdBu',256)))
imagesc(W(it,it)',[-0.5 .5]);
axis xy
fhandle = figure(fignum); 
matlab2tikzWrapper('WfbSortedByTau',fhandle);
fignum = fignum + 10;


figure(fignum); clf; 
colormap(flipud(cbrewer('div','RdBu',256)))
imagesc(W((end-params.N+1):end,it)',[-0.5 .5]);
axis xy
fhandle = figure(fignum); 
matlab2tikzWrapper('WpropSortedByTau',fhandle);



% imagesc(W([it,(end-params.N+1):end],it)',[-0.5 .5]); 
% %%% hard-codes BP rather than PB
% %%% hard-codes a color range to make the picture look interesting
% axis xy
% title('sorted by optimal lag')
% % colormap(flipud(cbrewer('div','RdBu',256)))

end
%-------------------------------------------------------------------------%

%-------------------------------------------------------------------------%
function doubleWeightSort(Ptau,taus,wts,params)
% Weight matrix connectivity
% Try to sort by both preferred lag and preferred position at once.  In
% practice, this means binning preferred positions, sorting these bins,
% then sorting *within* these bins by preferred lag.




%%%%% NB: 
% MAY WANT TO REWRITE TO WORK WITH PRUNED WEIGHT MATRIX, TAUS, PTAU$
%%%%%%%%%%%%%%%%%%%%%





% resort the neurons by (lagged) position tuning
[iTh,maxInds] = plotDelayedRFs(Ptau,params);

% bin the position tunings
Nbins = 8;
edges = linspace(1,max(maxInds)+0.01,Nbins+1);
[Ns,bins] = histc(maxInds,edges);

% sort by optimal lag length *within* each position bin
[sortedVals,sortedInds] = arrayfun(@(iBin)(sort(taus(bins==iBin)')),...
    1:Nbins,'UniformOutput',false);
trueInds = arrayfun(@(iBin)(find(bins==iBin)),1:Nbins,...
    'UniformOutput',false);
sortedTrueInds = arrayfun(@(iBin)(trueInds{iBin}(sortedInds{iBin})),...
    1:Nbins,'UniformOutput',false);
allSortedTrueInds = [sortedTrueInds{:}];

% plot (transpose to put visibles on x-axis, hiddens on y-axis)
figure(1)
colormap(flipud(cbrewer('div','RdBu',256)))
imagesc(wts{1}(allSortedTrueInds,allSortedTrueInds)',[-1 1]); 
axis xy

figure(2)
colormap(flipud(cbrewer('div','RdBu',256)))
imagesc(wts{1}((end-params.N+1):end,allSortedTrueInds)',[-1 1]); 
axis xy


% other things
foo = wts{1}(allSortedTrueInds,allSortedTrueInds)';
figure(4)
imagesc((foo - foo')/2,[-0.05,0.05]);
axis xy;
centerBlocks = arrayfun(@(iBin)(...
    wts{1}(sortedTrueInds{iBin},sortedTrueInds{iBin}))',...
    1:Nbins,'UniformOutput',false);
blockAsymmetries = arrayfun(@(iBin)(...
    centerBlocks{iBin}-centerBlocks{iBin}'),...
    1:Nbins,'UniformOutput',false);


totalBlockAsymmetries = arrayfun(@(iBin)(sum(arrayfun(@(ii)(...
    sum(blockAsymmetries{iBin}(ii,ii+1:end))),...
    1:size(blockAsymmetries{iBin},1)))),1:Nbins)

end
%-------------------------------------------------------------------------%