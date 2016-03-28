function [vishid,hidvis,hidbiases,visbiases,vishidinc,hidvisinc,...
    hidbiasinc,visbiasinc] = reinitializeRBMnonsym(i_rbm,numsUnits,wts)
% reinitializeRBM   Initialize weights, biases, and increments of RBM

%-------------------------------------------------------------------------%
% Adapted: 03/13/15
%   - changed for separate hid to vis weights
%   by KHPD
% Revised: ??/??/14
%   -changed allocation for GPU compatibility
% Created: ??/??/?? (very early)
%   by JGM
%-------------------------------------------------------------------------%


% init
Nvis = numsUnits(i_rbm);
Nhid = numsUnits(i_rbm+1);
[~,machine] = system('hostname');
machine = strtrim(machine);
if strcmp(machine,'domestica'),
    yrclass = 'gpuArray';
else
    yrclass = 'double';
end

% initialize symmetric weights and biases.
vishid      = 0.01*randn(Nvis, Nhid,yrclass);
hidvis      = 0.01*randn(Nhid, Nvis,yrclass); %changed for nonsymm
hidbiases   = 0*ones(1,Nhid,yrclass);
visbiases   = 0*ones(Nvis,1,yrclass);

%%% doesn't happen for 1Dinteg, so it's not updated for nonsymm case
% if i_rbm > 1
%     if Nhid==numsUnits(i_rbm-1)
%         numRBMs     = length(numsUnits)-1;
%         vishid      = wts{2*numRBMs-i_rbm+2}(1:end-1,:);% i.e. hidvis
%         hidbiases   = wts{2*numRBMs-i_rbm+2}(end,:);    % i.e. the visbiases
%         visbiases   = wts{i_rbm-1}(end,:)';             % i.e. the hidbiases
%     end
% end

vishidinc	= zeros(Nvis,Nhid,yrclass);
hidvisinc	= zeros(Nhid,Nvis,yrclass);
hidbiasinc  = zeros(1,Nhid,yrclass);
visbiasinc  = zeros(Nvis,1,yrclass);



%%% when did you ever use these??
% poshidmeans     = zeros(numcases,numhid);
% neghidstates    = zeros(numcases,numhid);
% posprods        = zeros(numdims,numhid);
% negprods        = zeros(numdims,numhid);


end
