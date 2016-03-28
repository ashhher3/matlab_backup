function [mw,mvb,mw2,mhb,b,k,Ts] = lrnrates(VISFXN,HIDFXN,params)
% LRNRATES  RBM learning rate adjuster
%
%   USAGE:
%
%	[mw,mvb,mhb,b,k,Ts] = lrnrates(VISFXN,HIDFXN,params)
%
%   LRNRATES sets the learning rates for the rbm based on unit types.  This
%   file was created just to unclutter EFH.m

%-------------------------------------------------------------------------%
% Adapted: 03/13/15
%   - changed for separate hid to vis weights
%   by KHPD
% Revised: 10/01/14
%   -rewritten in terms of mass-spring-damper parameters
% Cribbed: 01/28/11
%   -from rbm.m
%   by JGM
%-------------------------------------------------------------------------%


% default
mw = params.mw;
mvb = params.mvb;
mw2 = params.mw2;
mhb = params.mhb; %changed for nonsymm
b = params.b;
k = params.k;
Ts = params.Ts;

switch VISFXN
    case 'Binomial'
        mass = 1e10; % ???
        mw = mass;
        mvb = mass;
        mhb = mass;
        k = 0.0001; % 0.001
    case 'Bernoulli'
        mw = 5;
        mvb = 1.2;
        mhb = 1.2;
        b = 2.5;
        k = 0.0001; % ?
    case {'GP','Gaussian'}
        mw = mw*100;
        mvb = mvb*100;
        mhb = mhb*100;
        %     case {'BP'}
        %         rates = 200e-4;
        %         ew = rates;
        %         evb = rates;
        %         ehb = rates;
        % elseif strcmp(VISFXN,'Poisson') || strcmp(VISFXN,'PB')
        %     rates = .05e-4; % 1000e-5;
        %     ew =  rates;
        %     evb = rates;
        %     ehb = rates;
        %     wtcost = 2e-4 % 0.001;
        %     etaf = 0.8;
        % elseif strcmp(VISFXN,'GB')
        %     rates = 100e-4;
        %     ew = rates;
        %     evb = rates;
        %     ehb = rates;
        %     wtcost = 0.001; % 0.001
        %     etaf = 0.9;
end

% hidden layer of size 3/2 the input layer:
% 15 for 160 epochs; got to 9.75 with modes; ~5.3 w/means (secretTest)


if strcmp(params.machine,'domestica')
    mw = gpuArray(mw);
    mvb = gpuArray(mvb);
    mhb = gpuArray(mhb);
    b = gpuArray(b);
    k = gpuArray(k);
end


end








