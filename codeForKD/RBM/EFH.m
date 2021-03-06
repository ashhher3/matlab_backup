% EFH   Exponential Family Harmonium training
%   EFH trains an exponential-family harmonium.  It requires the following
%   to have been defined outside this script:
%
%   Nhid        -- number of hidden units
%   batchdata   -- the data that are divided into batches
%                   (Ncases Ndims Nbatches)
%   restart     -- set to 1 if learning starts from beginning
%   params      -- the network parameters
%   datagenargs -- some extra parameters for the creation of data; you can
%                   initialize this to an empty cell if there's nothing

%-------------------------------------------------------------------------%
% Adapted: 03/13/15
%   - changed for separate hid to vis weights
%   by KHPD
% Revised: 01/15/15
%   -pulled out testEFHdecoding as a separate function
% Revised: 09/18/14
%   -renamed rbm.m -> EFH.m
% Revised: 05/14/13
%   -cleaned up
% Revised: 05/10/13
%   -added options for recurrency
% Revised: 09/14/10
%   -fixed tester to work with DBNs rather than just RBMs
% Revised: 05/25/10
%   -made number of contrastive-divergence steps variable
%   -incorporated rbmhidlinear
% Revised: 05/24/10
%   -fixed formatting
% modified by JGM
%-------------------------------------------------------------------------%



%---------------------------------- INIT ---------------------------------%
% for displaying
DISP = zeros(1,6);
if sum(DISP)
    close; figure(110); colormap(gray);
    rows = 4; cols = 5; space = 0;
    ax = getCustomAxesPos(rows,cols,space);
    indices = ceil(Nhid*rand(rows*cols,1));
end

% (re)init
if restart == 1
    
    % extract/set params
    HIDFXN = params.typeUnits{i_rbm+1};
    VISFXN = params.typeUnits{i_rbm};
    [mw,mvb,mw2,mhb,b,k,Ts] = lrnrates(VISFXN,HIDFXN,params); %changed for nonsymm
    maxepoch = params.DBNmaxepoch;
    amass = params.amass;
    
    datagenargs = [datagenargs,{'dbndepth',i_rbm,'dbnwts',wts}];
    
    % initialize weights - changed for nonsymm
    [vishid,hidvis,hidbiases,visbiases,vishidinc,hidvisinc,hidbiasinc, ...
        visbiasinc] = reinitializeRBMnonsym(i_rbm,params.numsUnits,wts);
    
    initialDiff=max(max(hidvis-vishid'));
    

    batchposhidmeans = zeros(Ncases,Nhid,Nbatches,'like',mw);
    allErrors = zeros(maxepoch,1,'like',mw);
    restart=0; epoch=1; erravg=0; tErravg=inf; trErravg=inf; counter=0;
    
    % plot errors
    if TESTDECODING
        
        % NB!!  If params.swing is 100%, the decoding error computed on
        % these data will start *increasing* after some point (e.g., epoch
        % 20).  Don't be fooled!  It may well still be decreasing, which
        % you can see by changing Rtest so that it is generated with (e.g.)
        % params.swing = 0.  This has happened before.
        
        if isfield(params,'dynamics')
            if ~exist('testData','var'), testData = getLDSdata(params); end
        else
            [Rtest,Stest] = DATAGENPP(Nbatches,params,datagenargs{:});
            [testData.R,testData.S] = longdata(Rtest,Stest);
            clear Rtest Stest
        end
        
        yvar = []; vvar = [];
        if usejava('desktop')
            setColors;
            figure(2014); clf; hold on;
            subplot(1,2,1); hold on;
            plotHandle(1) = plot(NaN,NaN);
            hold off;
            subplot(1,2,2); hold on;
            plotHandle(2) = plot(NaN,NaN);
            plotHandle(3) = plot(NaN,NaN);
            hold off;
        end
    end
    
end
%-------------------------------------------------------------------------%





%---------------------------------- LOOP ---------------------------------%
% cycle through training data
for epoch = epoch:maxepoch
    errsum = 0;
    
    % mass updates
    mw = amass*mw;
    mvb = amass*mvb;
    mw2 = amass*mw2;
    mhb = amass*mhb;
    % mw = 5000 - (5000 - 500)*amass^-epoch;
    % mvb = 5000 - (5000 - 500)*amass^-epoch;
    % mhb = 5000 - (5000 - 500)*amass^-epoch;
    
    
    % every Ntest epochs, make new data and maybe test the network
    if mod(epoch,Ntest)==1
        fprintf('.');
        [batchdata,Strain] = DATAGENPP(Nbatches,params,datagenargs{:});
        if TESTDECODING
            [yvar,tErravg] = testEFHDecoding(vishid,hidvis,hidbiases,...
                visbiases',yvar,testData,params); %changed for nonsymm
        end
    end
    [Ncases,Ndims,Nbatches] = size(batchdata);    %%% ?? why redo?
    
    
    for iBatch = 1:Nbatches
        fprintf('.');
        TRAJINIT = 0;
        
        % positive phase
        posdata = batchdata(:,:,iBatch);
        if isfield(params,'dynamics')&&(i_rbm==1)
            %%% ASSUMES RECURRENCE ON LEFT (e.g, BP RATHER THAN PB)
            if sum(sum(abs(posdata(:,1:Nhid)))) == 0
                TRAJINIT = 1;
            else
                posdata(:,1:Nhid) = poshidstates; % poshidmeans;
            end
        end
        poshidmeans = feedforward(posdata,vishid,hidbiases,HIDFXN,params);
        poshidstates = sampler(poshidmeans,HIDFXN,params);
        batchposhidmeans(:,:,iBatch) = poshidmeans;	%%% really use means for RBM2 data??
        posprods = posdata'*poshidstates;
        poshidact = sum(poshidstates);
        posvisact = sum(posdata);
        
        % negative phase
        [negvisstates,neghidstates] = CDstepperNONSYM(poshidstates,vishid,...
            hidvis,visbiases',hidbiases,HIDFXN,VISFXN,params);
        if TRAJINIT, negvisstates(:,1:Nhid) = 0; end
        negprods  = negvisstates'*neghidstates;
        neghidact = sum(neghidstates);
        negvisact = sum(negvisstates);
                
        % for printing (only)
        err = sum(sum((posdata - negvisstates).^2))/Ncases; % CHANGED
        errsum = err + errsum;
        
        
        % update weight/biases and their velocities simultaneously
        [vishid,vishidinc] = secondOrderWeightUpdate(vishid,...
            vishidinc,(posprods - negprods)/Ncases,mw,b,k,VISFXN,params);
        [hidbiases,hidbiasinc] = secondOrderWeightUpdate(hidbiases,...
            hidbiasinc,(poshidact - neghidact)/Ncases,mhb,b,k,HIDFXN,params);
        [visbiases,visbiasinc] = secondOrderWeightUpdate(visbiases,...
            visbiasinc,(posvisact - negvisact)'/Ncases,mvb,b,k,VISFXN,params);
        
        % "transpose" update for down weights
        [hidvis,hidvisinc] = secondOrderWeightUpdate(hidvis,...
            hidvisinc,(posprods - negprods)'/Ncases,mw2,b,k,VISFXN,params);
        
        % display
        if sum(DISP)
            EFHdisp(DISP,posdata,neghidstates,negvismeans,neghidprobs,...
                vishid,rows,cols,indices,params)
        end
        
    end
    fprintf('\n');
    erravg = errsum/Nbatches;
    allErrors(epoch) = erravg;
    
    % say error
    fprintf('epoch %4i error %6.4e trerror %6.4e terror %6.4e \n',...
        epoch,erravg,trErravg,tErravg);
    save rbmwts vishid hidvis hidbiases params epoch %changed for nonsymm
    
    if TESTDECODING&&usejava('desktop')
        set(0,'CurrentFigure',figure(2014));
        hold on;
        set(plotHandle(1),'XData',1:epoch,'YData',allErrors(1:epoch),'color','k');
        if mod(epoch,Ntest)==1
            set(plotHandle(2),'XData',find(mod(1:epoch,Ntest)==1),'YData',yvar,'color','r');
            % set(plotHandle(3),'XData',(Ntest:Ntest:epoch)-Ntest,'YData',vvar,'color','g');
        end
        hold off;
        %%% set(subplotHandle,'XLim',[0 Max]);
        %%% if you want to keep the x or y limits fixed over the animation....
    end
    
    
end

finalDiff=max(max(hidvis-vishid'));

%-------------------------------------------------------------------------%
%-------------------------------------------------------------------------%





%-------------------------------------------------------------------------%
% Version 1.000
%
% Code provided by Geoff Hinton and Ruslan Salakhutdinov
%
% Permission is granted for anyone to copy, use, modify, or distribute this
% program and accompanying programs and documents for any purpose, provided
% this copyright notice is retained and prominently displayed, along with a
% note saying that the original programs are available from our web page.
% The programs and documents are distributed without any warranty, express
% or implied.  As the programs were written for research purposes only,
% they have not been tested to the degree that would be advisable in any
% important application.  All use of these programs is entirely at the
% user's own risk.

% This program trains Restricted Boltzmann Machine in which visible,
% binary, stochastic pixels are connected to hidden, binary, stochastic
% feature detectors using symmetrically weighted connections. Learning is
% done with 1-step Contrastive Divergence. The program assumes that the
% following variables are set externally:
%
%   maxepoch    -- maximum number of epochs
%   Nhid      -- number of hidden units
%   batchdata   -- the data that is divided into batches
%   (Ncases Ndims Nbatches)
%   restart     -- set to 1 if learning starts from beginning
%-------------------------------------------------------------------------%
