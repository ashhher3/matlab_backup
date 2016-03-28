function [cntrsOfMass, ttlSpks] = GTPNsuffstats(R,params)
% GTPNsuffstats     Sufficient statistics of Gaussian-tuned Poisson neurons
%
%   USAGE:
%       [cntrsOfMass, ttlSpks] = GTPNsuffstats(R,params);
%
% Given a matrix of GTPN populations (Nexample x Nmods*N^Ndims) and the
% params structure, compute the centers of mass (Nexamples x Ndims x Nmods)
% and the total spike counts (Nexamples x Nmods).

%-------------------------------------------------------------------------%
% Revised: 07/09/14
%   -transposed output ttlSpks
% Revised: 07/07/14
%   -rewrote from scratch to work on tensors rather than using loops
% Revised: 01/03/14
%   -fixed addition from last revision
% Revised: 12/18/13
%   -now works on data PPCs containing more than one "modality."  NB: this
%   changes the shape of the structure fields, Shat and CvrnMat!!!
% Revised: 07/02/13
%   -added input argument tuningCov, which improves the logic of the
%   function (it shouldn't have itself to decide to use the prop cov).
% Renamed: 07/02/13
%   -from PPC.m to GTPNsuffstats.m
% Cribbed: 06/28/13
%   from KF4PPC
%   by JGM
%-------------------------------------------------------------------------%

% Ns
Nmods = params.Nmods;
Ndims = params.Ndims;
N = params.N;
Nexamples = size(R,1);
Nlattice = N^Ndims;

% (Nexamples x Nmods*N^Ndims) -> (Nmods*Nexamples x N^Ndims)
Rtall = reshape(permute(reshape(R,[Nexamples,Nlattice,Nmods]),[3,1,2]),...
    [Nmods*Nexamples,Nlattice]);
ttlSpksTall = Rtall*ones(Nlattice,1,'like',R);

% construct "preferred directions" on an Ndims lattice
latticePts = ones(1,1,'like',R):N;
lattices1Dcell = mat2cell(latticePts(ones(Ndims,1,'like',R),:),...
    ones(Ndims,1),N);
lattice = ndgrid(lattices1Dcell{:});
latticePDs = NaN(Nlattice,Ndims,'like',R);
for iDim = 1:Ndims
    thisLattice = shiftdim(lattice,iDim-1);
    latticePDs(:,iDim) = thisLattice(:);
end

% for data on the Ndims-torus, center the data first
if isfield(params,'dynamics')
    if strcmp(params.dynamics.walls,'wrapping')
        [Rtall,gridShifts] = centerNoisyHills(Rtall,params.N); 
    end
else
    gridShifts = zeros(Nmods*Nexamples,Ndims);
end

% calculate center of mass
cntrsOfMassTall = bsxfun(@rdivide,Rtall*latticePDs,ttlSpksTall) - gridShifts;

% put back into standard tensor shapes (Nexamples x Ndims x Nmods)
cntrsOfMass = permute(reshape(cntrsOfMassTall,[Nmods,Nexamples,Ndims]),[2 3 1]);
ttlSpks = reshape(ttlSpksTall,[Nmods,Nexamples])';

% in case of no spikes!
badCoMs = isnan(cntrsOfMass);
cntrsOfMass(badCoMs) = (N-1)*rand(sum(badCoMs(:)),1,'like',R) + 1;


if sum(badCoMs(:))>0
    fprintf('warning: found %d bad centers of mass; ',sum(badCoMs(:)));
    fprintf('replacing with guesses -- jgm\n');
end

% convert into real-world units
for iMod = 1:Nmods
    cntrsOfMass(:,:,iMod) = grid2world(cntrsOfMass(:,:,iMod)',...
        [params.smin(:,iMod) params.smax(:,iMod)],params)';
end



end