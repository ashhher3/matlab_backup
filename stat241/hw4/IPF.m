function [ psi, psi_init, p, pTilde ] = IPF( data, adj, states)
% Implementation of IPF algorithm. Takes data and adjacency matrix and
% returns potentials for the model 
%        p(x) = 1/Z prod_{(i,j) in E} psi_{i,j}(x_i, x_j)
%
% Maybe if I was more clever at Matlab this wouldn't be required, but you
% also have to pass in a big matrix showing all the states the graph can
% be in. (can't figure out how to do it without a variable # of for loops)
%
% Potentials and empirical distribs for (x_i, x_j) are saved as 2x2 
% matrices of the shape: 
%              [ i=0 j=0  ,  i=0 j=1;
%                i=1 j=0  ,  i=1 j=1 ]
%

% get the sizes of things
numV = size(adj, 1);
numE = 0.5*sum(sum(adj));
numD = size(data,1);
numS = size(states,1);

% randomly initialize potentials
[psi,numP]=initializePsi(numV,adj);
psi_init = psi;

% calculate empirical distribution
pTilde=empirical(states,numS,data,numD);

% run IPF
d = inf;
c = 0;
count = 0;
while count<10000
    newPsi = psi;
    c = mod(c, numE)+1;
    
    p = probFromPotential(states,numS,psi,numP);
    
    i=psi(c).pair(1);
    j=psi(c).pair(2);
    
    joint = pairFromFullDistrib(p,numS,i,j);
    jointTilde = pairFromFullDistrib(pTilde,numS,i,j);
    
    newPsi(c).potential=psi(c).potential .* (jointTilde ./ joint);
    
    psi = newPsi;
    count=count+1;
end

end

function d = Dkl(dist1,dist2,numS)
% THIS DOESN'T WORK AND I DON'T KNOW WHY
d=0;
for s=1:numS
    d=d+(dist1(s).prob*log(dist1(s).prob/dist2(s).prob));
end

end


function joint = pairFromFullDistrib(dist,numS,i,j)
if j<i % check that we're following the i,j convention to match psi
    error('xi must be the smaller index')
end

joint = zeros(2,2); % initalize empty matrix

for s=1:numS % add up the states
    state = [dist(s).state];
    joint(state(i)+1,state(j)+1)=dist(s).prob+joint(state(i)+1,state(j)+1);
end

end


function p = probFromPotential(states,numS,psi,numP)
% compute full prob distribution from model and store in struct
p = struct('state',[],'prob',[]); 

for s=1:numS % input all the states
    p(s).state = states(s,:);
    p(s).prob = 1; % initialize to 1 bc we will multiply
    
    for e=1:numP % product over edges
        i = psi(e).pair(1); % get indices of the edge
        j = psi(e).pair(2);
        
        % get state of the edge
        xi = p(s).state(i);
        xj = p(s).state(j);
        
        % look up val in potential table and multiply
        p(s).prob = p(s).prob*psi(e).potential(xi+1, xj+1);
        
    end
end

% normalize
Z = sum([p(:).prob]);
for s=1:numS
    p(s).prob = (1/Z)*p(s).prob;
end

end

function [psi, numP]= initializePsi(numV,adj)
% initialize psi
psi = struct('pair',[],'potential',[]);
idx = 0;
for i = 1:numV % for each edge, enter it and make array of right size
    for j = i+1:numV
        if adj(i,j)==1
            idx = idx+1;
            psi(idx).pair = [i,j];
            psi(idx).potential = rand(2,2);
        end
    end
end
numP = numel(psi);
end

function pTilde = empirical(states,numS,data,numD)
% compute complete empirical distribution and store in struct
pTilde = struct('state',[],'prob',[]); 

for s=1:numS % input all the states
    pTilde(s).state = states(s,:);
    pTilde(s).prob = 0; % initialize to zero bc we will add
end

for t=1:numD % count how many times each state appears in data
    [~,s]=ismember([data(t,:)],states,'rows');
    pTilde(s).prob = 1 + pTilde(s).prob;
end

for s=1:numS % normalize
    pTilde(s).prob = (1/numD)*pTilde(s).prob;
end

end