function probs = SumProduct(tree, root, potentials, evidence)
% SumProduct.m - implemenbts sum-product algorithm for undirected trees
%
% Note that I'm assuming all variables take on the same number of possible
% values, as is the case in the example, but this code should still work
% if you make the matrices big enough for the variable that takes on the
% most values and then put NaN in all extra spaces.
%
% In all comments that follow, n is the # of nodes, and k is the # of 
% values a node can have.
%
% INPUTS:  tree - an nxn matrix whose i,j-th entry is 1 if there exists an 
%                 edge between nodes i and j, and 0 if no edge exists
%          root - an integer indicating the index of the root node in tree
%          potentials - an nxnxkxk matrix defining the clique potentials.
%                 Singleton potentials are stored as diagonal matrices.
%                 If singleton potentials don't exist, then entries are NaN
%          evidence - a vector of length n giving the values taken by any
%                 node in the evidence set. Nodes not in the evidence set
%                 must have a value of NaN.
% OUTPUTS: rootProbs - a vector of length k giving p(root|evidence)
%
% Kate Derosier, October 1, 2015 for CS281A/STAT241A

% do some basic input checking for my own good while debugging
a=size(tree,1);
b=size(tree,2);
c=size(potentials,1);
d=size(potentials,2);
e=size(evidence,1);
if (root>size(tree,1)) || ~(a==b&&a==c&&a==d&&a==e) ...
        || (size(evidence,2)~=1) 
    error('Something is wrong with input dimensions.')
elseif ~(isnan(evidence(root)))
    error('Root node is in evidence set.')
end
clear a b c d e

% turn inputs into global variables
global T
global P
global E
T=tree;
P=potentials;
E=evidence;

% compute and initialize other variables
global M
M=nan(size(T,1),size(T,1),size(P,3)); % nxnxk matrix to hold messages
probs=nan(size(T,1),size(P,3)); % nxk matrix to hold local marginals
nbrs=find(T(:,root)==1); % neighbors of root

% following the pseudocode from p.14 of ch.4 of textbook:
for i=1:numel(nbrs) % for neighbors of root
    Collect(nbrs(i),root);
end
for i=1:numel(nbrs)% for neighbors of root
    Distribute(root,nbrs(i));
end
for i=1:size(T,1) % for all nodes
    probs(i,:)=ComputeMarginal(i);
end

end

function []=Collect(j,i) % collect message from j to i
global T
nbrs = find(T(:,j)==1); % get all neighbors of j and remove i
nbrs = nbrs(nbrs~=i);

for h=1:numel(nbrs)
    Collect(nbrs(h),j) % collect messages from neighbors to j
end

SendMessage(j,i,nbrs) % send a message from j to i

end

function []=Distribute(i,j) % distribute messages from i to j
global T
nbrs = find(T(:,i)==1); % get all neighbors of i and remove j
nbrs = nbrs(nbrs~=j);

SendMessage(i,j,nbrs) % send message from i to j

nbrs = find(T(:,j)==1); % get all neighbors of j and remove i
nbrs = nbrs(nbrs~=i);

for h=1:numel(nbrs)
    Distribute(j,nbrs(h)) % distribute messages from j to neighbors
end
end

function []=SendMessage(j,i,nbrs) % send message from j to i
global M % we need to access messges, evidence, and potentials
global E
global P

% initialize vector for message from j to i for all values of i
sumVec = zeros(size(P,4),1); 

for v=1:size(P,4) % sum over all values of j
    if ( isnan(E(j)) )||( v==E(j) ) % if j in E, only use its true value
        prod = 1;
        for h=1:numel(nbrs)
            prod = prod * M(nbrs(h),j,v); % message from h to j when j has val v
        end
        
        if ~(isnan(P(j,j,v,v))) % singleton potential, if it exists
            prod = prod * P(j,j,v,v);
        end
        
        prodVec = prod * squeeze(P(i,j,v,:)); % edge potentials when j has value v
        
        sumVec = sumVec + prodVec;
        
    end
end

M(j,i,:)=sumVec; % record messages from j to i for all values of i

end

function p=ComputeMarginal(i)
global M %we need access to messages, evidence, potentials, and tree
global E
global P
global T

nbrs = find(T(:,i)==1); % get neighbors of i

prod=ones(size(P,4),1);
for h=1:numel(nbrs); % product over all neighbors of i
    prod = prod .* squeeze(M(nbrs(h),i,:));
end

if ~(isnan(P(i,i,1,1))) % singleton potential, if it exists
    prod = squeeze(P(i,i,:,:)) * prod;
end

if ~(isnan(E(i))) % if i in evidence
    for h=1:size(P,4)
        if h~=E(i)
            prod(h)=0; % throw out everything for all other values of i
        end
    end
end

p=prod/sum(prod); % normalize to make it a probability

end











