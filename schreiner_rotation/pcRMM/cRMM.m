function A = cRMM(s, prd, ccpr)
% A = cRMM(s,prd, ccpr)
%
% Implementation of the Raster Marginals Model with coupling terms, 
% as introduced in "Diverse coupling of neurons to populations in
% sensory cortex", Nature 2015.
%
% Input: s,prd,ccpr are row vectors of non-negative integers, such that
% sum(s) == sum((0:length(prd)-1).*prd), the idea being that s represents the
% total no. of spikes emitted be each unit, and prd is the population rate
% distribution
% ccpr provides the inner product of each unit with the population rate (PR)
% (see README.m for an example), thus in particular it must hold that:
% sum(prd.*(0:length(prd)-1).^2) = sum(ccpr)
%  
% Output: A is a random 0/1 matrix of dimensions length(s) by sum(prd)
% such that sum(A(i,:)) == s(i) and the number of columns that contain
% exactly (j-1) 1s is rd(j). In addition the constraint on inner products
% is satisfied (up to a small error).
% If no such matrix exists, A = [].
% 

% Copyright (c) 2015 Michael Okun, michael.okun@mail.huji.ac.il
% License: GPLv2

N = length(s);

% ---- Make the canonical matrix with given marginals
[ss, ix] = sort(s, 'descend'); % ix is permutation that takes s to ss.
ix = sortrows([ix; 1:length(s)]')'; 
ix = ix(2,:); % now ix is the permutation that takes ss to s

A = Ryser(ss, prd); 
if isempty(A)
  return;
end;
if sum(prd.*(0:length(prd)-1).^2) ~= sum(ccpr)
  error('cRMM: bad/incompatible inner product parameter');
end;
A = A(ix, :)'; % the canonical matrix with the provided constraints

% ---- Now shuffle
A=Kshuffle(A);

% ---- Finally adjust so that each unit respects also its inner product
%      with PR constraint

prevStep = size(A,1)*size(A,2)+7;
problemCounter = 0;
while true  
  if problemCounter > N
    error('cRMM: failed...');    
  end;
  PR = sum(A, 2);
  ip = PR'*A; % current inner product of units & PR
  assert(problemCounter> 0 || prevStep > sum(abs(ip-ccpr))); % we should be improving (i.e. getting ip closer to ccpr, which is the target) with time
  prevStep = sum(abs(ip-ccpr));
  poor = find(ccpr > ip+N); % these units should be more correlated (with PR) than they are (differences of N and less are really small to bother)
  rich = find(ccpr < ip-N); % these units should be less correlated than they are
  if isempty(poor) || isempty(rich)
     break; % if no one is really poor, no one can be really rich (and vice versa), so we're done.
  end;
  rich = rich(randperm(length(rich))); rich = rich(1); % pick some rich guy
  poor = poor(randperm(length(poor))); 
  for i = 1:length(poor)+1
    if any(A(:,poor(i)) & ~A(:,rich)) % we don't want a poor unit that spikes only when the rich one does. Such poor unit is useless
      break;
    end;
  end;
  if i == length(poor)+1
    problemCounter = problemCounter + 1;
    continue;
  end;
  poor = poor(i); % pick a useful poor unit
  
  others = setdiff(1:N, [poor rich]); % all the other units
  
  shiftRows = ceil(min([ccpr(poor)-ip(poor) ; ip(rich)-ccpr(rich)])/N); 
  % How many rows we want to shift. Each row will reduce the "unequality"
  % by at least 1, so we will fix at least 1/N fraction...
  
  % from the 'perspective of the poor'
  looseValues = (sum(A(:,others), 2)+1) .* ~A(:,rich) .* A(:,poor); % columns with spike of poor & no spike of rich, and how much each will loose
  gainValues = (sum(A(:,others), 2)+1) .* A(:,rich) .* ~A(:,poor);   
  
  cutValue = min(looseValues(looseValues >0)); % poor will loose values upto this (we expect it to be 1 except for extreme cases)
  
  loosePositions = find(looseValues > 0 & looseValues <= cutValue); 
  gainPositions = find(gainValues > cutValue);
 
  if isempty(gainPositions)
    % we are in a "strange" situation where the rich unit has only columns
    % with small no. of 1s, whereas the poor unit already has only columns
    % with big no. of 1s
    problemCounter = problemCounter+1;
    continue;
  end;    
  
  shiftRows = min([shiftRows; length(loosePositions); length(gainPositions)]);
  
  loosePositions = loosePositions(randperm(length(loosePositions)));
  gainPositions = gainPositions(randperm(length(gainPositions)));  
  loosePositions = loosePositions(1:shiftRows);
  gainPositions = gainPositions(1:shiftRows);  
  
  A(loosePositions, rich) = true;
  A(loosePositions, poor) = false;
  A(gainPositions, rich) = false;
  A(gainPositions, poor) = true;  
end;

assert(~any(abs(sum(A) - s)));

A = A';


