function A = RMM(s, rd)
% A = RMM(s,rd)
%
% Implementation of the Raster Marginals Model, as introduced in
% "Population rate dynamics and multineuron firing patterns in 
%  sensory cortex", Journal of Neuroscience.
%
% Input: s,rd are row vectors of non-negative integers, such that
% sum(s) == sum((0:length(rd)-1).*rd)
%
% Output: A is a random 0/1 matrix of dimensions length(s) by sum(rd)
% such that sum(A(i,:)) == s(i) and the number of columns that contain
% exactly (j-1) 1s is rd(j). If no such matrix exists, A = [].
% 

% Copyright (c) 2012 Michael Okun, michael.okun@mail.huji.ac.il
% License: GPLv3


[ss, px] = sort(s, 'descend'); % px is permutation that takes s to ss.
px = sortrows([px; 1:length(s)]')'; 
px = px(2,:); % now px is the permutation that takes ss to s

A = Ryser(ss, rd); 
if isempty(A)
  return;
end;
A = A(px, :)'; % the canonical matrix with the provided constraints

% Now shuffle
A = Kshuffle(A)';

