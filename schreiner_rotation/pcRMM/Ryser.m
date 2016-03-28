function A = Ryser(s, rd)
% A = Ryser(s,rd)
%
% Input: s,rd are row vectors of non-negative integers, such that:
%    (i) sum(s) == sum((0:length(rd)-1).*rd)
%   (ii) s is nonincreasing
%
% Output: A is the 'canonical' 0/1 matrix of dimensions length(s) by sum(rd)
% such that sum(A(i,:)) == s(i) and the number of columns that contain
% exactly (j-1) 1s is rd(j). If no such matrix exists, A = [].
%
% Example: Ryser([6 2 2 1], [1 1 2 2]) means we want a 4-by-6 0/1 matrix,
% whose rows sum to 6,2,2,1 and its columns sum to 0,1,2,2,3,3. 
% Such a matrix does not exist, so the output is []. On the other hand, 
% Ryser([4 3 2 2], [1 1 2 2]) produces an output matrix.

% The matrix A is constructed using Ryser's recursive algorithm, e.g.,
% see "Algorithms for constructing (0, 1)-matrices with
% prescribed row and column sum vectors", R. A. Brualdi,
% Discrete Mathematics 306 (2006) 3054-3062.
%
% Copyright (c) 2012 Michael Okun, michael.okun@mail.huji.ac.il
% License: GPLv3



if nargin < 2 || sum(s) ~= sum((0:length(rd)-1).*rd) || ...
    sum(s<0) > 0 || sum(rd<0) > 0  || sum(diff(s) > 0) > 0 || ...
    size(s, 1) > 1 || size(rd, 1) > 1 
  error('Ryser: bad input parameters');
end;

n = length(s); % no. of rows in A
m = sum(rd); % no. of columns in A

if length(rd) < n+1
  rd = [rd, zeros(1, n + 1 - length(rd))];
end;

% Let R denote the vector sum(A), i.e., the sum of columns in the matrix A
% we construct. R's conjugate (terminology as in Brualdi's paper) is:
Rc = cumsum(rd(end:-1:2)); Rc = Rc(end:-1:1);
Rc = [Rc, zeros(1, n-length(Rc))];

% Rc must dominate s, which is both required and sufficient for A to exist.
if  sum(cumsum(Rc)-cumsum(s) < 0) > 0
  A = [];
  return;
end;

A = false(n,m);
% First we build matrix with the required column sums, with 1s in their
% topmost positions. Columns with more 1s are to the left.
for i = length(rd):-1:2
  A(1:i-1, sum(rd(i+1:end))+1:sum(rd(i:end))) = true;
end;

% Next, we fix the row sums, recursively, one row a time (from bottom up).
% Once a row is fixed, we can "forget" about it, and focus on the submatrix
% without this row.
for r = n:-1:2
  d = s(r) - sum(A(r,:)); % how many 1s we need to move from above to row r
  %assert(d >= 0);
  
  rdC = histc(sum(A(1:r,:), 1), 0:1:r); 
  % rdC is the 'current rd', i.e., for the submatrix of A formed by the top
  % r rows, rdC tells how many columns sum up to 0, 1, ..., r-2, r-1, r.
  
  for j = r-1:-1:1
    % first we will move 1s from columns that sum to r-1, then (if
    % necessary) from columns that sum to r-2, and so on...
    if rdC(j+1) >= d
      %assert(sum(A(j, sum(rdC(j+1:end))-d+1:sum(rdC(j+1:end)))) == d);
      %assert(sum(A(r, sum(rdC(j+1:end))-d+1:sum(rdC(j+1:end)))) == 0);
      A(j, sum(rdC(j+1:end))-d+1:sum(rdC(j+1:end))) = false;
      A(r, sum(rdC(j+1:end))-d+1:sum(rdC(j+1:end))) = true;            
      break; % we are done ==> exit the loop
    else
      %assert(sum(A(j, sum(rdC(j+2:end))+1:sum(rdC(j+1:end)))) == rdC(j+1));
      %assert(sum(A(r, sum(rdC(j+2:end))+1:sum(rdC(j+1:end)))) == 0);
      A(j, sum(rdC(j+2:end))+1:sum(rdC(j+1:end))) = false;
      A(r, sum(rdC(j+2:end))+1:sum(rdC(j+1:end))) = true;            
      d = d - rdC(j+1); % the no. of 1s that we still have to take down 
    end;    
  end;  
end;

% Finally we check the output satisfies the requirements:
assert(sum(abs(sum(A, 2) - s')) == 0); % A fulfills the condition on rows
rdA=histc(sum(A), 0:1:n); 
assert(sum(abs(rd-rdA)) == 0); % A fulfills the condition on columns

