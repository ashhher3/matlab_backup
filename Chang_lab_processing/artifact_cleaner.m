% [Y, Z, P, E] = artifact_cleaner(X, i, a)
%
% X  input data matrix, (n by t) for t >> n
% i  indices (or logical indices) into X to clean (eg artifact times)
% a  how aggressively to clean artifact, integer from 1 to n
%
% Y  data with the artifact cleaned (n by t)
% Z  data projected into artifact-space (n by t)
% P  projection matrix (n by n)
% E  eigenvalues of projection matrix (n by 1)
%
% JE ODoherty v3.1

function [Y, Z, P, E] = artifact_cleaner(X, i, a)

    if nargin < 3
        disp('  ');
        disp('  Y  = artifact_cleaner(X, i, a)');
        disp('  ');
        disp('  X  input data matrix, (n by t) for t >> n');
        disp('  i  indices (or logical indices) into X to clean (eg artifact times)');
        disp('  a  how aggressively to clean artifact, integer from 1 to n');
        disp('  ');
        disp('  Y  data with the artifact cleaned (n by t)');
		disp('  ');
		disp(' JE ODoherty v3.1');
        return;
    end
    
    [n, t] = size(X);
    
    if n >= t
        error('data matrix, X, must be n by t for t >> n');
    end

    if isempty(i)
        error('indices, i, cannot be empty');
    end
    
    if islogical(i)
        if length(i) ~= t
            error('logical indices i are not of same length as X');
        end
        elseif (max(i) > t) || (min(i) < 1)
            error('indices, i, are invalid');
    end
    
    if length(a) > 1
        error('parameter, a, must be a scalar');
    end
    
    if (a > n) || (a < 1)
        error('parameter, a, must be a positive integer less than n)'); 
    end

	%% compute prefilter means	
	mx = nanmean(X, 2);

	%% compute covariances
	Cx = cov(X');
	Ca = cov(X(:,i)');

	%% get projection matrix
	[P, E] = eig(Cx, Ca, 'vector'); % nb, second output is necessary!

    %% get dimensions to kill
    a = floor(a);
	j = 1:a;
	G = eye(n);
	G( (j-1)*(n+1)+1 ) = 0;
    
	%% project
	Z = P' * X;
    
	%% clean and backproject
	Y = P' \ (G * Z);

	%% compute postfilter means
	my = nanmean(Y, 2); 

	% subtract postfilter means and add prefilter means
	m = mx - my;
	Y = bsxfun(@plus, Y, m);
	
end
