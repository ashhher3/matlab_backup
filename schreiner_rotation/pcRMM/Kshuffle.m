function A = Kshuffle(A)
%Kshuffle - a function for shuffling a raster
%  Input: A - a 0/1 matrix of time-by-N dimensions (N=no. of spike trains) 
% Output: A - the same matrix shuffled
% The function performs what is called "spike exchange across neurons"
% shuffling. 

% It is called Kshuffling because it is explained in fig. 9K 
% of S. Grun's JNeurophysiol'09 paper which discusses different methods

% Copyright (c) 2015 Michael Okun, michael.okun@mail.huji.ac.il
% License: GPLv2


for i = 1:10*nchoosek(size(A,2),2)  
  c = ceil(rand(1,2)*size(A,2)); % two randomly selected columns
  
  I = A(:,c(1)) + A(:,c(2)) == 1; % where the 2 columns don't coincide
  cA = A(I,[c(1) c(2)]); % a copy of the part that matters, to make it run faster
  i01 = find(cA(:,1) == 0);
  i10 = find(cA(:,1) == 1);  
    
  toFlip = ceil(min(length(i01), length(i10))/2); % how many 01s & 10s to flip
  
  i01 = i01(randperm(length(i01)));
  i01 = i01(1:toFlip);
  i10 = i10(randperm(length(i10)));
  i10 = i10(1:toFlip);
  
  % the flip itself:
  cA(i01,1) = true; cA(i01,2) = false;
  cA(i10,1) = false; cA(i10,2) = true;
  A(I,[c(1) c(2)]) = cA;  
end;


