function pc = computePC(raster)
% function to compute population coupling values for given raster matrix
% based on pcRMM code from Okun
%
% rows in raster should be units, columns should be timebins
%
% (KD 2015/05/20)

n=size(raster, 1);
pc=zeros(1,n);

for i=1:n %loop over all neurons
    %pop rate is sum of all other neurons
    popRate=sum([raster(1:i-1,:);raster(i+1:end,:)]);
    
    %(discrete) pop coupling is inner product w pop rate, normalized by FR
    pc(i)=(popRate*raster(i,:)')/(sum(raster(i,:)));
end

end


