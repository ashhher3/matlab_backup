function [ bCoh ] = getCohNSamp( chan1, chan2, Fs, N )

bCoh=zeros(1,floor(length(chan1)/N));
params.Fs=Fs;
params.fpass=[13 30];

for i=1:N:length(chan1)-N
    C=coherencyc(chan1(i:i+N-1),chan2(i:i+N-1),params);
    bCoh((i-1)/N+1)=mean(C);
end


end