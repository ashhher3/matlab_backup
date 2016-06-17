function [ bCoh ] = getCoh10Sec( chan1, chan2, Fs)

w=10*Fs; % 10sec is w samples long

bCoh=zeros(1,floor(length(chan1)/w));
params.Fs=Fs;
params.fpass=[13 30];

for i=1:w:length(chan1)-w
    C=coherencyc(chan1(i:i+w-1),chan2(i:i+w-1),params);
    bCoh((i-1)/w+1)=mean(C);
end


end