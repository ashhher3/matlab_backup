function [ pwr ] = getPwr10Sec( chan, Fs )

w=10*Fs; % 10sec is w samples long

pwr=zeros(1,floor(length(chan)/w));

for i=1:w:length(chan)-w
    pwr((i-1)/w+1)=mean(chan(i:i+w-1));
end


end