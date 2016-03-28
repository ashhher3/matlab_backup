function data = getPeakTimes( data )
time=0:1/data.FsB(1):length(data.contact(1).beta)/data.FsB(1)-1/data.FsB(1);

for i=1:length(data.contact)
    burstBeta=data.contact(i).betaPwr([data.contact(i).burst]==1);
    burstTimes=time([data.contact(i).burst]==1);
    
    % find peaks of beta bursts
    [pks,locs]=findpeaks(burstBeta);
    peakTimes=burstTimes(locs);
    
    data.contact(i).burstPeaks=pks;
    data.contact(i).burstTimes=peakTimes;
    
end
end