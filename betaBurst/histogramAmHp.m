function lags = histogramAmHp( AmTimes, HpTimes, coh )

lags=[];

win=1.5;

for h=1:length(HpTimes)
    hTime=HpTimes(h);
    
    % check coherence at hTime (high>=0.455, low<0.455)
%     sample=hTime*220;
%     width=10*220;
%     bin=floor(sample/width)+1;
%     if coh(bin)<0.455
        % find all Am peaks within win
        pre=find(AmTimes>=hTime-win,1,'first');
        post=find(AmTimes<=hTime+win,1,'last');

        for a=pre:post 
            lags = [lags, AmTimes(a)-hTime];
        end
        
%     end
end

figure()
hist(lags)
    
end