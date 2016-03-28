%% all times Am above threshold given Hp above threshold

% start 0.5 seconds in so we can look 0.5 sec in each direction
step = 1.5*Fs;
lagsHA=[];

for i=step+1:length(betaTime)-step
    if hp.contact(1).burst(i)==1 % there's a Hp burst
        for j=i-step:i+step % check 0.5s on either side
            if am.contact(1).burst(j)==1 % if there's an Am burst
                lagsHA = [lagsHA, betaTime(j)-betaTime(i)]; % record lag
            end
        end
    end
end

%% all Am bursts within +-1.5s of each Hp burst

allHA=[];

for i=1:length(hp.contact(1).startTimes)
    h=hp.contact(1).startTimes(i);
    
    start=find([am.contact(1).startTimes]>=h-1.5,1,'first');
    stop=find([am.contact(1).startTimes]<=h+1.5,1,'last');
    
    
    for j=start:stop
        allHA = [allHA, am.contact(1).startTimes(j)-h ];
    end
    
end

%% closest Am burst to each Hp burst

closeHA=[];

for i=1:length(hp.contact(1).startTimes)
    h=hp.contact(1).startTimes(i);
    
    % find closest preceeding burst
    b=am.contact(1).startTimes(find([am.contact(1).startTimes]<h,1,'last'));
    
    % find closest following burst
    a=am.contact(1).startTimes(find([am.contact(1).startTimes]>h,1,'first'));
    
    if abs(b-h)<(a-h)
        a=b;
    end
    
    closeHA = [closeHA, a-h];
end


%% closest Hp burst to each Am burst

closeAH=[];

for i=1:length(am.contact(1).startTimes)
    a=am.contact(1).startTimes(i);
    
    % find closest preceeding burst
    b=hp.contact(1).startTimes(find([hp.contact(1).startTimes]<a,1,'last'));
    
    % find closest following burst
    h=hp.contact(1).startTimes(find([hp.contact(1).startTimes]>a,1,'first'));
    
    if abs(b-a)<(h-a)
        h=b;
    end
    
    closeAH = [closeAH, h-a];
end

