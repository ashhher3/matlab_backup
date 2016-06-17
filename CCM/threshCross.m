function crossings = threshCross(data,X,thresh)

testX=X(1):X(end);
crossings = nan(size(data,1),1);

for d = 1:size(data,1)
    f = fit(X', data(d,:)','power2');
    c=find(feval(f,testX')>=thresh,1,'first');
    if ~isempty(c)
        crossings(d,1)=c;
    end
end

end

