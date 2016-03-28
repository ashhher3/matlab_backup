function d = Dkl(dist1,dist2,numS)
d=0;
for s=1:numS
    d=d+dist1(s).prob*log(dist1(s).prob/dist2(s).prob);
end

end