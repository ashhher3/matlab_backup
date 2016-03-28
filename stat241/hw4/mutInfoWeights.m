function weights = mutInfoWeights(adj, pTilde)
numV=size(adj,1);
numS=size(pTilde,2);

weights = zeros(size(adj));
for i=1:numV
    for j=i+1:numV
        if adj(i,j)==1
            joint = pairFromFullDistrib(pTilde,numS,i,j);
            pi=singleFromFullDistrib(pTilde,numS,i);
            pj=singleFromFullDistrib(pTilde,numS,j);
            
            w00=joint(1,1)*(log(joint(1,1) /(pi(1)*pj(1)) ));
            w01=joint(1,2)*(log(joint(1,2) /(pi(1)*pj(2)) ));
            w10=joint(2,1)*(log(joint(2,1) /(pi(2)*pj(1)) ));
            w11=joint(2,2)*(log(joint(2,2) /(pi(2)*pj(1)) ));
            
            weights(i,j)=w00+w01+w10+w11;
            weights(j,i)=weights(i,j);
        end
    end
end


end

function marginal = singleFromFullDistrib(dist,numS,i)
marginal = [0 0];
for s=1:numS % add up the states
    state = [dist(s).state];
    marginal(state(i)+1)=dist(s).prob+marginal(state(i)+1);
end

end

function joint = pairFromFullDistrib(dist,numS,i,j)
if j<i % check that we're following the i,j convention to match psi
    error('xi must be the smaller index')
end

joint = zeros(2,2); % initalize empty matrix

for s=1:numS % add up the states
    state = [dist(s).state];
    joint(state(i)+1,state(j)+1)=dist(s).prob+joint(state(i)+1,state(j)+1);
end

end