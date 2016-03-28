function spiketrain = integrateAndFireSolution( inputs, threshold )
% Sample solution for Integrate and Fire problem at UCSF NS Bootcamp coding
% day 2015.  Takes neural data in a matrix with neurons as columns and
% times as rows and threshold as a scalar.  Returns a vector of 0s and 1s,
% where 1 indicates a spike, for each time step.
%
% KD & NH 7/23/2015

spiketrain = zeros(size(inputs,1),1); % initialize an empty spike train

for i = 1:size(inputs,1) % for all timesteps
    rawInput = sum(inputs(i,:)); % sum all the inputs at that timestep
    if rawInput > threshold % if we've passed the threshold
        spiketrain(i) = 1; % then we should spike
    end
end

% matlab will return spiketrain without us telling it to

end

