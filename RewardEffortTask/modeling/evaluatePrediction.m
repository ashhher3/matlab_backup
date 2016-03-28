function err = evaluatePrediction(pred,choice)
if length(pred) ~= length(choice)
    error('prediction and choice are different lenghts')
end

n=length(pred);
err=0;
for i=1:n
    err=err+abs(choice(i)-pred(i));
end
err=err/n;

end