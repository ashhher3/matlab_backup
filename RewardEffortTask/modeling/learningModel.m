function [alpha, beta, Xhat, iteration]=learningModel(C,V,X)

Xhat=X;
err=inf;
iteration=0;
progress=inf;

while (progress > 0.0001)
    if mod(iteration,2)==0
        [alpha, beta]=fitLearningRates(C,V,Xhat);
        Vhat=estimateValue(alpha,beta,V);
    else
        Vhat=estimateValue(alpha,beta,V);
        Xhat=logistic(Vhat,C);
    end
    pred=1./(1+exp(-Vhat*Xhat));
    newErr=evaluatePrediction(pred,C);
    progress=abs(err-newErr);
    err=newErr;
    iteration=1+iteration;
end

end
