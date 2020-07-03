function rts = eg_random(TauA,TauB,MuC,SigmaC,Lambda,SOA,nTrials)
    % Generate random values from the model.
    A = exprnd(TauA,nTrials,1);
    B = exprnd(TauB,nTrials,1) + SOA;
    C = normrnd(MuC,SigmaC,nTrials,1);
    Bwins = B < A;
    C(Bwins) = C(Bwins) + Lambda;
    rts = B + C;
end
