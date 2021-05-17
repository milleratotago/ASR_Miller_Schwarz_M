function [rts, errcodes] = simRTs(nRTs, TauA,TauB,MuC,SigmaC,SOA, varargin)
    % Simulate nRTs and possibly correct/error booleans for a given set of model parameters.
    % TauA and TauB are the mean finishing times of the exponential processes
    % MuC and SigmaC are 2-position vectors with position 1 corresponding
    %  to the situation where A finishes first (excitation or inhibition is suppressed)
    %  and position 2 corresponding to the situation where B finishes first
    %  so there is excitation or inhibition.
    % varargin is PC, which is a 2-position vector indicating the probability
    %  of a correct response when A or B finishes first.
    
    A = exprnd(TauA,nRTs,1);
    B = exprnd(TauB,nRTs,1) + SOA;
    Awins = A < B;
    nAwins = sum(Awins);

    C = zeros(nRTs,1);
    C(Awins) = randn(nAwins,1) * SigmaC(1) + MuC(1);
    C(~Awins) = randn(nRTs-nAwins,1) * SigmaC(2) + MuC(2);
    rts = B + C - SOA;
    
    if numel(varargin) > 0
        PCs = varargin{1};
        errcodes = false(nRTs,1);
        unif = rand(nRTs,1);
        errcodes(Awins) = unif(Awins) > PCs(1);
        errcodes(~Awins) = unif(~Awins) > PCs(2);
    end
        
end
