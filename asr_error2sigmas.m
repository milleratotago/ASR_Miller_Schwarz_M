function thislike = asr_error2sigmas(x,rtscon,rtsinc,varargin)
    % Compute likelihood of congruent & incongruent RTs as f(parameters in x).
    % This version uses two SigmaC's for each condition.
    % Note x can have 6 or 8 parameters:
    %   6 : TauA, TauB, MuC, SigmaC, LambdaInh, SigmaCinh (assume lambdaExc=0 and SigmaCexc=SigmaC)
    %   8 : TauA, TauB, MuC, SigmaC, LambdaExc, SigmaCexc, LambdaInh, SigmaCinh
    % If LambdaExc is used, negative values mean that there is facilitation (speeding) of RT
    % Note using abs() of TauA, TauB, MuC, SigmaC's but not Lambda's
    if numel(varargin)==0
        SOA = 0;
    else
        SOA = varargin{1};
    end
    TauA = abs(x(1));
    TauB = abs(x(2));
    MuC = abs(x(3));
    SigmaC = abs(x(4));
    if numel(x) == 6
        LambdaExc = 0;
        SigmaCexc = SigmaC;
        LambdaInh = x(5);
        SigmaCinh = abs(x(6));
    elseif numel(x) == 8
        LambdaExc = x(5);
        SigmaCexc = abs(x(6));
        LambdaInh = x(7);
        SigmaCinh = abs(x(8));
    else
        error('Wrong number of parameters');
    end
    Alpha = 1/TauA;
    Beta = 1/TauB;
    pdfscon = f_rt(rtscon,Alpha,Beta,MuC,SigmaC,LambdaExc,SigmaCexc,SOA);
    pdfsinc = f_rt(rtsinc,Alpha,Beta,MuC,SigmaC,LambdaInh,SigmaCinh,SOA);
    thislike = -sum(log(pdfscon)) - sum(log(pdfsinc));
end

