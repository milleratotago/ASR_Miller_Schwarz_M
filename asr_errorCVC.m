function thislike = asr_errorCVC(x,rtscon,rtsinc,varargin)
    % Compute likelihood of congruent & incongruent RTs as f(parameters in x).
    % This version uses a single stage C _coefficient of variance parameter CVC.
    % Note using abs() of TauA, TauB, MuC, CVC but not Lambda's
    % Note x can have 5 or 6 parameters:
    %   if 5, assume lambdaExc=0;
    %   if 6, LambdaExc is x(5) and lambdaInh is x(6)
    % If LambdaExc is used, negative values mean that there is facilitation (speeding) of RT
    if numel(varargin)==0
        SOA = 0;
    else
        SOA = varargin{1};
    end
    TauA = abs(x(1));
    TauB = abs(x(2));
    MuC = abs(x(3));
    CVC = abs(x(4));
    if numel(x) == 6
        LambdaExc = x(5);
        LambdaInh = x(6);
    else
        LambdaExc = 0;
        LambdaInh = x(5);
    end
    SigmaC = MuC*CVC;
    SigmaCexc = (MuC+LambdaExc)*CVC;
    SigmaCinh = (MuC+LambdaInh)*CVC;
    Alpha = 1/TauA;
    Beta = 1/TauB;
    pdfscon = f_rt(rtscon,Alpha,Beta,MuC,SigmaC,LambdaExc,SigmaCexc,SOA);
    pdfsinc = f_rt(rtsinc,Alpha,Beta,MuC,SigmaC,LambdaInh,SigmaCinh,SOA);
    thislike = -sum(log(pdfscon)) - sum(log(pdfsinc));
end

