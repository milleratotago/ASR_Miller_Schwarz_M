function thislike = asr_error1or2sigmas(parm,rtscon,rtsinc,SOA,TwoSigmas)
    % Compute likelihood of congruent & incongruent RTs as f(parameters in parm)
    %  for 4 different ASR models depending on:
    %    estimate_excitation or not (indicated by the number of parameters in parm)
    %    whether or not SigmaC changes when there is excitation/inhibition associated with lambda (indicated by TwoSigmas)
    % TwoSigmas is a Boolean indicating whether to use 1 or 2 sigmas:
    % TwoSigmas FALSE version uses a single SigmaC.
    %   parm can have 5 or 6 parameters:
    %   5: TauA, TauB, MuC, SigmaC, LambdaInh,    (assume lambdaExc=0)
    %   6: TauA, TauB, MuC, SigmaC, LambdaExc, LambdaInh
    % TwoSigmas TRUE version uses two SigmaC's for each condition.
    %   parm can have 6 or 8 parameters:
    %   6: TauA, TauB, MuC, SigmaC, LambdaInh, SigmaCinh (assume lambdaExc=0 and SigmaCexc=SigmaC)
    %   8: TauA, TauB, MuC, SigmaC, LambdaExc, SigmaCexc, LambdaInh, SigmaCinh
    %
    % If LambdaExc is used, negative values mean that there is facilitation (speeding) of RT
    % Note using abs() of TauA, TauB, MuC, SigmaC's but not Lambda's
    
    TauA = abs(parm(1));
    TauB = abs(parm(2));
    MuC = abs(parm(3));
    SigmaC = abs(parm(4));
    
    if TwoSigmas
        if numel(parm) == 6
            LambdaExc = 0;
            SigmaCexc = SigmaC;
            LambdaInh = parm(5);
            SigmaCinh = abs(parm(6));
        elseif numel(parm) == 8
            LambdaExc = parm(5);
            SigmaCexc = abs(parm(6));
            LambdaInh = parm(7);
            SigmaCinh = abs(parm(8));
        else
            error('Wrong number of parameters: asr_error1or2sigmas requires exactly 6 or 8 parameters when TwoSigmas is true.');
        end
    else
        % One Sigma
        SigmaCexc = SigmaC;
        SigmaCinh = SigmaC;
        if numel(parm) == 6
            LambdaExc = parm(5);
            LambdaInh = parm(6);
        elseif numel(parm) == 5
            LambdaExc = 0;
            LambdaInh = parm(5);
        else
            error('Wrong number of parameters: asr_error1or2sigmas requires exactly 5 or 6 parameters when TwoSigmas is false.');
        end
    end
    
    Alpha = 1/TauA;
    Beta = 1/TauB;
    pdfscon = f_rt(rtscon,Alpha,Beta,MuC,SigmaC,LambdaExc,SigmaCexc,SOA);
    pdfsinc = f_rt(rtsinc,Alpha,Beta,MuC,SigmaC,LambdaInh,SigmaCinh,SOA);
    thislike = -sum(log(pdfscon)) - sum(log(pdfsinc));
end

