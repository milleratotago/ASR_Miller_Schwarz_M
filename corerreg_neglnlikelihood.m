function thislike = corerreg_neglnlikelihood(parm,corrtscon,corrtsinc,errrtscon,errrtsinc,varargin)
    % Compute likelihood of congruent/incongruent parm correct/error RTs as f(parameters in parm).
    % This version uses a single SigmaC.
    % Note using abs() of TauA, TauB, MuC, SigmaC but not Lambda's, and
    %    transform PCsup, PCexc, PCinh into 0-1 scale.
    % Note parm can have 8 or 9 parameters:
    %   if 8: TauA, TauB, MuC, SigmaC, LambdaInh, PCsup, PCexc, PCinh  (assume lambdaExc=0)
    %   if 9: TauA, TauB, MuC, SigmaC, LambdaExc, LambdaInh, PCsup, PCexc, PCinh
    % If LambdaExc is used, negative values mean that there is facilitation (speeding) of RT
    % Derived from errorDPs.m
    %
    % Equations: Likelihood is computed on the basis of relationships like this
    % for Cor and Err:
    %  Pr(RT=t & Cor) = Pr(RT=t & Awins & Cor) + Pr(RT=t & Bwins & Cor)
    %    = Pr(Awins) * Pr(Cor|Awins) * Pr(RT=t|Awins)
    %    + Pr(Bwins) * Pr(Cor|Bwins) * Pr(RT=t|Bwins)
    if numel(varargin)==0
        SOA = 0;
    else
        SOA = varargin{1};
    end
    TauA = abs(parm(1));
    TauB = abs(parm(2));
    MuC = abs(parm(3));
    SigmaC = abs(parm(4));
    if numel(parm) == 9
        LambdaExc = parm(5);
        LambdaInh = parm(6);
        RTparms = 6;
    else
        LambdaExc = 0;
        LambdaInh = parm(5);
        RTparms = 5;
    end
    PCsup = maprto01(parm(RTparms+1));  % A wins so there is no exc or inhib
    PCexc = maprto01(parm(RTparms+2));      % B wins on a congruent trial
    PCinh = maprto01(parm(RTparms+3));    % B wins on an incongruent trial
    
    A = ExponenMn(TauA);
    if SOA==0
        B = ExponenMn(TauB);
    else
        B = AddTrans(ExponenMn(TauB),SOA);
    end
    PrExc = PrXGTY(A,B);  % This is Pr(Bwins)
    BifWinner = ConditXLTY(B,A);
    BifLoser  = ConditXGTY(B,A);
    
    Cwithout = Normal(MuC,SigmaC);  % Post-perceptual controlled processing time with no effect of A
    CwithExc = Normal(MuC+LambdaExc,SigmaC);  % Post-perceptual controlled processing time with excitation from A
    CwithInh = Normal(MuC+LambdaInh,SigmaC);  % Post-perceptual controlled processing time with inhibition from A
    
    RTwithout = Convolution(Cwithout,BifLoser);  % Faster reversing order of C and B
    RTwithExc = Convolution(CwithExc,BifWinner);  % Faster reversing order of C and B
    RTwithInh = Convolution(CwithInh,BifWinner);  % Faster reversing order of C and B
    
    pdfwithoutcorrtscon = RTwithout.PDF(corrtscon);
    pdfwithcorrtscon = RTwithExc.PDF(corrtscon);
    corpdfscon = pdfwithoutcorrtscon .* (1-PrExc) * PCsup + pdfwithcorrtscon .* PrExc * PCexc;
    
    pdfwithouterrrtscon = RTwithout.PDF(errrtscon);
    pdfwitherrrtscon = RTwithExc.PDF(errrtscon);
    errpdfscon = pdfwithouterrrtscon .* (1-PrExc) * (1 - PCsup) + pdfwitherrrtscon .* PrExc * (1 - PCexc);
    
    pdfwithoutcorrtsinc = RTwithout.PDF(corrtsinc);
    pdfwithcorrtsinc = RTwithInh.PDF(corrtsinc);
    corpdfsinc = pdfwithoutcorrtsinc .* (1-PrExc) * PCsup + pdfwithcorrtsinc .* PrExc * PCinh;
    
    pdfwithouterrrtsinc = RTwithout.PDF(errrtsinc);
    pdfwitherrrtsinc = RTwithInh.PDF(errrtsinc);
    errpdfsinc = pdfwithouterrrtsinc .* (1-PrExc) * (1 - PCsup) + pdfwitherrrtsinc .* PrExc * (1 - PCinh);
    
    thislike = -sum(log(corpdfscon)) - sum(log(corpdfsinc)) - sum(log(errpdfscon)) - sum(log(errpdfsinc));
%     [thislike TauA TauB MuC SigmaC LambdaInh PCsup PCexc PCinh]
end

