function thislike = corerreg_neglnlikelihood(parm,corrtscon,corrtsinc,errrtscon,errrtsinc,varargin)
    % Compute likelihood of congruent/incongruent parm correct/error RTs as f(parameters in parm).
    % This version uses a single SigmaC.
    % Note: Use abs() of TauA, TauB, MuC, SigmaC but not Lambda's, and
    %    transform PCsup, PCexc, PCinh into 0-1 scale.
    % Note: parm can have 8 or 9 parameters:
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
    PCsup = maprto01(parm(RTparms+1));   % Pr(cor|A wins), so there is no exc or inhib
    PCexc = maprto01(parm(RTparms+2));   % Pr(cor|B wins) on a congruent trial
    PCinh = maprto01(parm(RTparms+3));   % Pr(cor|B wins) wins on an incongruent trial
    
%     Nconcor = numel(corrtscon);
%     Ninccor = numel(corrtsinc);
%     Nconerr = numel(errrtscon);
%     Nincerr = numel(errrtsinc);
    
    A = ExponenMn(TauA);
    if SOA==0
        B = ExponenMn(TauB);
    else
        B = AddTrans(ExponenMn(TauB),SOA);
    end
    PrBwins = PrXGTY(A,B);
    PrAwins = 1 - PrBwins;
    BifWinner = ConditXLTY(B,A);
    BifLoser  = ConditXGTY(B,A);
    
    C_Awins = Normal(MuC,SigmaC);             % Post-perceptual controlled processing time with no effect of A
    CwithExc = Normal(MuC+LambdaExc,SigmaC);  % Post-perceptual controlled processing time with excitation from A
    CwithInh = Normal(MuC+LambdaInh,SigmaC);  % Post-perceptual controlled processing time with inhibition from A
    
    RT_Awins = Convolution(C_Awins,BifLoser);     % Faster reversing order of C and B
    RTwithExc = Convolution(CwithExc,BifWinner);  % Faster reversing order of C and B
    RTwithInh = Convolution(CwithInh,BifWinner);  % Faster reversing order of C and B
    
    pdf_Awinscorrtscon = RT_Awins.PDF(corrtscon);
    pdfwithcorrtscon = RTwithExc.PDF(corrtscon);
    corpdfscon = pdf_Awinscorrtscon .* PrAwins * PCsup + pdfwithcorrtscon .* PrBwins * PCexc;
    
    pdf_Awinserrrtscon = RT_Awins.PDF(errrtscon);
    pdfwitherrrtscon = RTwithExc.PDF(errrtscon);
    errpdfscon = pdf_Awinserrrtscon .* PrAwins * (1 - PCsup) + pdfwitherrrtscon .* PrBwins * (1 - PCexc);
    
    pdf_Awinscorrtsinc = RT_Awins.PDF(corrtsinc);
    pdfwithcorrtsinc = RTwithInh.PDF(corrtsinc);
    corpdfsinc = pdf_Awinscorrtsinc .* PrAwins * PCsup + pdfwithcorrtsinc .* PrBwins * PCinh;
    
    pdf_Awinserrrtsinc = RT_Awins.PDF(errrtsinc);
    pdfwitherrrtsinc = RTwithInh.PDF(errrtsinc);
    errpdfsinc = pdf_Awinserrrtsinc .* PrAwins * (1 - PCsup) + pdfwitherrrtsinc .* PrBwins * (1 - PCinh);
    
    thislike = -(sum(log(corpdfscon)) + sum(log(corpdfsinc)) + sum(log(errpdfscon)) + sum(log(errpdfsinc)));
    
%     PCcon = PCsup*PrAwins + PCexc*PrBwins;
%     PCinc = PCsup*PrAwins + PCinh*PrBwins;
% 
%     thislike = -thislike - Nconcor*log(PCcon) - Nconerr*log(1-PCcon) - Ninccor*log(PCinc) - Nincerr*log(1-PCinc);
%     [thislike TauA TauB MuC SigmaC LambdaInh PCsup PCexc PCinh]
end

