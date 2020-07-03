function [RTcon, RTinc, PrBwins] = ASRDistsEG(TauA,TauB,MuC,SigmaC,LambdaExc,LambdaInh,SOA)
    % Create predicted RTcon & RTinc distributions for ASR model with exponential
    %   finishing times for processes A & B and Gaussian finishing times for process C.
    % ALL TIMES MUST BE IN THE SAME UNITS (E.G., MILLISECONDS).
    %
    % INPUT PARAMETERS:
    %
    % TauA: mean finishing time for activation suppression process
    %
    % TauB: mean finishing time for relevant info extraction process
    %
    % MuC & SigmaC: mean & std deviation of stage C finishing times
    %   without influence from irrelevant activation
    %
    % LambdaExc: Change in mean stage C finishing times with excitation in
    %   congruent trials (i.e., irrelevant activation has not been suppressed).
    %   NEGATIVE NUMBER INDICATES SPEEDUP
    %
    % LambdaInh: Change in mean stage C finishing times with inhibition in
    %   incongruent trials (i.e., irrelevant activation has not been suppressed).
    %   POSITIVE NUMBER INDICATES SLOWDOWN
    %
    % SOA: time from the start of process A to the start of process B
    %   (0 if they start simultaneously).
    %
    % OUTPUTS:
    %
    % RTcon: Distribution of RTs in congruent trials.
    %
    % RTinc: Distribution of RTs in incongruent trials.
    %
    % PrBwins: Probability that B wins the A/B race so that stage C is affected by excitation or inhibition.
    
    % Time 0 is the onset of A at the beginning of this function.
    % Time 0 will be changed to the onset of B at the end.

    A = ExponenMn(TauA);
    if SOA == 0
        B = ExponenMn(TauB);
    else
        B = AddTrans(ExponenMn(TauB),SOA);
    end
    
    PrBwins = PrXGTY(A,B);
    BifWinner = ConditXLTY(B,A);
    BifLoser  = ConditXGTY(B,A);
    
    % Distribution of post-race processing time...
    Cwithout = Normal(MuC,SigmaC);            % ... with no effect of irrelevant attribute
    CwithExc = Normal(MuC+LambdaExc,SigmaC);  % ... with excitation from irrelevant attribute
    CwithInh = Normal(MuC+LambdaInh,SigmaC);  % ... with inhibition from irrelevant attribute
    
    % Distributions of RTs in various cases.
    % (The convolutions are faster with C first.)
    RTwithout = Convolution(Cwithout,BifLoser);   % RT with no effect of irrelevant attribute
    RTwithExc = Convolution(CwithExc,BifWinner);  % RT with excitation from irrelevant attribute
    RTwithInh = Convolution(CwithInh,BifWinner);  % RT with inhibition from irrelevant attribute
    
    RTcon = Mixture(PrBwins,RTwithExc,RTwithout);
    RTinc = Mixture(PrBwins,RTwithInh,RTwithout);
    
    % Now adjust Time 0 to the onset of B if SOA was not zero.
    if SOA ~= 0
        RTcon = AddTrans(RTcon,-SOA);
        RTinc = AddTrans(RTinc,-SOA);
    end

    RTcon.StringName = 'RTcon';
    RTinc.StringName = 'RTinc';
    
end
