function [RTcon, RTinc, PrBwins] = ASRDistsArb(A,B,Cwithout,CwithExc,CwithInh,SOA)
    % Create predicted RTcon & RTinc distributions for ASR model with arbitrary
    %   distributions of process finishing times.
    %
    % ALL TIMES MUST BE IN THE SAME UNITS (E.G., MILLISECONDS).
    %
    % INPUT PARAMETERS:
    %
    % A: distribution of finishing times for activation suppression process
    %
    % B: distribution of finishing times for relevant info extraction process
    %
    % Cwithout: distribution of stage C finishing times without influence from
    %   irrelevant activation (i.e., irrelevant activation has been suppressed).
    %
    % CwithExc: distribution of stage C finishing times with excitation in
    %   congruent trials (i.e., irrelevant activation has not been suppressed).
    %
    % CwithInh: distribution of stage C finishing times with inhibition in
    %   incongruent trials (i.e., irrelevant activation has not been suppressed).
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

    if SOA == 0
        B2 = B;
    else
        B2 = AddTrans(B,SOA);
    end
    
    PrBwins = PrXGTY(A,B2);
    BifWinner = ConditXLTY(B2,A);
    BifLoser  = ConditXGTY(B2,A);
    
    RTwithout = Convolution(Cwithout,BifLoser);   % Faster reversing order of C and B
    RTwithExc = Convolution(CwithExc,BifWinner);  % Faster reversing order of C and B
    RTwithInh = Convolution(CwithInh,BifWinner);  % Faster reversing order of C and B
    
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
