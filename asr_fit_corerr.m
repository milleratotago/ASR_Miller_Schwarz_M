function results = asr_fit_corerr(corcongruent_rts,corincongruent_rts,errcongruent_rts,errincongruent_rts,soa,varargin)
    % Function to fit ASR model to observed correct/error X congruent/incongruent RTs.
    % REQUIRED INPUTS:
    %   corcongruent_rts, corincongruent_rts, errcongruent_rts, errincongruent_rts: vectors of the correct & error rts from the two conditions
    %   soa: the time from the onset of the irrelevant attribute to the onset of the relevant, in msec.
    % OPTIONAL INPUTS (in any order following the required inputs):
    %   estimate_excitation: a boolean indicating whether the lambda_exc parameter is free to vary (default: false)
    %   sigmaProps: a vector used in determining the starting parameters for multiple runs of fminsearch.
    %      The numbers in this vector are proportions (eps,1-eps) of the overall sd of RT allocated to the normal component.
    %   optparms: a MATLAB structure, produced by optimset, passed to fminsearch to set its options.
    %   'TwoSigmas' : fit model with separate sigmas depending on whether A or B wins the race.
    %   ***NOTE: TwoSigmas OPTION NOT YET SUPPORTED.
    % OUTPUT:
    %   results: a struct with field of
    %     bestParms: the best parameters found across all searches
    %     bestNegLn: the lowest value of -log(likelihood) found across all searches
    %     sigmaProps: the vector of sigmaProps values used to choose starting parameters
    %     searchResults: best parameter estimates found in fminsearches with each of the different starting points
    %     searchBests: the lowest value of -log(likelihood) found in fminsearches with each of the different starting points
    
    % set defaults for optional parameters
    estimate_excitation = false;
    sigmaProps = [0.1 0.3 0.5 0.7 0.9];  % Proportion of grandsd RT in sigma
    useOptParms = false;
    TwoSigmas = false;
    
    % Check varargin to see if any optional parameters are over-ridden.
    % The optional parameters all have different types so we can discriminate them that way.
    for iParm=1:numel(varargin)
        thisParm = varargin{iParm};
        if (ischar(thisParm) && strcmpi(thisParm,'TwoSigmas'))
            TwoSigmas = true;
        elseif islogical(thisParm)
            estimate_excitation = thisParm;
        elseif isfloat(thisParm)
            sigmaProps = thisParm;
        elseif isstruct(thisParm) && isfield(thisParm,'InitTrustRegionRadius')
            % Also check some obscure optimset field that hopefully won't be present in other structures
            useOptParms = true;
            optparms = thisParm;
        else
            error('Unrecognized argument');
        end
    end
    
    % Choose how many parameters are to be estimated.
    if estimate_excitation && TwoSigmas
        nParms = 11;  % TauA, TauB, MuC, SigmaC, LambdaExc, SigmaCexc, LambdaInh, SigmaCinh, PCsup, PCexc, PCinh
    elseif ~estimate_excitation && TwoSigmas
        nParms = 9;  % TauA, TauB, MuC, SigmaC, LambdaInh, SigmaCinh, PCsup, PCexc, PCinh (assume lambdaExc=0 and SigmaCexc=SigmaC)
    elseif estimate_excitation && ~TwoSigmas
        nParms = 9;  % TauA, TauB, MuC, SigmaC, LambdaExc, LambdaInh, PCsup, PCexc, PCinh
    elseif ~estimate_excitation && ~TwoSigmas
        nParms = 8;  % TauA, TauB, MuC, SigmaC, LambdaInh, PCsup, PCexc, PCinh
    end
    
    
    % Initialize arrays to hold results of multiple fminsearch's
    nTries = numel(sigmaProps);  % Number of times to run fminsearch with different starting values.
    holdparms = zeros(nTries,nParms);
    holdbest = zeros(nTries,1);
    
    % Compute summary stats on observed RTs to use in setting initial parameters.
    congruent_rts = [corcongruent_rts(:); errcongruent_rts(:)];
    incongruent_rts = [corincongruent_rts(:); errincongruent_rts(:)];
    meanCon = mean(congruent_rts);
    meanInc = mean(incongruent_rts);
    grandmean = (meanCon + meanInc)/2;
    grandsd = std([congruent_rts; incongruent_rts]);
    
    % define the error function for fminsearch to minimize
    if TwoSigmas
        errfun = @(x) corerreg_neglnlikelihood2sigmas(x,corcongruent_rts,corincongruent_rts,errcongruent_rts,errincongruent_rts,soa);  % NOT YET SUPPORTED
    else
        errfun = @(x) corerreg_neglnlikelihood(x,corcongruent_rts,corincongruent_rts,errcongruent_rts,errincongruent_rts,soa);
    end

    % Get starting values for PC parameters (I am just guessing here):
    % My intuition is that incongruence reduces accuracy more than congruence increases it, relative to a control condition with no irrelevant information
    PCcongru = numel(corcongruent_rts) / ( numel(corcongruent_rts) + numel(errcongruent_rts) );
    PCincongru = numel(corincongruent_rts) / ( numel(corincongruent_rts) + numel(errincongruent_rts) );
    PCsup = 0.75*PCcongru + 0.25*PCincongru;
    PCexc = PCcongru;
    PCinh = PCincongru;
    rPCsup = map01tor(PCsup);
    rPCexc = map01tor(PCexc);
    rPCinh = map01tor(PCinh);
        
    % run fminsearch with different sets of starting parameters.
    for iTry=1:nTries
        % get some guesses for starting parameters:
        TauA = sigmaProps(iTry) * grandsd;
        TauB = TauA;
        MuC = grandmean - TauB;
        SigmaC = sqrt((grandsd^2) - (TauB^2));
        LambdaInh = 2 * (meanInc - meanCon);

        if estimate_excitation && TwoSigmas
            startparms = [TauA, TauB, MuC, SigmaC, 0, SigmaC, LambdaInh, SigmaC, rPCsup, rPCexc, rPCinh];  % Start lambdaExc at 0
        elseif ~estimate_excitation && TwoSigmas
            startparms = [TauA, TauB, MuC, SigmaC, LambdaInh, SigmaC, rPCsup, rPCexc, rPCinh]; % (assume lambdaExc=0 and SigmaCexc=SigmaC)
        elseif estimate_excitation && ~TwoSigmas
            startparms = [TauA,TauB,MuC,SigmaC,0,LambdaInh, rPCsup, rPCexc, rPCinh];  % Start lambdaExc at 0
        elseif ~estimate_excitation && ~TwoSigmas
            startparms = [TauA,TauB,MuC,SigmaC,LambdaInh, rPCsup, rPCexc, rPCinh];
        end
        
        if useOptParms
            [holdparms(iTry,:), holdbest(iTry)] = fminsearch(errfun,startparms,optparms);
        else
            [holdparms(iTry,:), holdbest(iTry)] = fminsearch(errfun,startparms);
        end
    end
    
    % eg_neglnlikelihood uses abs() of these parms:
    holdparms(:,1) = abs(holdparms(:,1));  % tauA
    holdparms(:,2) = abs(holdparms(:,2));  % tauB
    holdparms(:,3) = abs(holdparms(:,3));  % mu
    holdparms(:,4) = abs(holdparms(:,4));  % sigma
    if estimate_excitation && TwoSigmas
        holdparms(:,6) = abs(holdparms(:,6));  % sigmaCexc
        holdparms(:,8) = abs(holdparms(:,8));  % sigmaCinh
    elseif ~estimate_excitation && TwoSigmas
        holdparms(:,6) = abs(holdparms(:,6));  % sigmaCexc
    end
    % map the three real PC parms back to 0/1:
    holdparms(:,end-2) = maprto01(holdparms(:,end-2));
    holdparms(:,end-1) = maprto01(holdparms(:,end-1));
    holdparms(:,end) = maprto01(holdparms(:,end));
    
    % See which fminsearch gave the best overall results
    [minbest, minpos] = min(holdbest);
    
    results.bestParms = holdparms(minpos,:);
    results.bestNegLn = minbest;
    results.sigmaProps = sigmaProps;
    results.searchResults = holdparms;
    results.searchBests = holdbest;
    
end

