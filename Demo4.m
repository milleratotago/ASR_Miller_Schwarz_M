% This demo provides an example of using the Cupid probability distribution toolbox
% to fit observed RTs with versions of the ASR model having processing time
% distributions other than the exponential and normal.

% The arbitrary model that is used as an illustration in this demo is defined
% in the asf_error function at the end of the file.

% CAUTION: This process is quite slow.  Running this script takes several minutes.

% Generate some random data just to use as an example.
% You would have your own data.
nTrials = 100;
obsRTcon = ExGauMn(300,30,120).Random(nTrials,1);
obsRTinc = ExGauMn(350,30,100).Random(nTrials,1);

% Set MATLAB fminsearch options if you want to adjust them.
SearchOptions = optimset('TolX',0.1,'TolFun',0.05);

% Provide some starting guesses for parameter values.
% These parameters correspond to the model defined below.
startparams = [50, 60, 300, 40, 50];  % TauA, TauB, mu, sigma, lambdaInh
startparams = sqrt(startparams);  % See note on parameter translation

ftominimize = @(x) asf_error(x,obsRTcon,obsRTinc);

[EndingVals,fval,exitflag,output] = fminsearch(ftominimize,startparams,SearchOptions);
EndingVals = EndingVals.^2;  % See note on parameter translation

function err = asf_error(params,obsRTcon,obsRTinc)
    % A: RNGammaMn(2,TauA)
    % B: RNGammaMn(3,TauB)
    % C: RNGammaMS(mu,sigma) if A<B
    %   or RNGammaMS(mu+LambdaInh,sigma)  if A>B in incongruent trial
    % The 5 parameters are thus TauA, TauB, mu, sigma, LambdaInh

    TauA = params(1)^2;  % See note on parameter translation
    TauB = params(2)^2;
    mu = params(3)^2;
    sigma = params(4)^2;
    LambdaInh = params(5)^2;
    
    A = RNGammaMn(2,TauA);
    B = RNGammaMn(3,TauB);
    Cwithout = RNGammaMS(mu,sigma);
    CwithExc = Cwithout;
    CwithInh =  RNGammaMS(mu+LambdaInh,sigma);
    SOA = 0;
    
    % Get the predicted RT distributions associated with these A,B,C stage finishing times.
    [RTcon, RTinc] = ASRDistsArb(A,B,Cwithout,CwithExc,CwithInh,SOA);
    pdfCon = RTcon.PDF(obsRTcon);
    pdfInc = RTinc.PDF(obsRTinc);
    lnlik = sum(log(pdfCon)) + sum(log(pdfInc));
    err = - lnlik;
end

% Note on parameter translation:
% ------------------------------
% fminsearch can potentially try both positive and negative numbers as parameter values,
% but we only want positive values for the parameters of this model.
% Thus, to make sure that our model is only given positive values to try,
% asf_error squares all of the suggested parameter values that it gets from fminsearch.
% In view of that, we should start with the square roots of the starting
% parameter values we really want, and at the end we must also square the
% final parameter values that fminsearch provides.