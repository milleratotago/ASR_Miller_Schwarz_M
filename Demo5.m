% Demo of fitting both RTs and errors by minimizing the function
% corerreg_neglnlikelihood(parm,corrtscon,corrtsinc,errrtscon,errrtsinc,varargin)
%
% Warning: This is quite slow!

%% *********** This section generates some example data: You would have your own data to fit.

% True parameter values:
TauA = 75;       % Mean time for Automatic exponential racer
TauB = 90;       % Mean time for Controlled exponential racer
MuC = 310;       % Mean time of Gaussian component of controlled exGaussian, ignoring Lambdas
SigmaC = 35;     % SD of Gaussian component of controlled exGaussian, independent of Lambdas.
LambdaExc =  0;  % RT change with excitation for congruent trials if Controlled racer wins the race.
LambdaInh = 60;  % RT change with inhibition for incongruent trials if Controlled racer wins the race.
PCsup = 0.95;    % Prob correct, congruent & incongruent trials, A wins so irrelevant info is suppressed
PCexc = 0.98;    % Prob correct, congruent trials, B wins so facilitation (i.e., irrelevant info is NOT suppressed)
PCinh = 0.80;    % Prob correct, incongruent trials, B wins so inhibition (i.e., irrelevant info is NOT suppressed)

SOA =  0;

% Generate some simulated RTs with the above parameters:
nRTs = 500;
[rtscon, errcodescon] = simRTs(nRTs, TauA,TauB,[MuC, MuC+LambdaExc],[SigmaC SigmaC],SOA, [PCsup, PCexc]);
[rtsinc, errcodesinc] = simRTs(nRTs, TauA,TauB,[MuC, MuC+LambdaInh],[SigmaC SigmaC],SOA, [PCsup, PCinh]);

% Select out the RTs of correct & incorrect trials:
corrtscon = rtscon(~errcodescon);
corrtsinc = rtsinc(~errcodesinc);
errrtscon = rtscon( errcodescon);
errrtsinc = rtsinc( errcodesinc);


%% Here is the parameter estimation:
errfun = @(x) corerreg_neglnlikelihood(x,corrtscon,corrtsinc,errrtscon,errrtsinc);
options = optimset('Display','iter');
[holdparms, holdbest] = fminsearch(errfun,parm,options);
