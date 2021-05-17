% This demo provides the simplest (and fastest) example of fitting observed RTs
% with the exG version of the ASR model.

SOA = 0;  % Determined by the experimental design

% Generate some random data just to use as an example.
% You would have your own data.
nTrials = 1000;
trueTauA = 100;
trueTauB = 140;
trueMu = 400;
trueSigma = 40;
trueLambdaExc = -50;
trueLambdaInh = 80;
obsRTcon = eg_random(trueTauA,trueTauB,trueMu,trueSigma,trueLambdaExc,SOA,nTrials);
obsRTinc = eg_random(trueTauA,trueTauB,trueMu,trueSigma,trueLambdaInh,SOA,nTrials);


% Now fit four versions of the model:

% one sigma, lambda_Excitation is NOT allowed to vary:
results_no_excitation = asr_fit(obsRTcon,obsRTinc,SOA);

% one sigma, lambda_Excitation is allowed to vary:
results_with_excitation = asr_fit(obsRTcon,obsRTinc,SOA,true);

% two sigmas, lambda_Excitation is NOT allowed to vary:
results_no_excitation_2sig = asr_fit(obsRTcon,obsRTinc,SOA,'TwoSigmas');

% two sigmas, lambda_Excitation is allowed to vary:
results_with_excitation_2sig = asr_fit(obsRTcon,obsRTinc,SOA,true,'TwoSigmas');

% increase the number of function evals
optparms = optimset('MaxFunEvals',2000);

% lambda_Excitation is allowed to vary & use user-defined optimset parameters:
results_with_excitation2 = asr_fit(obsRTcon,obsRTinc,SOA,true,optparms);

% Use different starting points for fminsearch:
mySigmaProps = [0.1 0.5 0.9];
results_no_excitation2 = asr_fit(obsRTcon,obsRTinc,SOA,mySigmaProps);


