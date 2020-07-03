% This demo shows how to:
%
% o Create predicted RTcon and RTinc distributions for a given set of
%   distributions and parameter values.
% o Plot the predicted PDFs and CDFs.
% o Plot the predicted delta plots.
%
% There are separate examples with exponential A, B and Gaussian C,
%  and with arbitrary distributions for A, B, and C.

%% Example with exponential A, B and Gaussian C, as discussed in the paper.

% Specify the parameters:
TauA = 80;
TauB = 100;
MuC = 400;
SigmaC = 50;
LambdaExc =  0;
LambdaInh = 70;
SOA = 0;

% Create the distributions

[RTcon, RTinc, PrInhib] = ASRDistsEG(TauA,TauB,MuC,SigmaC,LambdaExc,LambdaInh,SOA);

figure;
cumProbs = [0.01 0.05 0.1:0.1:0.9 0.95 0.99];
ASRplots(RTcon,RTinc,cumProbs);


%% Example with arbitrary distributions of finishing times for A, B, and C.

% Specify the distributions and their parameters:
%
% For more information about what distributions can be used
% and what their parameters are, check the Cupid documentation.
%
% For simplicity, this example uses Gamma's for all 5 of these distributions
% just for simplicity, but they could be any set of 5 cupid distributions---
% i.e., all different if you want.


A = RNGamma(2,.01);
B = RNGamma(3,.01);
Cwithout = RNGamma(40,.1);
CwithExc =  RNGamma(40,.12);
CwithInh =  RNGamma(40,.08);
SOA = 0;

% Get the predicted RT distributions associated with these finishing times.
[RTcon2, RTinc2, PrInhib2] = ASRDistsArb(A,B,Cwithout,CwithExc,CwithInh,SOA);

figure;
cumProbs2 = [0.01 0.05 0.1:0.1:0.9 0.95 0.99];
ASRplots(RTcon2,RTinc2,cumProbs2);

