function [phat,relativeError] = defineThresholdPrice(params,nt)
%  defineThresholdPrice calculates threshold price using linear
%  inter/extrapolation and guesses

% nt - number of time steps in the procedure
threshold = params.threshold;
X_0 = params.X_0;

phat1 = X_0*(1-threshold);
params.phat = phat1;
tic
[S,~,~,~,zGridArray,~,~] = solveViaMarkovChain1(params,nt);
toc
S01 = S(zGridArray(:,1,1) == 0,1,1);
clear S zGridArray

phat2 = 0.99*phat1;
params.phat = phat2;
tic
[S,~,~,~,zGridArray,~,~] = solveViaMarkovChain1(params,nt);
toc
S02 = S(zGridArray(:,1,1) == 0,1,1);
clear S zGridArray

phat = (phat2*S01-phat1*S02)/((phat2-phat1)/(1-threshold)+(S01-S02));
params.phat = phat;
tic
[S,~,~,~,zGridArray,~,~] = solveViaMarkovChain1(params,nt);
toc
S03 = S(zGridArray(:,1,1) == 0,1,1);
relativeError = ((S03-phat)/S03-threshold)/threshold;
end

